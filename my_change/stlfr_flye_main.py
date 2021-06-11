# -*- coding: utf-8 -*-
# @Time : 2021/6/3 11:02
# @Author : Jiangzhesheng
# @File : stlfr_flye_main.py
# @Software: PyCharm
import sys
import argparse
import logging

from my_change import subgraph
from flye.repeat_graph.repeat_graph import *
from my_change import solve_repeat

BWA="/share/app/bwa/0.7.17/"
SAMTOOLS="/share/app/samtools/1.11/bin/"

def chech_flye(flye_dir: str,logger):
    """
    用于检查flye的中间结果完整性，返回repeat_graph和edges_fasta这两个关键文件路径
    :param flye_dir:
    :param logger:
    :return:repeat_graph和edges_fasta路径
    """
    edge_file=os.path.join(flye_dir,'repeat_graph_edges.fasta')
    repeat_graph=os.path.join(flye_dir,'repeat_graph_dump')
    if not (os.path.isfile(edge_file) and os.path.isfile(repeat_graph)):
        logger.error("flye_dir is incomplete")
        sys.exit(1)
    return edge_file,repeat_graph

def aln_stlfr_flye(read1, read2, ref,
                   thread, outputdir, logger):
    """
    比对模块，将stLFR_read比对到ref上
    过滤掉了多比对结果（-F 2048）
    :param read1:read1路径
    :param read2:read2路径
    :param ref:ref路径
    :param thread:
    :param outputdir:
    :param logger:
    :return:最终bam文件路径
    """
    thread=str(thread)
    outputdir=os.path.join(outputdir,"bwa_result")
    os.mkdir(outputdir)
    os.symlink(src=ref, dst=os.path.join(outputdir, os.path.basename(ref)))
    ref=os.path.join(outputdir, os.path.basename(ref))

    bwa_log=os.path.join(outputdir,'bwa_log')
    bwa_log = open(bwa_log, mode='ab+', buffering=0)
    #1、bwa建索引
    logger.info('start bwa index')
    cmd=['bwa','index', ref]
    logger.debug(' '.join(cmd))
    pp=subprocess.Popen(cmd,stdout=bwa_log,stderr=bwa_log)
    pp.wait()
    assert pp.returncode==0,'bwa index failing'
    #2、bwa比对
    logger.info('start bwa mem')
    samfile=os.path.join(outputdir, 'stlfr_flye.sam')
    cmd=['bwa','mem',
         '-t', thread,
         '-o', samfile,
         ref, read1, read2]
    logger.debug(' '.join(cmd))
    pp=subprocess.Popen(cmd,stdout=bwa_log,stderr=bwa_log)
    pp.wait()
    assert pp.returncode==0,'bwa mem failing'


    samtools_log=os.path.join(outputdir,'samtools_log')
    samtools_log = open(samtools_log, mode='ab+', buffering=0)
    # 3、转成bam，过滤掉多比对
    logger.info('start samtools view')
    bamfile=os.path.join(outputdir, 'stlfr_flye.bam')
    cmd=['samtools','view',
         '-bS','-@',thread,'-F','2048',
         '-o',bamfile,samfile]
    logger.debug(' '.join(cmd))
    pp=subprocess.Popen(cmd,stdout=samtools_log,stderr=samtools_log)
    pp.wait()
    assert pp.returncode==0,'samtools view failing'

    logger.info('start samtools sort')
    sortbamfile = os.path.join(outputdir, 'stlfr_flye.sort.bam')
    cmd=['samtools','sort','-@',thread,'-o',sortbamfile,bamfile]
    logger.debug(' '.join(cmd))
    pp=subprocess.Popen(cmd,stdout=samtools_log,stderr=samtools_log)
    pp.wait()
    assert pp.returncode==0,'samtools sort failing'

    logger.info('start samtools index')
    cmd=['samtools','index','-@',thread,sortbamfile]
    logger.debug(' '.join(cmd))
    pp=subprocess.Popen(cmd,stdout=samtools_log,stderr=samtools_log)
    pp.wait()
    assert pp.returncode==0,'samtools index failing'

    return sortbamfile

def show_subgs(subgs:List[RepeatGraph],outputdir:str)->None:
    """
    可选流程，展示所有子图
    :param subgs:
    :param outputdir:
    :return:None
    """
    outputdir=os.path.join(outputdir,'subgs')
    os.mkdir(outputdir)
    for i in range(len(subgs)):
        subgs[i].output_svg(prefixname=os.path.join(outputdir,'subg%d'%i))


def stlfr_flye_main(read1:str, read2:str, flye_dir:str,
                    outputdir:str, thread:str,logger,bam)->None:
    """
    大致流程框架，
    1、首先检查flye完整性
    2、如果没有了已经index和sort后的bam则进行比对模块
    3、读入flye解完的repeat_graph
    4、拆成子图
    5、对于每个子图进行解repeat
    :param read1:路径
    :param read2:路径
    :param flye_dir:20-repeat文件夹路径
    :param outputdir:输出文件夹路径
    :param thread:
    :param logger:
    :param bam:bam文件路径
    :return:
    """

    edge_file,graph=chech_flye(flye_dir,logger)
    if not bam:
        os.environ["PATH"]=os.environ["PATH"]+os.pathsep+BWA+os.pathsep+SAMTOOLS
        bam=aln_stlfr_flye(read1=read1, read2=read2, outputdir=outputdir,
                           ref=edge_file,
                           thread=thread, logger=logger)
        logger.info('align finish')

    repeat_graph = RepeatGraph(fp.read_sequence_dict(edge_file))
    repeat_graph.load_from_file(graph)

    subgs=subgraph.graph2subgraph(repeat_graph=repeat_graph)
    # show_subgs(subgs=subgs,outputdir=outputdir)

    for i in range(0,len(subgs),2):#间隔为2是因为两个子图结构是一模一样的，只是edge_id取了反，解其中一个另一个自动解决
        logger.info("processing subg%d"%i)
        subgdir=os.path.join(outputdir,'subg%d'%i)
        os.mkdir(subgdir)
        subg=subgs[i]
        repeat_graph=solve_repeat.solve_repeat_main(subgraph=subg, bamfile=bam, origraph=repeat_graph,
                                                    outputdir=subgdir,logger=logger)

    repeat_graph.output_all(outdir=outputdir)

def main(argv):
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-r1', '--stlfr_read1',required=True)
    parser.add_argument('-r2','--stlfr_read2',required=True)
    parser.add_argument('-flye','--flye_20_repeat',required=True)
    #该参数影响流程
    parser.add_argument('-bam','--sorted_index_bamfile',default='')

    parser.add_argument('-o', '--outputdir',default='stlfr_flye_result')
    parser.add_argument('-t','--thread',default=8)

    args = parser.parse_args(argv[1:])

    read1 = args.stlfr_read1
    read2 = args.stlfr_read2
    flye_dir=args.flye_20_repeat
    outputdir = args.outputdir
    thread= args.thread
    bam=args.sorted_index_bamfile

    try:
        os.mkdir(outputdir)
    except FileExistsError as e:
        shutil.rmtree(outputdir)
        os.mkdir(outputdir)
        print(outputdir, ' is exist, files in it may be overwritten')

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    sh = logging.StreamHandler()
    shfmt = logging.Formatter('%(asctime)s-%(message)s')
    sh.setFormatter(fmt=shfmt)
    sh.setLevel(logging.INFO)
    logger.addHandler(sh)

    fh = logging.FileHandler(filename=os.path.join(outputdir, 'totallog'), mode='w')
    fhfmt = logging.Formatter('%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s')
    fh.setFormatter(fmt=fhfmt)
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)

    logger.debug(args)

    stlfr_flye_main(read1=os.path.abspath(read1),
                    read2=os.path.abspath(read2),
                    flye_dir=os.path.abspath(flye_dir),
                    outputdir=os.path.abspath(outputdir),
                    logger=logger,
                    thread=thread,
                    bam=bam)

if __name__ == '__main__':
    main(sys.argv)