# -*- coding: utf-8 -*-
# @Time : 2021/4/28 13:19
# @Author : Jiangzhesheng
# @File : bwa_stLFR_flye.py
# @Software: PyCharm
import sys
import os
import subprocess
from typing import List, Dict, Set
from bioclass.sam import GeneralSamFile, Sam
from flye.repeat_graph.repeat_graph import RepeatGraph, RgEdge, RgNode, EdgeSequence
import flye.utils.fasta_parser as fp
import my_change.subgraph as subgraph
import shutil

def bamfilter(bam:str, regions:List[str],outputdir:str):
    """
    从bam文件中抽提出比对到ref：start-end的read信息
    :param bam:需要过滤的bam文件路径
    :param regions:提取区域列表
    例如：['edge1:0-10000','edge2']
    :return: 过滤后的Sam文件路径
    """
    outfile=os.path.join(outputdir,'tmp.sam')
    assert 'sort' in bam,'bam file must be sort'
    assert os.path.isfile(bam+'.bai'),'index file must be exist'

    with open(os.path.join(outputdir,'regions'),mode='w') as r:
        for i in regions:
            r.write(i+'\n')

    cmd=['/share/app/samtools/1.11/bin/samtools','view','-F','2048','-@','8','-o',outfile,bam]
    cmd=cmd+regions
    p=subprocess.Popen(cmd)
    p.wait()
    return outfile

def subgbam(subg:RepeatGraph,bampath:str,outputdir:str,flanking:int=-1):
    """
    对于一个子图subg，构造需要得到的read的regions序列，输出sam文件
    如果是repeat边，所有比对上去的read都要
    如果是unique边，只需要取 和repeat相接的一段序列区域(默认全取)
    :param subg:子图
    :param bampath:所有read比对到graph的总bam文件路径
    :param flanking:unique边取的长度，如果为-1则全取，默认全取
    :return:sam文件路径
    """

    regions=[]

    def edge_region2str(edge:RgEdge,start:int=-1,end:int=-1)->List[str]:
        ans=[]
        for seq in edge.edge_sequences:
            chr_name = seq.edge_seq_name[1:]
            region=chr_name+":%d-%d" % (start, end) if start>=0 and end>start else chr_name
            ans.append(region)
        return ans

    for edge in subg.edges.values():
        if edge.repetitive==True:
            regions+=edge_region2str(edge=edge)
            for i in edge.node_left.in_edges:
                if i.repetitive==False:
                    if i.edge_id<0:
                        start = 0 if flanking != -1 else -1
                        end = flanking if flanking != -1 else -1
                        regions += edge_region2str(edge=i, start=start, end=end)
                    else:
                        start=i.length()-flanking if flanking!=-1 else -1
                        end=i.length() if flanking!=-1 else -1
                        regions+=edge_region2str(edge=i,start=start,end=end)
            for i in edge.node_right.out_edges:
                if i.repetitive == False:
                    if i.edge_id<0:
                        start = i.length() - flanking if flanking != -1 else -1
                        end = i.length() if flanking != -1 else -1
                        regions += edge_region2str(edge=i, start=start, end=end)
                    else:
                        start=0 if flanking!=-1 else -1
                        end=flanking if flanking!=-1 else -1
                        regions+=edge_region2str(edge=i,start=start,end=end)
    sam=bamfilter(bam=bampath,regions=regions,outputdir=outputdir)
    return sam

    print(d)
    return d

def main(argv):
    workdir = '/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/step_by_step'
    bamfile='/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/bwa_stLFR_flye/chr19/chr19.sort.bam'
    edge_file = os.path.join(workdir, '20-repeat-copy/repeat_graph_edges.fasta')
    graph = os.path.join(workdir, '20-repeat-copy/repeat_graph_dump_after_rr')
    repeat_graph = RepeatGraph(fp.read_sequence_dict(edge_file))
    repeat_graph.load_from_file(graph)

    outputdir="/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/bwa_stLFR_flye/tmp"

    try:
        os.mkdir(outputdir)
    except FileExistsError:
        shutil.rmtree(outputdir)
        os.mkdir(outputdir)

    subgs=subgraph.graph2subgraph(repeat_graph=repeat_graph)

    for i in [6,13]:
        print('processing',i)
        subg_dir=os.path.join(outputdir,'subg%d'%i)
        subgs[i].output_all(outdir=subg_dir)
        # ans=subgs[i].get_unbranching_paths()
        # print(ans)
        sam=subgbam(subg=subgs[i],bampath=bamfile,outputdir=subg_dir,flanking=10000)
        subgs[i].updatebarcodes(samfile=sam)
        barcode_path_d,path_barcode_d=subgs[i].compute_dict()

        with open(os.path.join(subg_dir,'path_barcode'),mode='w') as f:
            for i in path_barcode_d:
                # f.write("%s\t%d\n" % (i, len(path_barcode_d[i])))
                f.write("%s\t%d\t%s\n"%(i,len(path_barcode_d[i]),path_barcode_d[i]))
        with open(os.path.join(subg_dir, 'barcode_path'), mode='w') as f:
            for i in barcode_path_d:
                f.write("%s\t%s\n"%(i,str(barcode_path_d[i])))

        # new_gv=os.path.join(subg_dir,"repeat_graph.gv")
        # addbarcode2gv(new_gv,barcodes=ans)
        # cmd = ['dot', '-Tpng', new_gv, '-o', os.path.join(subg_dir, 'repeat_graph.png')]
        # a = subprocess.Popen(cmd)
        # a.wait()

if __name__ == '__main__':
    main(sys.argv)