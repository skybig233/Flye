# -*- coding: utf-8 -*-
# @Time : 2021/5/19 15:23
# @Author : Jiangzhesheng
# @File : solve_repeat.py
# @Software: PyCharm
import sys
import logging
from typing import List, Any, Dict, Union

from flye.repeat_graph.repeat_graph import *
from my_change import subgraph
from my_change import subgbam
import flye.trestle.graph_resolver as tres_graph

class Path:
    def __init__(self,path:List[RgEdge]=[]) -> None:
        self.edges=path

    def __add__(self, other):
        self.edges+=other.edges

    def __str__(self) -> str:
        return '->'.join([str(edge.edge_id) for edge in self.edges])

    def __len__(self)->int:
        return len(self.edges)

    def set_barcodes(self):
        ans=self.edges[0].barcodes.keys()
        for edge in self.edges:
            ans=ans&edge.barcodes.keys()

        self.barcodes=ans

    def is_end_in_repeat(self):
        return self.edges[-1].repetitive

    def have_repeat(self):
        for i in self.edges:
            if i.repetitive==True:
                return True
        return False

def get_candidate_paths(subg:RepeatGraph, start_edge:RgEdge, max_extend_length:int)->List[Path]:
    """
    给定子图，起始边，最大延伸长度，返回所有候选路径
    :param subg:子图
    :param start_edge:起始边
    :param max_extend_length:最长延申长度，不包括start_edge长度，
    即：max_extend_length<0，至少有start_edge本身
    max_extend_length=1，至少有start_edge本身和它所有的outedge
    :return:候选路径
    """
    assert start_edge in subg.edges.values(),'start_edge not in subg'
    ans: List[RgEdge]=[]
    all: List[Path]=[]
    def dfs(subg:RepeatGraph, start_edge:RgEdge, max_extend_length:int,ans:List[RgEdge]):
        """
        该dfs不能重复遍历某条edge
        深度优先遍历子图，以start_edge
        :param ans:遍历到该edge时，已有的路径，作为值传递到下一个递归中
        """
        visited=ans.copy()
        visited.append(start_edge)
        # 递归出口：1、延伸长度不够了    2、遍历到的边不在子图中
        #条件2出现的原因是抽子图时，边是直接从原图抽取
        #对于某条边a，可能在子图中已经终止，但在原图中还有后续的节点和边
        if start_edge not in subg.edges.values():
            all.append(Path(path=ans))
            return
        if max_extend_length < start_edge.length():
            all.append(Path(path=visited))
            return
        for out_edge in start_edge.node_right.out_edges:
            if out_edge not in visited:
                dfs(subg,out_edge,max_extend_length-start_edge.length(),visited)

    # def dfs_repetitive_visit(subg:RepeatGraph, start_edge:RgEdge, max_extend_length:int,ans:List[RgEdge]):
    #     """
    #     该dfs可以重复遍历某条edge，直到没有extendlength
    #     深度优先遍历子图，以start_edge
    #     :param ans:遍历到该edge时，已有的路径，作为值传递到下一个递归中
    #     """
    #     tmp=ans.copy()
    #     tmp.append(start_edge)
    #     # 递归出口：1、延伸长度不够了    2、遍历到的边不在子图中
    #     #条件2出现的原因是抽子图时，边是直接从原图抽取
    #     #对于某条边a，可能在子图中已经终止，但在原图中还有后续的节点和边
    #     if start_edge not in subg.edges.values():
    #         # all.append(ans)
    #         all.append(Path(path=ans))
    #         return
    #     if max_extend_length < start_edge.length():
    #         # all.append(tmp)
    #         all.append(Path(path=tmp))
    #         return
    #     for out_edge in start_edge.node_right.out_edges:
    #         dfs_repetitive_visit(subg,out_edge,max_extend_length-start_edge.length(),tmp)

    dfs(subg,start_edge,max_extend_length+start_edge.length(),ans)
    return all

def subg_candidate_paths(subg:RepeatGraph,max_extend_length:int)->List[Path]:
    """
    找候选路径，会过滤掉落入repeat的path
    :param subg:
    :param max_extend_length:
    :return:
    """
    ans: List[Path]=[]
    for edge in subg.edges.values():
        if edge.repetitive==False and len(subgraph.get_inedges(subg=subg,edge=edge))==0:
            ans+=get_candidate_paths(subg,edge,max_extend_length)

    # ans=[i for i in ans if i.edges[-1].repetitive==False]
    return ans

def paths_chooser(paths: List[Path],topk:int)->List[Path]:
    """
    从所有路径中选出最优路径，可以通过设置path最少需要的barcode数目限定，也可以通过topk来选
    :param paths: 
    :return: 
    """
    # best_path=[]
    # for path in paths:
    #     if len(path.barcodes)>max_barcodes:
    #         max_barcodes=len(path.barcodes)
    #         best_path=path
    # return best_path
    qualified_paths=[i for i in paths if i.have_repeat()==True]
    max_barcodes=100
    qualified_paths=[i for i in qualified_paths if len(i.barcodes)>max_barcodes]
    qualified_paths=qualified_paths[:topk]
    return qualified_paths

def analyse_path(path:Path)->List[List[Dict[str, List[int]]]]:
    """
    遍历path的所有repeat边的所有seq
    利用path.barcodes，给出所有seq的pathbarcode信息
    后续可以用于repeat挑出最好的seq
    :param path:
    :return:返回列表信息如下：
    len(list)=len(path) 列表元素个数等于path长度
    len(list[0])=len(edge.seq) 每个元素代表一条repeat edge，长度等于edge的seq个数，
    """
    #对于重复edge的每条seq都要进行打分(pathbarcode出现次数)，选取最好的seq
    edges=[]#用于记录path的每条edge情况
    for edge in path.edges:
        seqs=[]#用于记录edge的每条seq情况
        if edge.repetitive==True:
            for seq in edge.edge_sequences:
                d_barcode_locations: Dict[str, List[int]]={}#用于记录seq中的pathbarcode出现位置
                for barcode in path.barcodes:
                    if barcode in seq.barcodes:
                        d_barcode_locations[barcode]=seq.barcodes[barcode]
                seqs.append(d_barcode_locations)
        edges.append(seqs)
    return edges

def seq_splitter(path:Path,path_info)->List[str]:
    """
    根据pathbarcode质量选对应seq
    :param path:
    :param path_info:
    :return:
    """
    repeat_seq=[]
    for edge_index in range(1,len(path)-1):
        assert path.edges[edge_index].repetitive==True
        max_barcode=0#用于存储当前最大barcode数目
        choose_seq_index=0
        for seq_index in range(len(path_info[edge_index])):#遍历每个seq，如果support_barcode大于当前最大，更新当前最大，并记录seq下标
            if len(path_info[edge_index][seq_index]) > max_barcode:
                choose_seq_index=seq_index
                max_barcode=len(path_info[edge_index][seq_index])

        support_barcode_proportion=len(path_info[edge_index][choose_seq_index])/len(path.barcodes)
        if support_barcode_proportion > 0.8:
            choose_seq=path.edges[edge_index].edge_sequences[choose_seq_index].edge_seq_name
            repeat_seq.append(choose_seq)
    return repeat_seq

def seq_seq_generator(repeat_seq:List[str],subg:RepeatGraph)->str:
    d={'A':'T','T':'A','C':'G','G':'C'}
    ans=[]
    for i in repeat_seq:
        if i[0]=='+':
            ans.append(subg.edges_fasta[i[1:]])
        elif i[0]=='-':
            dna_string=subg.edges_fasta[i[1:]]
            reverse_dna_string=''.join([d[i] for i in dna_string])
            reverse_dna_string=reverse_dna_string[::-1]
            ans.append(reverse_dna_string)
        else:
            ans.append(subg.edges_fasta[i])
    return ''.join(ans)

def updategraph(path:Path, rg:RepeatGraph,logger):
    """
    采用separate_path方法更新repeat_graph
    :param path:
    :param rg:
    :return:
    """
    assert path.edges[0].repetitive==False
    assert path.edges[-1].repetitive == False


    path_info=analyse_path(path=path)

    repeat_seq=seq_splitter(path=path,path_info=path_info)
    assert repeat_seq!=[]
    path_id=[edge.edge_id for edge in path.edges]
    new_seq_id="barcode_solved_repeat_%s_%s"%(str(path.edges[0].edge_id),str(path.edges[-1].edge_id))
    new_seq_seq=seq_seq_generator(repeat_seq=repeat_seq, subg=rg)
    assert new_seq_seq!=""
    if repeat_seq==[] or new_seq_seq=="":
        logger.warning("cannot choose good seq for path %s"%str(path))
        return rg
    rg.separate_path(graph_path=path_id, new_seq_id=new_seq_id, new_seq_seq=new_seq_seq)

    return rg

def solve_mini_subg(paths,subg,origraph,logger):
    for path in paths:
        path.set_barcodes()
    paths.sort(key=lambda x:len(x.barcodes),reverse=True)
    in_egdes=[i for i in subg.edges.values() if subgraph.get_inedges(subg=subg,edge=i)==[]]
    topk=len(in_egdes)
    good_paths = paths_chooser(paths=paths, topk=topk)
    logger.info("%d good paths choosed: %s" % (len(good_paths),str([str(path) for path in good_paths])))
    good_paths=conflict_paths_resolver(good_paths,logger)
    logger.info("%d good paths with no conflict: %s" % (len(good_paths), str([str(path) for path in good_paths])))

    edges_to_remove = set()
    for path in good_paths:
        #path最终不落入repeat中才更新子图
        if path!=[] and path.edges[-1].repetitive==False:
            path_id=[i.edge_id for i in path.edges]
            edges_to_remove.update(path_id[1:-1])
            updategraph(path=path, rg=origraph,logger=logger)

    # for edge_id in edges_to_remove:
    #     try:
    #         edge = origraph.edges[edge_id]
    #         if not edge.self_complement:
    #             origraph.remove_edge(origraph.complement_edge(edge))
    #         origraph.remove_edge(edge)
    #     except KeyError:continue


def conflict_paths_resolver(paths:List[Path],logger):
    new_paths=[]
    visited_unique=set()
    for path in paths:
        if not ((path.edges[0].edge_id in visited_unique)
                and (path.edges[-1].edge_id in visited_unique)
                and (path.edges[-1].repetitive==False)):
            new_paths.append(path)

            visited_unique.add(path.edges[0].edge_id)
            visited_unique.add(path.edges[-1].edge_id)
        else:
            logger.info("conflict paths %s"%str(path))
    return new_paths

def solve_big_subg(paths,subg,origraph,logger):
    # TODO 添加大图处理
    solve_mini_subg(paths,subg,origraph,logger)

def solve_repeat_main(subgraph:RepeatGraph, bamfile:str,
                      outputdir:str, origraph:RepeatGraph,
                      logger:logging.Logger):
    """
    根据子图解原图
    :param subgraph:子图
    :param bamfile:stLFR比对到全图数据
    :param outputdir:
    :return:利用每个子图解完repeat后的原图
    """
    #根据子图选出所有可能候选路径
    paths=subg_candidate_paths(subg=subgraph, max_extend_length=40000)
    if not len(paths):
        logger.info("no candidate paths")
        return origraph
    logger.info("find %d candidate paths"%len(paths))
    for path in paths:
        logger.debug(str(path))
    #挑出子图的比对sam
    sam=subgbam.subgbam(subg=subgraph, bampath=bamfile, outputdir=outputdir, flanking=10000)
    #根据sam更新子图barcode
    logger.info("update barcode")
    subgraph.updatebarcodes(sam)

    if len(paths)<10:
        solve_mini_subg(subg=subgraph,paths=paths,origraph=origraph,logger=logger)
    else:
        solve_big_subg(subg=subgraph,paths=paths,origraph=origraph,logger=logger)
    return origraph

def main(argv):
    workdir = '/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/step_by_step'
    bamfile='/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/bwa_stLFR_flye/chr19/chr19.sort.bam'
    edge_file = os.path.join(workdir, '20-repeat-keephaplotype/repeat_graph_edges.fasta')
    graph = os.path.join(workdir, '20-repeat-keephaplotype/repeat_graph_dump_after_rr')
    repeat_graph = RepeatGraph(fp.read_sequence_dict(edge_file))
    repeat_graph.load_from_file(graph)

    outputdir="/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/bwa_stLFR_flye/tmp"

    subgs=subgraph.graph2subgraph(repeat_graph=repeat_graph)

    for i in range(len(subgs)):
        print("processing subg",i)
        subg=subgs[i]
        rg=solve_repeat_main(subgraph=subg, bamfile=bamfile, outputdir=outputdir, origraph=repeat_graph)

    rg.output_all(outputdir)

if __name__ == '__main__':
    main(sys.argv)