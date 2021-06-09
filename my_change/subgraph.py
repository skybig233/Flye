# -*- coding: utf-8 -*-
# @Time : 2021/3/22 10:58
# @Author : Jiangzhesheng
# @File : subgraph.py
# @Software: PyCharm
from flye.repeat_graph.repeat_graph import *
import flye.utils.fasta_parser as fp
import os
import matplotlib.pyplot as plt
from typing import List,Set
import math

#只用作类型提示，相等于typedef
RepeatCluster=Set[RgEdge]

def add_edges_fasta(edges:RepeatCluster,all_edges_fasta:Dict[str,str]):
    """
    从所有的fasta序列中找到给定edges集合中包含的fasta序列
    :param edges:
    :param all_edges_fasta:
    :return:
    """
    all_keys=set()
    for edge in edges:
        # edge_keys=[i.edge_seq_name[1:] for i in edge.edge_sequences]
        edge_keys=set()

        for seq in edge.edge_sequences:
            key=seq.edge_seq_name[1:]
            edge_keys.add(key)

        all_keys.update(edge_keys)

    values=[all_edges_fasta[i] for i in all_keys]

    return edges,{k:v for k,v in zip(all_keys,values)}

def get_repeat_clusters(repeat_graph:RepeatGraph)->List[RepeatCluster]:
    """
    获取repeat_clusters
    repeat_clusters是由多个repeat_cluster构成
    每个repeat_cluster是由一个或多个repeat边构成
    如果A、B两条邻近边都被标记为repeat，那么A、B两条会出现在同一个repeat_cluster中
    从repeat_graph中获取repeat_clusters，和C++文件里的outputdot算法思路一样
    :return: List[set[RgEdge]]，一个set即为一个cluster，一个cluster中含有不重复的RgEdge
    """
    # TODO:能否以边为单位遍历
    #  这里是以节点为单位进行遍历，判断节点的出度入度边是否是repetitive
    #  但也可以直接遍历所有边，如果边是repetitive，拉入所有adjcent边，直到stack没有重复边
    dfsstack = []
    visited = []
    repeat_clusters = []
    for node in repeat_graph.nodes:
        if node in visited: continue
        dfsstack.insert(0, node)
        repeat_cluster=set()
        while dfsstack != []:
            curNode = dfsstack.pop()
            if curNode in visited: continue

            visited.append(curNode)

            adjEdge: RgEdge
            for adjEdge in curNode.out_edges:
                if adjEdge.repetitive:
                    repeat_cluster.add(adjEdge)
                    if adjEdge.node_right not in visited:
                        dfsstack.append(adjEdge.node_right)

            for adjEdge in curNode.in_edges:
                if adjEdge.repetitive:
                    repeat_cluster.add(adjEdge)
                    if adjEdge.node_left not in visited:
                        dfsstack.append(adjEdge.node_left)

        if len(repeat_cluster) != 0:
            repeat_clusters.append(repeat_cluster)


    return repeat_clusters

def add_unique_edge(repeat_cluster:RepeatCluster)->RepeatCluster:
    """
    对repeat_cluster中的边进行一次扩充，为了更清晰的看到repeat结构
    :param repeat_cluster:不包含unique边的集合
    :return:包含unique边的集合
    """
    ans=repeat_cluster.copy()
    for edge in repeat_cluster:
        if edge.repetitive:
            for adjedge in edge.adjacentEdges():
                ans.add(adjedge)
    return ans

def repeat_cluster2repeat_graph(repeat_cluster:RepeatCluster,edges_fasta:Dict[str,str])->RepeatGraph:
    """
    单独的repeat_cluster即可构成一个repeatgraph，是大的graph下的一个子图
    :param repeat_cluster: 由RgEdge构成的set
    :return: 由repeat_cluster构成的子图
    """
    # 定义一个空的repeatgraph
    tmp=RepeatGraph(edges_fasta=edges_fasta)
    node_set=[]
    edge_set={}
    # 对repeat_cluster中的边进行遍历
    for edge in repeat_cluster:
        node_set.append(edge.node_right)
        node_set.append(edge.node_left)
        edge_set[edge.edge_id]=edge
    node_set=list(set(node_set))
    tmp.edges=edge_set
    tmp.nodes=node_set
    return tmp

def graph2subgraph(repeat_graph:RepeatGraph)->List[RepeatGraph]:
    """
    输入一个repeat_graph，根据repeat情况，将其拆成许多子图
    :param repeat_graph:
    :return:各个subgraph
    """
    ans=[]
    # 首先获取repeat_graph中所有repeat边，相连接的repeat放在一起
    repeat_clusters=get_repeat_clusters(repeat_graph)

    # 将repeat边的集合进行一个unique_edge的扩充
    repeat_clusters=list(map(add_unique_edge,repeat_clusters))

    # 对扩充后的边集找fasta序列
    for i in range(len(repeat_clusters)):
        tmp=add_edges_fasta(edges=repeat_clusters[i],
                            all_edges_fasta=repeat_graph.edges_fasta)
        repeat_clusters[i]=tmp


    for repeat_cluster,edges_fasta in repeat_clusters:
        subgraph=repeat_cluster2repeat_graph(repeat_cluster=repeat_cluster,edges_fasta=edges_fasta)
        ans.append(subgraph)
    return ans

def draw_SD(repeat_graph:RepeatGraph, outputdir:str):
    """
    给定一个repeat_graph，输出SD信息和length信息
    :param repeat_graph:
    :param outputdir:
    :return:
    """
    repeat_clusters = get_repeat_clusters(repeat_graph)
    sd_data=[len(repeat_cluster) for repeat_cluster in repeat_clusters]
    edges=[i for i in repeat_graph.edges.values()]
    repeat_edge=filter(lambda x:x.repetitive==True,edges)
    length_data=[i.length() for i in repeat_edge]
    # 写频次文件
    sd_file=os.path.join(outputdir,'sd_info.txt')
    length_file = os.path.join(outputdir, 'length_info.txt')

    with open(sd_file,mode='w') as sdf,open(length_file,mode='w') as lf:
        d={}
        sdf.write(str(sd_data))
        lf.write(str(length_data))
    # 绘制频次图
    plt.hist(x=sd_data,bins=max(sd_data))
    plt.xlabel(xlabel='SD_complexity')
    plt.ylabel(ylabel='Frequency')
    plt.title(label='SD_INFO')
    plt.savefig(os.path.join(outputdir,'sd_info.png'))
    plt.show()

    log_10_length_data=[math.log10(i) for i in length_data]
    # bins=int(max(log_10_length_data))+1
    plt.hist(x=log_10_length_data)
    plt.xlabel(xlabel='log_length')
    plt.ylabel(ylabel='Frequency')
    plt.title(label='length_INFO')
    plt.savefig(os.path.join(outputdir, 'length_info.png'))
    plt.show()

def get_inedges(subg:RepeatGraph,edge:RgEdge)->List[RgEdge]:
    """
    获取子图的起始边
    :param subg:
    :param edge:
    :return:
    """
    assert edge in subg.edges.values()
    ans=[]
    for in_edge in edge.node_left.in_edges:
        if in_edge in subg.edges.values():
            ans.append(in_edge)
    return ans

def main():
    workdir='/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/step_by_step'
    edge_file= os.path.join(workdir,'20-repeat-copy/repeat_graph_edges.fasta')
    graph=os.path.join(workdir,'20-repeat-copy/repeat_graph_dump_after_rr')
    repeat_graph = RepeatGraph(fp.read_sequence_dict(edge_file))
    repeat_graph.load_from_file(graph)

    outputdir=os.path.join(workdir,'22-subgraph-after-rr')

    subgs=graph2subgraph(repeat_graph=repeat_graph)

    #以下是数据测试代码
    # r_string=[repeat_graph.edges_fasta[i.edge_seq_name[1:]] for i in r]
    # in_string=[repeat_graph.edges_fasta[i.edge_seq_name[1:]] for i in in_edge]
    # out_string=[repeat_graph.edges_fasta[i.edge_seq_name[1:]] for i in out_edge]
    #
    # # new_contig=in_string[0]+r_string[0]+out_string[0]+r_string[1]+in_string[1]
    # # a=in_string[0]+r_string[0]+out_string[0]+r_string[1]+r_string[2]+in_string[1]
    # a=in_string[0]
    # b=out_string[0]
    # c=in_string[1]
    # d=out_string[1]
    # r0=r_string[0]
    # r1=r_string[1]
    # r2=r_string[2]
    # print('a=%d\nb=%d\nc=%d\nd=%d\nr0=%d\tr1=%d\tr2=%d'
    #       %(len(a),len(b),len(c),len(d),len(r0),len(r1),len(r2)))
    # # test=b+r0+a+r1+c
    # test1=a+r0+b[0:10000]
    # test2=b[20000:30000]+r2+b[50000:60000]
    # contig=[test1,test2]
    #
    # outputdir=os.path.join(workdir,'10-test-consensus')
    # with open(os.path.join(outputdir,'10-test.fasta'),mode='w') as f:
    #     cnt=0
    #     for i in contig:
    #         cnt+=1
    #         f.write('>disjointig_%d\n'%cnt)
    #         f.write(i+'\n')

if __name__ == '__main__':
    main()