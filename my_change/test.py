# -*- coding: utf-8 -*-
# @Time : 2021/5/10 11:13
# @Author : Jiangzhesheng
# @File : test0.py
# @Software: PyCharm
import sys
import os

from flye.repeat_graph.repeat_graph import RepeatGraph
import flye.utils.fasta_parser as fp

from my_change.ond_algo import *

from bioclass import Sequence

def get_edgeseq_info(repeat_graph:RepeatGraph):
    unique=[]
    repeat=[]
    for edge in repeat_graph.edges.values():
        seq_num=len(edge.edge_sequences)
        try:
            cmp=edge.edge_sequences[0].edge_seq_name
            for seq in edge.edge_sequences:
                assert seq.edge_seq_name[0]==cmp[0],'+-error in %d'%edge.edge_id
            if edge.repetitive:
                if seq_num==1:
                    repeat.append(edge.edge_id)
                    raise AssertionError('repeat %d has %d seq'%(edge.edge_id,seq_num))
            else:
                if seq_num>1:
                    unique.append((edge.edge_id,seq_num))
                    raise AssertionError('unique %d has %d seq'%(edge.edge_id,seq_num))

        except AssertionError as e:
            print(e)
            pass
    return repeat,unique

def main(argv):
    workdir='/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/step_by_step'
    edge_file= os.path.join(workdir,'20-repeat-keephaplotype/repeat_graph_edges.fasta')

    before_graph=os.path.join(workdir,'20-repeat-keephaplotype/repeat_graph_dump_before_rr')
    before_repeat_graph = RepeatGraph(fp.read_sequence_dict(edge_file))
    before_repeat_graph.load_from_file(before_graph)
    repeat,unique=get_edgeseq_info(repeat_graph=before_repeat_graph)

    repeat,unique=len(repeat),len(unique)
    sum_repeat=len([i for i in before_repeat_graph.edges.values() if i.repetitive==True])
    sum_unique = len([i for i in before_repeat_graph.edges.values() if i.repetitive == False])

    print("before rr:")
    print("1 seq repeat:",repeat/sum_repeat)
    print(">2 seq unique:",unique/sum_unique)

    # seqs=repeat_graph.find_edgeid_string(id=326)
    # v=ond_algo(A=seqs[0],B=seqs[1])
    # print(get_edit_distance(v))

    # after_graph = os.path.join(workdir, '20-repeat-copy/repeat_graph_dump_after_rr')
    # after_repeat_graph = RepeatGraph(fp.read_sequence_dict(edge_file))
    # after_repeat_graph.load_from_file(after_graph)
    # repeat,unique = get_edgeseq_info(repeat_graph=after_repeat_graph)
    #
    # repeat,unique=len(repeat),len(unique)
    # sum_repeat=len([i for i in after_repeat_graph.edges.values() if i.repetitive==True])
    # sum_unique = len([i for i in after_repeat_graph.edges.values() if i.repetitive == False])
    #
    # print("after rr:")
    # print("1 seq repeat:", repeat / sum_repeat)
    # print(">2 seq unique:", unique / sum_unique)

    #
    # print(len(before_repeat_id),len(after_repeat_id))
    #
    # for i in before_repeat_id:
    #     if i not in after_repeat_id:
    #         print('error with id ',i)

if __name__ == '__main__':
    main(sys.argv)