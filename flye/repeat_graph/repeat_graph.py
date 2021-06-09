#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
This module provides repeat graph parsing/serializing functions,
(as output by repeat graph construction module)
as well as some basic operations
"""

from __future__ import division
from typing import Dict, List, Set, Any
import os
import shutil
import subprocess
from collections import defaultdict

import flye.utils.fasta_parser as fp

from bioclass.sam import GeneralSamFile



class RgNode(object):
    __slots__ = ("in_edges", "out_edges")

    def __init__(self):
        self.in_edges = []
        self.out_edges = []

    def is_bifurcation(self):
        return len(self.in_edges) != 1 or len(self.out_edges) != 1

class EdgeSequence(object):
    __slots__ = ("edge_seq_name", "edge_seq_len", "orig_seq_id", "orig_seq_len",
                 "orig_seq_start", "orig_seq_end","barcodes")

    def __init__(self, edge_seq_name=None, edge_seq_len=0, orig_seq_id="*",
                 orig_seq_len=0, orig_seq_start=0, orig_seq_end=0):
        self.edge_seq_name = edge_seq_name
        self.edge_seq_len = edge_seq_len

        self.orig_seq_id = orig_seq_id
        self.orig_seq_len = orig_seq_len
        self.orig_seq_start = orig_seq_start
        self.orig_seq_end = orig_seq_end

        barcodes: Dict[str, List[str]]

        self.barcodes={}# key=barcode value=[positions]

    def addbarcode(self,barcode,position)->None:
        """
        对seq添加barcode
        如果barcode已经出现，则记录不同的出现位置
        :param barcode:
        :param position:
        :return:
        """
        self.barcodes.setdefault(barcode,[]).append(position)

class RgEdge(object):
    edge_sequences: List[EdgeSequence]
    __slots__ = ("node_left", "node_right", "edge_id", "repetitive",
                 "self_complement", "resolved", "mean_coverage",
                 "alt_group", "edge_sequences")

    def __init__(self, node_left:RgNode=None, node_right:RgNode=None,edge_id=None):
        self.node_left = node_left
        self.node_right = node_right
        self.edge_id = edge_id

        self.repetitive = False
        self.self_complement = False
        self.resolved = False
        self.mean_coverage = 0
        self.edge_sequences = []
        self.alt_group = -1

    @property
    def barcodes(self)->Dict[str,int]:
        """
        对于每个Edge而言，就是把每个seq的barocdes合并，相同barcodes数目相加
        :return:
        """
        d: Dict[str, List[List[int]]]={}
        for seq in self.edge_sequences:
            for k,v in seq.barcodes.items():
                d.setdefault(k,[]).append(v)
        return d

    def length(self):
        if not self.edge_sequences:
            return 0

        return sum([s.edge_seq_len
                    for s in self.edge_sequences]) // len(self.edge_sequences)

    def __repr__(self):
        return "(id={0}, len={1}, cov={2} rep={3})" \
                .format(self.edge_id, self.length(),
                        self.mean_coverage, self.repetitive)

    def adjacentEdges(self):
        """
        :return:返回一条边所有的邻近边的集合，不包括自身
        """
        edges=set()
        for e in self.node_left.in_edges:
            edges.add(e)
        for e in self.node_left.out_edges:
            edges.add(e)
        for e in self.node_right.in_edges:
            edges.add(e)
        for e in self.node_right.out_edges:
            edges.add(e)
        edges.remove(self)
        return edges

class RepeatGraph(object):
    __slots__ = ("nodes", "edges", "edges_fasta")

    def __init__(self, edges_fasta:Dict[str,str],
                 edges:Dict[int,RgEdge]={},
                 nodes:List[RgNode]=[]):
        self.nodes = nodes
        self.edges = edges #key = edge id
        self.edges_fasta = edges_fasta #这是一个字典，key是edge_disjointigid，value是fastastr

    def add_node(self)->RgNode:
        self.nodes.append(RgNode())
        return self.nodes[-1]

    def add_edge(self, edge):
        self.edges[edge.edge_id] = edge
        edge.node_left.out_edges.append(edge)
        edge.node_right.in_edges.append(edge)

    def remove_edge(self, edge):
        _remove_from_list(edge.node_left.out_edges, edge)
        _remove_from_list(edge.node_right.in_edges, edge)
        del self.edges[edge.edge_id]

    def complement_edge(self, edge:RgEdge)->RgEdge:
        if edge.self_complement:
            return edge
        return self.edges[-edge.edge_id]

    def complement_node(self,node):
        if node.out_edges!=[]:
            return self.complement_edge(node.out_edges[0]).node_right
        elif node.in_edges!=[]:
            return self.complement_edge(node.in_edges[0]).node_left
        return None

    def get_unbranching_paths(self):
        unbranching_paths = []
        visited_edges = set()

        for edge in self.edges.values():
            if edge in visited_edges:
                continue

            traversed = [edge]
            if not edge.self_complement:
                cur_node = edge.node_left
                while (not cur_node.is_bifurcation() and
                       len(cur_node.in_edges) > 0 and
                       cur_node.in_edges[0] not in visited_edges and
                       not cur_node.in_edges[0].self_complement):
                    traversed.append(cur_node.in_edges[0])
                    visited_edges.add(cur_node.in_edges[0])
                    cur_node = cur_node.in_edges[0].node_left

                traversed = traversed[::-1]
                cur_node = edge.node_right

                while (not cur_node.is_bifurcation() and
                       len(cur_node.out_edges) > 0 and
                       cur_node.out_edges[0] not in visited_edges and
                       not cur_node.out_edges[0].self_complement):
                    traversed.append(cur_node.out_edges[0])
                    visited_edges.add(cur_node.out_edges[0])
                    cur_node = cur_node.out_edges[0].node_right

            unbranching_paths.append(traversed)

        return unbranching_paths

    def load_from_file(self, filename):
        id_to_node = {}
        cur_edge = None
        with open(filename, "r") as f:
            for line in f:
                tokens = line.strip().split()
                if tokens[0] == "Edge":
                    (edge_id, left_node, right_node, repetitive,
                     self_complement, resolved, mean_coverage, alt_group) = tokens[1:]
                    if left_node not in id_to_node:
                        id_to_node[left_node] = self.add_node()
                    if right_node not in id_to_node:
                        id_to_node[right_node] = self.add_node()

                    cur_edge = RgEdge(id_to_node[left_node],
                                      id_to_node[right_node],
                                      _to_signed_id(int(edge_id)))
                    cur_edge.repetitive = bool(int(repetitive))
                    cur_edge.self_complement = bool(int(self_complement))
                    cur_edge.resolved = bool(int(resolved))
                    cur_edge.mean_coverage = int(mean_coverage)
                    cur_edge.alt_group = int(alt_group)
                    self.add_edge(cur_edge)

                elif tokens[0] == "Sequence":
                    (edge_seq_name, edge_seq_len, orig_seq_id,
                    orig_seq_len, orig_seq_start, orig_seq_end) = tokens[1:]
                    edge_seq = EdgeSequence(edge_seq_name, int(edge_seq_len),
                                            orig_seq_id, orig_seq_len,
                                            orig_seq_start, orig_seq_end)
                    cur_edge.edge_sequences.append(edge_seq)

                else:
                    raise Exception("Error parsing " + filename)

    def dump_to_file(self, filename):
        next_node_id = 0
        node_ids = {}
        for node in self.nodes:
            node_ids[node] = next_node_id
            next_node_id += 1

        with open(filename, "w") as f:
            for edge in self.edges.values():
                f.write("Edge\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n"
                    .format(_to_unsigned_id(edge.edge_id), node_ids[edge.node_left],
                            node_ids[edge.node_right], int(edge.repetitive),
                            int(edge.self_complement), int(edge.resolved),
                            int(edge.mean_coverage), int(edge.alt_group)))

                for seq in edge.edge_sequences:
                    f.write("\tSequence\t{0} {1} {2} {3} {4} {5}\n"
                        .format(seq.edge_seq_name, seq.edge_seq_len,
                                seq.orig_seq_id, seq.orig_seq_len,
                                seq.orig_seq_start, seq.orig_seq_end))

    def output_dot(self, filename):
        next_node_id = 0
        node_ids = {}
        for node in self.nodes:
            node_ids[node] = next_node_id
            next_node_id += 1

        with open(filename, "w") as f:
            f.write("digraph {\nnodesep = 0.5;\n"
                    "node [shape = circle, label = \"\", height = 0.3];\n")
            for edge in self.edges.values():
                f.write("{0} -> {1} [label = \"id {2}\\l{3} {4}x\", color = \"{5}\"]\n"
                    .format(node_ids[edge.node_left], node_ids[edge.node_right],
                            edge.edge_id, edge.length(),edge.mean_coverage,"red" if edge.repetitive else "black"))
            f.write("}")

    def output_svg(self,prefixname:str):
        dot_file=prefixname+'.gv'
        svg_file=prefixname+'.svg'
        self.output_dot(filename=dot_file)
        cmd=['dot','-Tsvg',dot_file,'-o',svg_file]
        a=subprocess.Popen(cmd)
        a.wait()

    def output_edgefasta(self,filename):
        fp.write_fasta_dict(filename=filename,fasta_dict=self.edges_fasta)

    def output_all(self,outdir):

        try:
            os.mkdir(outdir)
        except FileExistsError:
            pass

        self.dump_to_file(filename=os.path.join(outdir,'repeat_graph_dump'))
        self.output_edgefasta(filename=os.path.join(outdir,'repeat_graph_edges.fasta'))
        self.output_svg(prefixname=os.path.join(outdir,'repeat_graph'))

    def compute_dict(self)->Dict[str, Set[int]]:
        """
        给定一个repeatgraph，计算所有path包含的barcode个数
        :return:dict    key=repeatgraph中某个path  value=该path包含的barcode数目
        """
        barcode_path_d={}# key=barcode  value=该barcode所支持的路径
        for edge in self.edges.values():
            for barcode in edge.barcodes:
                if barcode in barcode_path_d:
                    barcode_path_d[barcode].update({edge.edge_id:edge.barcodes[barcode]})
                else:
                    barcode_path_d[barcode]={edge.edge_id:edge.barcodes[barcode]}
        barcode_path_d={k:v for k,v in barcode_path_d.items() if len(v)>1}

        path_barcode_d={}
        for i in barcode_path_d:
            key=str(list(barcode_path_d[i].keys()))
            path_barcode_d.setdefault(key,list()).append(i)
        # max_val=max([len(i) for i in ans.values()])
        # ans=dict(filter(lambda x:len(x[1])>=max_val,ans.items()))
        return barcode_path_d,path_barcode_d

    def separate_path(self, graph_path:List[int], new_seq_id:str, new_seq_seq:str):
        """
        Separates the path (and its complement) on the graph.
        First and last edges in the path are disconnected from the graph,
        and then connected by a new edge. For example,
        a path (A -> X -> Y -> ... -> Z -> B) will be transformed
        into a new unbranching path A -> N -> B,
        where N represents a new edge with the given sequence.
        The intermediate path edges remain in the graph (their mean coverage 
        is modified accordingly) and they acquire the attribute 'resolved'.
        Resolved edges could later be cleaned up by using XXX function.
        :param graph_path:路径包含的边id列表
        """

        def separate_one(edges_path, new_edge_id, new_edge_seq):
            left_node = self.add_node()
            _remove_from_list(edges_path[0].node_right.in_edges, edges_path[0])
            edges_path[0].node_right = left_node
            left_node.in_edges.append(edges_path[0])

            path_coverage = (edges_path[0].mean_coverage +
                             edges_path[-1].mean_coverage) // 2
            for mid_edge in edges_path[1:-1]:
                mid_edge.resolved = True
                mid_edge.mean_coverage -= path_coverage

            right_node = left_node
            if len(edges_path) > 2:
                right_node = self.add_node()
                new_edge = RgEdge(left_node, right_node, new_edge_id)
                self.add_edge(new_edge)
                new_edge.mean_coverage = path_coverage
                new_edge.edge_sequences.append(new_edge_seq)

            _remove_from_list(edges_path[-1].node_left.out_edges, edges_path[-1])
            edges_path[-1].node_left = right_node
            right_node.out_edges.append(edges_path[-1])

        if len(graph_path) < 2:
            raise Exception("Path is too short")

        fwd_edges = []
        rev_edges = []
        for e in graph_path:
            if e not in self.edges:
                raise Exception("Nonexistent edge")
            fwd_edges.append(self.edges[e])
            rev_edges.append(self.complement_edge(self.edges[e]))
        rev_edges = rev_edges[::-1]

        new_edge_seq = EdgeSequence("+" + new_seq_id, len(new_seq_seq))
        compl_edge_seq = EdgeSequence("-" + new_seq_id, len(new_seq_seq))
        self.edges_fasta[new_seq_id] = new_seq_seq

        new_edge_id = max(self.edges.keys()) + 1
        separate_one(fwd_edges, new_edge_id, new_edge_seq)
        separate_one(rev_edges, -new_edge_id, compl_edge_seq)

    def find_edgeid_string(self,id:int)->List[str]:
        ans=[]
        for edge in self.edges.values():
            if edge.edge_id==id:
                for seq in edge.edge_sequences:
                    ans.append(self.edges_fasta[seq.edge_seq_name[1:]])
        return ans

    def updatebarcodes(self,samfile:str):
        for sam in GeneralSamFile(path=samfile):
            barcode = sam.query_name.split('#')[1]
            edge_id=int(sam.ref_name.split('_')[1])
            edge_seq_id=int(sam.ref_name.split('_')[2])
            if edge_id in self.edges:
                self.edges[edge_id].edge_sequences[edge_seq_id].\
                    addbarcode(barcode=barcode,position=sam.position)
            if -edge_id in self.edges:
                self.edges[-edge_id].edge_sequences[edge_seq_id]. \
                    addbarcode(barcode=barcode, position=sam.position)

def _remove_from_list(lst, elem):
    lst = [x for x in lst if x != elem]


def _to_signed_id(unsigned_id):
    return -(unsigned_id + 1) // 2 if unsigned_id % 2 else unsigned_id // 2 + 1


def _to_unsigned_id(signed_id):
    unsigned_id = abs(signed_id) * 2 - 2
    return unsigned_id + int(signed_id < 0)