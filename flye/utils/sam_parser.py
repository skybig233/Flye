#(c) 2019 by Authors
#This file is a part of Flye program.
#Released under the BSD license (see LICENSE file)

"""
Provides multithreaded parser for SAM files
"""

from __future__ import absolute_import
from __future__ import division

import os
import re
import sys
from collections import namedtuple, defaultdict
import subprocess
import logging
import multiprocessing
import ctypes
import time
import gzip
import io
import random

#In Python2, everything is bytes (=str)
#In Python3, we are doing IO in bytes, but everywhere else strngs = unicode
if sys.version_info < (3, 0):
    from string import maketrans
    _STR = lambda x: x
    _BYTES = lambda x: x
else:
    maketrans = bytes.maketrans
    _STR = bytes.decode
    _BYTES = str.encode

from flye.six.moves import range
from flye.six import iteritems

import flye.utils.fasta_parser as fp

logger = logging.getLogger()

SAMTOOLS_BIN = "flye-samtools"
Alignment = namedtuple("Alignment", ["qry_id", "trg_id", "qry_start", "qry_end",
                                     "qry_sign", "qry_len", "trg_start",
                                     "trg_end", "trg_sign", "trg_len",
                                     "qry_seq", "trg_seq", "err_rate",
                                     "is_secondary"])


class AlignmentException(Exception):
    pass


class PafHit(object):
    """
    Stores paf alignment
    """
    __slots__ = ("query", "query_length", "query_start", "query_end",
                 "target", "target_length", "target_start", "target_end")
    def __init__(self, raw_hit):
        hit = raw_hit.split()

        self.query = hit[0]
        self.query_length = int(hit[1])
        self.query_start = int(hit[2])
        self.query_end = int(hit[3])

        self.target = hit[5]
        self.target_length = int(hit[6])
        self.target_start = int(hit[7])
        self.target_end = int(hit[8])

    def query_mapping_length(self):
        return self.query_end - self.query_start + 1

    def target_mapping_length(self):
        return self.target_end - self.target_start + 1

    def query_left_overhang(self):
        return self.query_start

    def query_right_overhang(self):
        return self.query_length - self.query_end + 1

    def target_left_overhang(self):
        return self.target_start

    def target_right_overhang(self):
        return self.target_length - self.target_end + 1


def read_paf(filename):
    """
    Streams out paf alignments
    """
    with open(filename, "rb") as f:
        for raw_hit in f:
            yield PafHit(_STR(raw_hit))


def read_paf_grouped(filename):
    """
    Outputs chunks of alignments for each (query, target)pair.
    Assumes that PAF alignment is already sorted by query.
    """
    prev_hit = None
    target_hits = defaultdict(list)
    for hit in read_paf(filename):
        if prev_hit is not None and hit.query != prev_hit.query:
            for trg in sorted(target_hits):
                yield target_hits[trg]
            target_hits = defaultdict(list)

        target_hits[hit.target].append(hit)
        prev_hit = hit

    if len(target_hits):
        for trg in sorted(target_hits):
            yield target_hits[trg]


class SynchronizedSamReader(object):
    """
    Parses SAM file in multiple threads.
    """
    def __init__(self, sam_alignment, reference_fasta,
                 max_coverage=None, use_secondary=False):
        #check that alignment exists
        if not os.path.exists(sam_alignment):
            raise AlignmentException("Can't open {0}".format(sam_alignment))

        #will not be changed during exceution, each process has its own copy
        self.aln_path = sam_alignment
        #self.ref_fasta = {_BYTES(h) : _BYTES(s)
        #                  for (h, s) in iteritems(reference_fasta)}
        self.fetch_list = [_BYTES(k) for k in reference_fasta.keys()]
        self.max_coverage = max_coverage
        self.use_secondary = use_secondary
        self.cigar_parser = re.compile(b"[0-9]+[MIDNSHP=X]")

        #will be shared between processes
        self.shared_manager = multiprocessing.Manager()
        self.shared_num_jobs = multiprocessing.Value(ctypes.c_int, 0)
        self.shared_lock = self.shared_manager.Lock()
        self.shared_eof = multiprocessing.Value(ctypes.c_bool, False)

        self.ref_fasta = self.shared_manager.dict()
        for (h, s) in iteritems(reference_fasta):
            self.ref_fasta[_BYTES(h)] = _BYTES(s)

        if len(self.fetch_list) == 0:
            self.shared_eof.value = True

    def close(self):
        pass

    def is_eof(self):
        return self.shared_eof.value

    def _parse_cigar(self, cigar_str, read_str, ctg_str, ctg_pos):
        #ctg_str = self.ref_fasta[ctg_name]
        trg_seq = []
        qry_seq = []
        trg_start = ctg_pos - 1
        trg_pos = ctg_pos - 1
        qry_start = 0
        qry_pos = 0

        left_hard = True
        left_soft = True
        hard_clipped_left = 0
        hard_clipped_right = 0
        soft_clipped_left = 0
        soft_clipped_right = 0
        for token in self.cigar_parser.findall(cigar_str):
            size, op = int(token[:-1]), token[-1:]
            if op == b"H":
                if left_hard:
                    qry_start += size
                    hard_clipped_left += size
                else:
                    hard_clipped_right += size
            elif op == b"S":
                qry_pos += size
                if left_soft:
                    soft_clipped_left += size
                else:
                    soft_clipped_right += size
            elif op == b"M":
                qry_seq.append(read_str[qry_pos : qry_pos + size].upper())
                trg_seq.append(ctg_str[trg_pos : trg_pos + size].upper())
                qry_pos += size
                trg_pos += size
            elif op == b"I":
                qry_seq.append(read_str[qry_pos : qry_pos + size].upper())
                trg_seq.append(b"-" * size)
                qry_pos += size
            elif op == b"D":
                qry_seq.append(b"-" * size)
                trg_seq.append(ctg_str[trg_pos : trg_pos + size].upper())
                trg_pos += size
            else:
                raise AlignmentException("Unsupported CIGAR operation: " + str(op))
            left_hard = False
            if op != b"H":
                left_soft = False

        trg_seq = b"".join(trg_seq)
        qry_seq = b"".join(qry_seq)
        matches = 0
        for i in range(len(trg_seq)):
            if trg_seq[i] == qry_seq[i]:
                matches += 1
        err_rate = 1 - matches / len(trg_seq)

        trg_end = trg_pos
        qry_end = qry_pos + hard_clipped_left
        qry_len = qry_end + hard_clipped_right
        qry_start += soft_clipped_left
        qry_end -= soft_clipped_right

        return (trg_start, trg_end, len(ctg_str), trg_seq,
                qry_start, qry_end, qry_len, qry_seq, err_rate)

    def get_chunk(self):
        """
        Gets a chunk - safe to use from multiple processes in parallel
        """
        ###
        job_id = None
        while True:
            with self.shared_lock:
                if self.shared_eof.value:
                    return None, []

                job_id = self.shared_num_jobs.value
                self.shared_num_jobs.value = self.shared_num_jobs.value + 1
                if self.shared_num_jobs.value == len(self.fetch_list):
                    self.shared_eof.value = True
                break

            time.sleep(0.01)

        parsed_contig = self.fetch_list[job_id]
        contig_str = self.ref_fasta[parsed_contig]
        chunk_buffer = []
        aln_file = subprocess.Popen(SAMTOOLS_BIN + " view '" + self.aln_path + "' '" + _STR(parsed_contig) + "'",
                                    shell=True, stdout=subprocess.PIPE).stdout
        for line in aln_file:
            chunk_buffer.append(line)
        ###

        #shuffle alignments so that they uniformly distributed. Needed for
        #max_coverage subsampling. Using the same seed for determinism
        random.Random(42).shuffle(chunk_buffer)

        sequence_length = 0
        alignments = []
        for line in chunk_buffer:
            tokens = line.strip().split()
            if len(tokens) < 11:
                #raise AlignmentException("Error reading SAM file")
                continue

            flags = int(tokens[1])
            is_unmapped = flags & 0x4
            is_secondary = flags & 0x100
            #is_supplementary = flags & 0x800    #allow supplementary
            #if is_unmapped or is_secondary: continue
            if is_unmapped: continue
            if is_secondary and not self.use_secondary: continue

            read_id = tokens[0]
            #read_contig = tokens[2]
            cigar_str = tokens[5]
            read_str = tokens[9]
            ctg_pos = int(tokens[3])
            is_reversed = flags & 0x16
            is_secondary = flags & 0x100

            if read_str == b"*":
                raise Exception("Error parsing SAM: record without read sequence")

            (trg_start, trg_end, trg_len, trg_seq,
            qry_start, qry_end, qry_len, qry_seq, err_rate) = \
                    self._parse_cigar(cigar_str, read_str, contig_str, ctg_pos)

            #OVERHANG = cfg.vals["read_aln_overhang"]
            #if (float(qry_end - qry_start) / qry_len > self.min_aln_rate or
            #        trg_start < OVERHANG or trg_len - trg_end < OVERHANG):
            aln = Alignment(_STR(read_id), _STR(parsed_contig),
                            qry_start, qry_end, "-" if is_reversed else "+", qry_len,
                            trg_start, trg_end, "+", trg_len,
                            _STR(qry_seq), _STR(trg_seq),
                            err_rate, is_secondary)
            alignments.append(aln)

            sequence_length += qry_end - qry_start
            if sequence_length // len(contig_str) > self.max_coverage:
                break

        #then, alignments by read and by score
        alignments.sort(key=lambda a: (a.qry_id, -(a.qry_end - a.qry_start)))

        #if parsed_contig is None:
        #    return None, []
        return _STR(parsed_contig), alignments


def _is_sam_header(line):
    return line[:3] in [b"@PG", b"@HD", b"@SQ", b"@RG", b"@CO"]
