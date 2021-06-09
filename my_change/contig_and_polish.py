# -*- coding: utf-8 -*-
# @Time : 2021/6/7 10:15
# @Author : Jiangzhesheng
# @File : contig_and_polish.py
# @Software: PyCharm
import sys
import subprocess
from flye.main import JobContigger,JobPolishing,JobFinalize,JobConfigure
import os
import argparse
import flye.config.py_cfg as cfg
FLYE="/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/debug/bin/"

def contigger_and_polish(args):
    # new_edges="/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/solve_only/repeat_graph_edges.fasta"
    # new_graph="/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/solve_only/repeat_graph_dump"
    # new_out="/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/solve_only/30-contigger"
    # new_log="/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/solve_only/flye.log"
    # cmd=[FLYE+"flye-modules","contigger",
    #      "--graph-edges",new_edges,
    #      "--reads","/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/data/human/chr19.lr.fa",
    #      "--out-dir",new_out,
    #      "--config","/ldfssz1/ST_OCEAN/USER/jiangzhesheng/software/anaconda3/lib/python3.7/site-packages/flye/config/bin_cfg/asm_raw_reads.cfg",
    #      "--repeat-graph",new_graph,
    #      "--graph-aln","/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/real_data_run/human/pb_ccs_chr19/20-repeat/read_alignment_dump",
    #      "--log",new_log,
    #      "--threads","8","--min-ovlp","5000"]
    # pp=subprocess.Popen(cmd)
    # pp.wait()

    work_dir=args.out_dir
    log_file=args.log_file

    jobconfig=JobConfigure(args, work_dir)
    jobconfig.run()

    #Contigger
    repeat_graph_edges="/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/solve_only/repeat_graph_edges.fasta"
    repeat_graph="/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/solve_only/repeat_graph_dump"
    reads_alignment="/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/real_data_run/human/pb_ccs_chr19/20-repeat/read_alignment_dump"
    jobcontig=JobContigger(args, work_dir, log_file, repeat_graph_edges,
                             repeat_graph, reads_alignment)
    raw_contigs = jobcontig.out_files["contigs"]
    scaffold_links = jobcontig.out_files["scaffold_links"]
    graph_file = jobcontig.out_files["assembly_graph"]
    gfa_file = jobcontig.out_files["gfa_graph"]
    final_graph_edges = jobcontig.out_files["edges_sequences"]
    repeat_stats = jobcontig.out_files["stats"]

    jobcontig.run()

    #Polishing
    contigs_file = raw_contigs
    polished_stats = None
    polished_gfa = gfa_file
    if args.num_iters > 0:
        jobpolish=JobPolishing(args, work_dir, log_file, raw_contigs,
                                 final_graph_edges, gfa_file)
        contigs_file = jobpolish.out_files["contigs"]
        polished_stats = jobpolish.out_files["stats"]
        polished_gfa = jobpolish.out_files["polished_gfa"]

    jobpolish.run()

    #Report results
    jobfinal=JobFinalize(args, work_dir, log_file, contigs_file,
                            graph_file, repeat_stats, polished_stats,
                            polished_gfa, scaffold_links)
    jobfinal.run()

def main():
    def check_int_range(value, min_val, max_val, require_odd=False):
        ival = int(value)
        if ival < min_val or ival > max_val:
            raise argparse.ArgumentTypeError("value should be in the "
                            "range [{0}, {1}]".format(min_val, max_val))
        if require_odd and ival % 2 == 0:
            raise argparse.ArgumentTypeError("should be an odd number")
        return ival

    parser = argparse.ArgumentParser \
        (description="Assembly of long reads with repeat graphs",
         formatter_class=argparse.RawDescriptionHelpFormatter)

    read_group = parser.add_mutually_exclusive_group(required=True)
    read_group.add_argument("--pacbio-raw", dest="pacbio_raw",
                        default=None, metavar="path", nargs="+",
                        help="PacBio raw reads")
    read_group.add_argument("--pacbio-corr", dest="pacbio_corrected",
                        default=None, metavar="path", nargs="+",
                        help="PacBio corrected reads")
    read_group.add_argument("--pacbio-hifi", dest="pacbio_hifi",
                        default=None, metavar="path", nargs="+",
                        help="PacBio HiFi reads")
    read_group.add_argument("--nano-raw", dest="nano_raw", nargs="+",
                        default=None, metavar="path",
                        help="ONT raw reads")
    read_group.add_argument("--nano-corr", dest="nano_corrected", nargs="+",
                        default=None, metavar="path",
                        help="ONT corrected reads")
    read_group.add_argument("--subassemblies", dest="subassemblies", nargs="+",
                        default=None, metavar="path",
                        help="high-quality contigs input")
    parser.add_argument("-g", "--genome-size", dest="genome_size",
                        metavar="size", required=False, default=None,
                        help="estimated genome size (for example, 5m or 2.6g)")
    parser.add_argument("-o", "--out-dir", dest="out_dir",
                        default=None, required=True,
                        metavar="path", help="Output directory")
    parser.add_argument("-t", "--threads", dest="threads",
                        type=lambda v: check_int_range(v, 1, 128),
                        default=1, metavar="int", help="number of parallel threads [1]")
    parser.add_argument("-i", "--iterations", dest="num_iters",
                        type=lambda v: check_int_range(v, 0, 10),
                        default=1, help="number of polishing iterations [1]",
                        metavar="int")
    parser.add_argument("-m", "--min-overlap", dest="min_overlap", metavar="int",
                        type=lambda v: check_int_range(v, 1000, 10000),
                        default=None, help="minimum overlap between reads [auto]")
    parser.add_argument("--asm-coverage", dest="asm_coverage", metavar="int",
                        default=None, help="reduced coverage for initial "
                        "disjointig assembly [not set]", type=int)
    parser.add_argument("--hifi-error", dest="hifi_error", metavar="float",
                        default=None, help="expected HiFi reads error rate (e.g. 0.01 or 0.001)"
                        " [0.01]", type=float)
    parser.add_argument("--plasmids", action="store_true",
                        dest="plasmids", default=False,
                        help="rescue short unassembled plasmids")
    parser.add_argument("--meta", action="store_true",
                        dest="meta", default=False,
                        help="metagenome / uneven coverage mode")
    parser.add_argument("--keep-haplotypes", action="store_true",
                        dest="keep_haplotypes", default=False,
                        help="do not collapse alternative haplotypes")
    parser.add_argument("--trestle", action="store_true",
                        dest="trestle", default=False,
                        help="enable Trestle [disabled]")
    parser.add_argument("--polish-target", dest="polish_target",
                        metavar="path", required=False,
                        help="run polisher on the target sequence")
    parser.add_argument("--resume", action="store_true",
                        dest="resume", default=False,
                        help="resume from the last completed stage")
    parser.add_argument("--resume-from", dest="resume_from", metavar="stage_name",
                        default=None, help="resume from a custom stage")
    parser.add_argument("--stop-after", dest="stop_after", metavar="stage_name",
                        default=None, help="stop after the specified stage completed")
    #parser.add_argument("--kmer-size", dest="kmer_size",
    #                    type=lambda v: check_int_range(v, 11, 31, require_odd=True),
    #                    default=None, help="kmer size (default: auto)")
    parser.add_argument("--debug", action="store_true",
                        dest="debug", default=False,
                        help="enable debug output")

    args = parser.parse_args()

    if args.asm_coverage and (args.genome_size is None):
        parser.error("--asm-coverage option requires genome size estimate (--genome-size)")

    if args.asm_coverage and args.meta:
        parser.error("--asm-coverage is incompatible with --meta")

    if args.hifi_error and not args.pacbio_hifi:
        parser.error("--hifi-error can only be used with --pacbio-hifi")

    #if not args.genome_size and not args.polish_target:
    #    parser.error("Genome size argument (-g/--genome-size) "
    #                 "is required for assembly")

    if args.pacbio_raw:
        args.reads = args.pacbio_raw
        args.platform = "pacbio"
        args.read_type = "raw"
    if args.pacbio_corrected:
        args.reads = args.pacbio_corrected
        args.platform = "pacbio"
        args.read_type = "corrected"
    if args.pacbio_hifi:
        args.reads = args.pacbio_hifi
        args.platform = "pacbio"
        args.read_type = "hifi"
    if args.nano_raw:
        args.reads = args.nano_raw
        args.platform = "nano"
        args.read_type = "raw"
    if args.nano_corrected:
        args.reads = args.nano_corrected
        args.platform = "nano"
        args.read_type = "corrected"
    if args.subassemblies:
        args.reads = args.subassemblies
        args.platform = "pacbio"    #arbitrary
        args.read_type = "subasm"

    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)
    args.out_dir = os.path.abspath(args.out_dir)

    args.log_file = os.path.join(args.out_dir, "flye.log")

    args.asm_config = os.path.join(cfg.vals["pkg_root"],
                                   cfg.vals["bin_cfg"][args.read_type])

    contigger_and_polish(args)

    return 0

if __name__ == '__main__':
    bin_absolute="/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/debug/bin"
    os.environ["PATH"] = bin_absolute + os.pathsep + os.environ["PATH"]
    main()