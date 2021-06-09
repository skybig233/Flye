#(c) 2016 by Authors
#This file is a part of ABruijn program.
#Released under the BSD license (see LICENSE file)

"""
Runs assemble binary
"""

from __future__ import absolute_import
import subprocess
import logging
import os

from flye.utils.utils import which

ASSEMBLE_BIN = "flye-modules"
logger = logging.getLogger()


class AssembleException(Exception):
    pass


def check_binaries():
    if not which(ASSEMBLE_BIN):
        raise AssembleException("Assemble binary was not found. "
                                "Did you run 'make'?")
    try:
        devnull = open(os.devnull, "w")
        subprocess.check_call([ASSEMBLE_BIN, "assemble", "-h"], stderr=devnull)
    except subprocess.CalledProcessError as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise AssembleException(str(e))
    except OSError as e:
        raise AssembleException(str(e))



def assemble(args, run_params, out_file, log_file, config_path):
    logger.info("Assembling disjointigs")
    logger.debug("-----Begin assembly log------")
    cmdline = [ASSEMBLE_BIN, "assemble", "--reads", ",".join(args.reads), "--out-asm", out_file,
               "--config", config_path, "--log", log_file, "--threads", str(args.threads)]
    # 2020/12/15 江喆圣：
    # 测试用例基础cmdline如下
    # ['flye-modules', 'assemble', '--reads',
    #  '/ldfssz1/ST_OCEAN/USER/jiangzhesheng/software/Flye/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz', '--out-asm',
    #  '/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/test_toy/00-assembly/draft_assembly.fasta', '--config',
    #  '/ldfssz1/ST_OCEAN/USER/jiangzhesheng/software/Flye/flye/config/bin_cfg/asm_corrected_reads.cfg', '--log',
    #  '/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/test_toy/flye.log', '--threads', '8']
    if args.debug:
        cmdline.append("--debug")
    if args.meta:
        cmdline.append("--meta")
    if args.genome_size:
        cmdline.extend(["--genome-size", str(args.genome_size)])
    #if args.kmer_size:
    #    cmdline.extend(["--kmer", str(args.kmer_size)])

    cmdline.extend(["--min-ovlp", str(run_params["min_overlap"])])
    if run_params["min_read_length"] > 0:
        cmdline.extend(["--min-read", str(run_params["min_read_length"])])

    if args.hifi_error:
        cmdline.extend(["--extra-params",
                        "assemble_ovlp_divergence={}".format(args.hifi_error)])

    #if args.min_kmer_count is not None:
    #    cmdline.extend(["-m", str(args.min_kmer_count)])
    #if args.max_kmer_count is not None:
    #    cmdline.extend(["-x", str(args.max_kmer_count)])

    try:
        logger.debug("Running: " + " ".join(cmdline))
        # 2020/12/15 江喆圣：
        # 传入C++的cmdline如下：
        # ['flye-modules', 'assemble',
        #  '--reads','/ldfssz1/ST_OCEAN/USER/jiangzhesheng/software/Flye/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz',
        #  '--out-asm', '/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/test_toy/00-assembly/draft_assembly.fasta',
        #  '--config','/ldfssz1/ST_OCEAN/USER/jiangzhesheng/software/Flye/flye/config/bin_cfg/asm_corrected_reads.cfg',
        #  '--log','/ldfssz1/ST_OCEAN/USER/jiangzhesheng/flye/test_toy/flye.log',
        #  '--threads', '8', '--genome-size', '500000','--min-ovlp', '1000']
        subprocess.check_call(cmdline)
    except subprocess.CalledProcessError as e:
        if e.returncode == -9:
            logger.error("Looks like the system ran out of memory")
        raise AssembleException(str(e))
    except OSError as e:
        raise AssembleException(str(e))