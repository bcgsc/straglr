#!/usr/bin/env python
import argparse
from src.ins import INSFinder
from src.tre import TREFinder
from src.version import __version__
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", type=str, help="bam file")
    parser.add_argument("genome_fasta", type=str, help="genome_fasta")
    parser.add_argument("out", type=str, help="output file")
    parser.add_argument("--reads_fasta", type=str, help="read indexed fasta file")
    parser.add_argument("--min_ins_size", type=int, default=300, help="minimum insertion size. Default:300")
    parser.add_argument("--exclude", type=str, help="bed file to exclude regions")
    parser.add_argument("--nprocs", type=int, help="number of processes", default=1)
    parser.add_argument("--chroms", type=str, nargs="+", help="chromosomes")
    parser.add_argument("--loci", type=str, help="bed file of loci for genotyping")
    parser.add_argument("--min_support", type=int, help="minimum number of supporting reads", default=2)
    parser.add_argument("--clustering", type=str, help="clustering method. Default:gmm (alternative:dbscan)", default='gmm')
    parser.add_argument("--eps", type=int, help="epsilon parameter in dbscan clustering. Default:50", default=50)
    parser.add_argument("--genotype_in_size", action="store_true", help="report genotype in size instead of copy numbers")
    parser.add_argument("--max_str_len", type=int, help="maximum STR length. Default:50", default=50)
    parser.add_argument("--min_str_len", type=int, help="minimum STR length. Default:2", default=2)
    parser.add_argument("--max_num_clusters", type=int, help="maximum number of clusters to try. Default:3", default=3)
    parser.add_argument("--max_cov", type=int, help="maximum allowed coverage for ins inspection. Default:100", default=100)
    parser.add_argument("--debug", action='store_true', help="debug mode i.e. keep trf output")
    parser.add_argument("--version", action='version', version=__version__)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    tre_finder = TREFinder(args.bam,
                           args.genome_fasta,
                           nprocs=args.nprocs,
                           reads_fasta=args.reads_fasta,
                           max_str_len=args.max_str_len,
                           min_str_len=args.min_str_len,
                           min_support=args.min_support,
                           clustering=args.clustering,
                           genotype_in_size=args.genotype_in_size,
                           eps=args.eps,
                           max_num_clusters=args.max_num_clusters,
                           remove_tmps=not args.debug)

    variants = []
    if not args.loci:
        ins_finder = INSFinder(args.bam,
                               args.genome_fasta,
                               args.min_ins_size,
                               reads_fasta=args.reads_fasta,
                               exclude=args.exclude,
                               chroms=args.chroms,
                               nprocs=args.nprocs,
                               min_support=args.min_support,
                               max_cov=args.max_cov,
                               )
        # find ins
        ins = ins_finder.find_ins()

        # find ins that are str
        if ins:
            variants = tre_finder.examine_ins(ins, min_expansion=args.min_ins_size)

    else:
        variants = tre_finder.genotype(args.loci)

    tre_finder.output(variants, args.out)

if __name__ == '__main__':
    main()
