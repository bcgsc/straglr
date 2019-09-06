import argparse
from ins import INSFinder
from tre import TREFinder
from version import __version__
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", type=str, help="bam file")
    parser.add_argument("genome_fasta", type=str, help="genome_fasta")
    parser.add_argument("out", type=str, help="output file")
    parser.add_argument("--min_ins_size", type=int, default=300, help="minimum insertion size. Default:300")
    parser.add_argument("--exclude", type=str, help="bed file to exclude regions")
    parser.add_argument("--nprocs", type=int, help="number of processes", default=1)
    parser.add_argument("--chroms", type=str, nargs="+", help="chromosomes")
    parser.add_argument("--loci", type=str, help="bed file of loci for genotyping")
    parser.add_argument("--max_str_len", type=int, help="maximum STR length. Default:50", default=50)
    parser.add_argument("--version", action='store_true')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    
    if args.version:
        print __version__
        sys.exit()

    tre_finder = TREFinder(args.bam, 
                           args.genome_fasta,
                           nprocs=args.nprocs,
                           max_str_len=args.max_str_len)

    if not args.loci:
        ins_finder = INSFinder(args.bam,
                               args.genome_fasta,
                               args.min_ins_size,
                               exclude=args.exclude,
                               chroms=args.chroms,
                               nprocs=args.nprocs,
                               )
        # find ins
        ins = ins_finder.find_ins()

        # find ins that are str
        variants = tre_finder.examine_ins(ins)

    else:
        tre_finder.check_split_alignments = True
        variants = tre_finder.genotype(args.loci)
        
    tre_finder.output(variants, args.out)

main()