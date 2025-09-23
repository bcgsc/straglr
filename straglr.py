#!/usr/bin/env python
import argparse
from src.ins import INSFinder
from src.tre import TREFinder
from src.version import __version__
import sys
import tempfile

def parse_args():
    trf_args_meta = ('Match', 'Mismatch', 'Delta', 'PM', 'PI', 'Minscore', 'MaxPeriod')
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", type=str, help="bam file")
    parser.add_argument("genome_fasta", type=str, help="genome_fasta")
    parser.add_argument("out_prefix", type=str, help="output prefix")

    scan = parser.add_argument_group('genome scan')
    scan.add_argument("--min_ins_size", type=int, default=100, help="minimum insertion size. Default:100")
    scan.add_argument("--min_str_len", type=int, help="minimum STR length. Default:2", default=2)
    scan.add_argument("--max_str_len", type=int, help="maximum STR length. Default:50", default=50)
    scan.add_argument("--min_support", type=int, help="minimum number of supporting reads for detecting expansion. Default:2", default=2)
    scan.add_argument("--trf_args", type=int, nargs=7, help="tandem repeat finder arguments. Default:2 5 5 80 10 10 500", metavar=trf_args_meta, default=[2,5,5,80,10,10,500])
    scan.add_argument("--chroms", type=str, nargs="+", help="chromosomes")
    scan.add_argument("--regions", type=str, help="bed file for scanning only specific regions")
    scan.add_argument("--exclude", type=str, help="bed file to exclude regions")
    scan.add_argument("--include_alt_chroms", action='store_true', help="include alternate chromosomes. By default, only chroms 1-22,X,Y are considered in genome scan")
    scan.add_argument("--max_cov", type=int, help="maximum allowed coverage for ins inspection. Default:100", default=100)
    scan.add_argument("--use_unpaired_clips", action='store_true', help="also examine unpaired clipped alignments in genome scan")

    genotyping = parser.add_argument_group('genotyping')
    genotyping.add_argument("--loci", type=str, help="bed file of loci for genotyping")
    genotyping.add_argument("--genotype_in_size", action="store_true", help="report genotype in size instead of copy numbers")
    genotyping.add_argument("--include_partials", action='store_true', help="detect and report reads only capturing partial repeats when genotyping")
    genotyping.add_argument("--use_mean", action="store_true", help="use mean instead of median cluster value for defining allele")
    genotyping.add_argument("--flank_size", type=int, default=80, help="flanking sequence length to add to loci when genotyping. Default:80")

    io = parser.add_argument_group('input/output')
    io.add_argument("--sample", type=str, help="sample name for VCF output", default='.')
    io.add_argument("--sex", type=str, help="sex, m(ale) or f(emale)", default='f')
    io.add_argument("--reads_fasta", type=str, nargs='+', help="read indexed fasta file(s)")
    io.add_argument("--nprocs", type=int, help="number of processes", default=1)
    io.add_argument("--tmpdir", type=str, help="directory to use for generating tmp files instead of system TEMP")
    io.add_argument("--debug", action='store_true', help="debug mode i.e. keep trf output")
    io.add_argument("--version", action='version', version=__version__)

    clustering = parser.add_argument_group('clustering (for defining alleles)')
    clustering.add_argument("--min_cluster_size", type=float, help="minimum number or fraction of supporting reads for allele clustering. Default:0.1", default=0.1)
    clustering.add_argument("--max_num_clusters", type=int, help="maximum number of clusters to try. Default:2", default=2)
    clustering.add_argument("--max_check_size", type=int, help="threshold of cluster (allele) size (median/mean) to skip filtering by --min_cluster_size. Default:5000", default=5000)
    clustering.add_argument("--min_cluster_d", type=int, help="minimum separation between clusters. Default:10", default=10)
    clustering.add_argument("--max_bad_cluster_size", type=int, help="maximum cluster size for checking quality. Default:5", default=5)

    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    if args.tmpdir:
        tempfile.tempdir = args.tmpdir

    #min_cluster_size = args.min_cluster_size if args.min_cluster_size < args.min_support else args.min_support

    sex = args.sex[0].lower()
    if sex != 'm':
        sex = 'f'
    tre_finder = TREFinder(args.bam,
                           args.genome_fasta,
                           args.flank_size,
                           nprocs=args.nprocs,
                           reads_fasta=args.reads_fasta,
                           max_str_len=args.max_str_len,
                           min_str_len=args.min_str_len,
                           min_support=args.min_support,
                           min_cluster_size=args.min_cluster_size,
                           max_check_size=args.max_check_size,
                           min_cluster_d=args.min_cluster_d,
                           genotype_in_size=args.genotype_in_size,
                           max_num_clusters=args.max_num_clusters,
                           max_bad_cluster_size=args.max_bad_cluster_size,
                           use_mean=args.use_mean,
                           trf_args=' '.join(map(str, args.trf_args + ['-d', '-h'])),
                           include_partials=args.include_partials,
                           sample=args.sample,
                           sex=sex,
                           debug=args.debug)

    variants = []
    if not args.loci:
        ins_finder = INSFinder(args.bam,
                               args.genome_fasta,
                               args.min_ins_size,
                               reads_fasta=args.reads_fasta,
                               exclude=args.exclude,
                               chroms=args.chroms,
                               include_alt_chroms=args.include_alt_chroms,
                               nprocs=args.nprocs,
                               min_support=args.min_support,
                               max_cov=args.max_cov,
                               use_unpaired_clips=args.use_unpaired_clips,
                               debug=args.debug,
                               )
        # find ins
        ins = ins_finder.find_ins(regions_bed_file=args.regions)

        # find ins that are str
        if ins:
            variants = tre_finder.examine_ins(ins, min_expansion=args.min_ins_size)

    else:
        tre_finder.min_cluster_size = args.min_cluster_size
        variants = tre_finder.genotype(args.loci)

    # output both bed and tsv
    tre_finder.output_bed(variants, '{}.bed'.format(args.out_prefix))
    tre_finder.output_tsv(variants, '{}.tsv'.format(args.out_prefix), cmd=' '.join(sys.argv))
    tre_finder.output_vcf(variants, '{}.vcf'.format(args.out_prefix))

if __name__ == '__main__':
    main()
