#!/usr/bin/env python
import argparse
import pysam
from collections import defaultdict
import re

def parse_tsv(tsv, loci=None):
    support = defaultdict(dict)
    with open(tsv, 'r') as ff:
        for line in ff:
            if line[0] == '#':
                continue
            cols = line.rstrip().split('\t')
            locus = cols[4]
            status = cols[14]
            read_name = cols[7]
            size = cols[10]
            read_start = cols[11]
            strand = cols[12]

            if status != 'full' or (loci is not None and locus not in loci):
                continue
            support[locus][read_name] = int(read_start), int(size), strand

    return support

def extract_repeats(bam, support, flank_size=10):
    seqs = {}
    for locus in support:
        seqs[locus] = []
        chrom, start, end = re.split('[:-]', locus)
        for aln in bam.fetch(chrom, int(start), int(end)):
            if aln.query_name in support[locus] and not aln.query_name in seqs:
                rlen = aln.infer_read_length()
                if support[locus][aln.query_name][2] == '+':
                    start, end = support[locus][aln.query_name][0], support[locus][aln.query_name][0] + support[locus][aln.query_name][1]
                else:
                    start = rlen - (support[locus][aln.query_name][0] + support[locus][aln.query_name][1])
                    end = start + support[locus][aln.query_name][1]

                try:
                    repeat_seq = aln.query_sequence[start:end].lower()
                    left = max(0, start - flank_size), start
                    right = end, min(end + flank_size, rlen)
                    left_seq = aln.query_sequence[left[0]:left[1]].upper()
                    right_seq = aln.query_sequence[right[0]:right[1]].upper()
                    seq = left_seq + repeat_seq + right_seq
                    seqs[locus].append((aln.query_name, support[locus][aln.query_name][1], seq))
                except:
                    print('problem extracting repeat from {}'.format(aln.query_name))

    return seqs

def report(seqs, out_fa):
    with open(out_fa, 'w') as out:
        for locus in sorted(seqs.keys()):
            for read_name, repeat_size, seq in seqs[locus]:
                out.write('>{} {} {}\n{}\n'.format(read_name, locus, repeat_size, seq))

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("tsv", type=str, help="Straglr tsv")
    parser.add_argument("bam", type=str, help="bam file") 
    parser.add_argument("out", type=str, help="output fasta")
    parser.add_argument("--locus", type=str, help="UCSC-format coordinate")
    parser.add_argument("--flank", type=int, default=100, help="flank size. Default:100")
    args = parser.parse_args()
    return args
    
def main():
    args = parse_args()
    support = parse_tsv(args.tsv, args.locus)
    bam = pysam.AlignmentFile(args.bam)
    seqs = extract_repeats(bam, support, flank_size=args.flank)
    report(seqs, args.out)

main()
