import argparse
from pybedtools import BedTool
from collections import defaultdict
from scipy import stats
import numpy as np
from operator import itemgetter
import sys

def create_bed(alleles):
    bed_str = ''
    loci = defaultdict(dict)
    max_num_alleles = 0
    for locus in alleles:
        if len(alleles[locus]) > max_num_alleles:
            max_num_alleles = len(alleles[locus])
        for allele in alleles[locus]:
            loci[locus][float(allele)] = ','.join(alleles[locus][allele])

    ncols = 5 + 2 * max_num_alleles
    for locus in loci:
        cols = list(locus)
        for allele in sorted(loci[locus]):
            cols.extend([str(allele), loci[locus][allele]])
        if len(cols) < ncols:
            cols.extend(['-'] * (ncols - len(cols)))
        bed_str += '{}\n'.format('\t'.join(cols))

    return BedTool(bed_str, from_string=True).sort()

def parse_straglr_tsv(tsv, use_size=True, skip_chroms=None):
    alleles = {}
    with open(tsv, 'r') as ff:
        for line in ff:
            if line[0] == '#':
                continue

            cols = line.rstrip().split('\t')
            if skip_chroms is not None and cols[0] in skip_chroms:
                continue
            locus = cols[:4]
            ref_size = int(cols[2]) - int(cols[1]) + 1
            if not use_size:
                ref_size = '{:.1f}'.format(ref_size / len(cols[3]))
            locus.append(str(ref_size))
            locus = tuple(locus)
            if not locus in alleles:
                alleles[locus] = defaultdict(list)
            allele = cols[-1]
            if allele == '-':
                continue
            if not use_size:
                allele = '{:.1f}'.format(float(allele) / len(cols[3]))
            size = cols[-3]
            copy_num = cols[-4]
            if use_size:
                alleles[locus][allele].append(size)
            else:
                alleles[locus][allele].append(copy_num)

    return create_bed(alleles)

def vs_each_parent(proband_bed, parent_bed, pval_cutoff, min_expansion=100, min_support=0):
    expanded_loci = {}
    new_loci = []
    common_loci = []
    for cols in proband_bed.intersect(parent_bed, loj=True, wao=True, f=0.9, F=0.9).saveas('abc.bed'):
        if cols[-2] == '.':
            new_loci.append(cols)
        else:
            common_loci.append(cols)

    for cols in new_loci:
        #locus = tuple(cols[:5])
        locus = (cols[0], int(cols[1]), int(cols[2]), cols[3], float(cols[4]))
        ref_size = locus[-1]
        proband_alleles = []
        supports = {}

        for i in range(5, 8, 2):
            if cols[i] == '-' or float(cols[i]) - ref_size < min_expansion:
                continue
            proband_calls = cols[i+1].split(',')
            if len(proband_calls) < min_support:
                continue
            proband_alleles.append(float(cols[i]))
            supports[float(cols[i])] = len(proband_calls)
        expanded_loci[locus] = proband_alleles, ['-'], '-', supports

    for cols in common_loci:
        proband_cols = cols[:9]
        locus = (proband_cols[0], int(proband_cols[1]), int(proband_cols[2]), proband_cols[3], float(proband_cols[4]))
        ref_size = locus[-1]

        expanded_alleles = []
        supports = {}
        parent_alleles = []
        pvals = []
        # loop thru proband allele
        for i in range(6, 9, 2):
            if proband_cols[i] == '-':
                continue
            proband_calls = list(map(float, proband_cols[i].split(',')))
            if len(proband_calls) < min_support:
                continue

            proband_allele = float(proband_cols[i-1])
            if float(proband_cols[i-1]) - ref_size < min_expansion:
                continue

            # loop thru parents
            not_expanded = False
            parent_cols = cols[9:9+10]
            parent_pvals = []

            # loop thru parent alleles
            for j in range(-3, 0, 2):
                if parent_cols[j] == '-':
                    continue
               
                parent_allele = float(parent_cols[j-1])
                if not parent_allele in parent_alleles:
                    parent_alleles.append(parent_allele)
                parent_calls = list(map(float, parent_cols[j].split(',')))

                result = stats.ttest_ind(np.array(parent_calls), np.array(proband_calls), alternative='less') 

                if result.pvalue > pval_cutoff or proband_allele - parent_allele + 1 < min_expansion:
                    not_expanded = True
                else:
                    parent_pvals.append(np.format_float_scientific(result.pvalue, precision=2))

            if not not_expanded:
                #print('expanded', locus, proband_allele, proband_calls, parent_cols, parent_pvals, type(parent_pvals[0]), proband_allele, parent_alleles, proband_allele-max(parent_alleles))
                expanded_alleles.append(proband_allele)
                supports[proband_allele] = len(proband_calls)
                pvals.append(','.join(parent_pvals))
            #else:
                #print('not_expanded', locus, proband_allele, proband_calls, parent_cols)
            
        if expanded_alleles:
            expanded_loci[locus] = expanded_alleles, [','.join(list(map(str, parent_alleles)))], ';'.join(pvals), supports

    return expanded_loci

def vs_all_parents(vs_parents):
    if len(vs_parents) == 1:
        return vs_parents[0]

    expanded_loci = {}
    for locus in set(vs_parents[0].keys()) & set(vs_parents[1].keys()):
        p1, p2 = vs_parents[0][locus], vs_parents[1][locus]
        expanded_alleles = list(set(p1[0]) & set(p2[0]))
        if not expanded_alleles:
            continue

        supports = []
        for allele in expanded_alleles:
            if allele in p1[3]:
                supports.append(p1[3][allele])
            elif allele in p2[3]:
                supports.append(p2[3][allele])

        parent_alleles = []
        pvals = []
        for p in (p1,p2):
            if p[1]:
                parent_alleles.append(p[1][0])
                pvals.append(p[2])
            else:
                parent_alleles.append('-')
                pvals.append('-')
        expanded_loci[locus] = list(expanded_alleles), parent_alleles, pvals, ','.join(map(str, list(supports)))

    return expanded_loci

def output(expanded_loci, out_file):
    header = ('chrom', 'start', 'end', 'repeat', 'ref_size', 'proband_allele', 'proband_allele_support', 'parent_alleles', 'pvals')
    with open(out_file, 'w') as out:
        out.write('#{}\n'.format('\t'.join(header)))
        for locus in sorted(expanded_loci.keys(), key=itemgetter(0,1,2)):
            expanded_alleles = ','.join(list(map(str, expanded_loci[locus][0])))
            parent_alleles = ';'.join(expanded_loci[locus][1])
            if parent_alleles == '-;-':
                parent_alleles = '-'
            pvals = ';'.join(expanded_loci[locus][2])
            if pvals == '-;-':
                pvals = '-'
            cols = list(map(str, locus)) + [expanded_alleles, expanded_loci[locus][3], parent_alleles, pvals]
            out.write('{}\n'.format('\t'.join(cols)))

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("proband", type=str, help="proband straglr results")
    parser.add_argument("parents", type=str, nargs='+', help="parent straglr results (maximum=2)")
    parser.add_argument("output", type=str, help="output")
    parser.add_argument("--use_size", action='store_true', help="use size")
    parser.add_argument("--min_expansion", type=int, default=0, help="minimum expansion")
    parser.add_argument("--min_support", type=int, default=0, help="minimum support")
    parser.add_argument("--skip_chroms", type=str, nargs='+', help="skip chromosomes")
    parser.add_argument("--pval_cutoff", type=float, default=0.001,  help="p-value cutoff for testing T-test hypothesis")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()

    if len(args.parents) > 2:
        print('more than 2 parent results given, only first 2 will be used')
        parents = args.parents[:2]
    else:
        parents = args.parents
    proband_bed = parse_straglr_tsv(args.proband, use_size=args.use_size, skip_chroms=args.skip_chroms)

    vs_parents = []
    for parent in parents:
        parent_bed = parse_straglr_tsv(parent, use_size=args.use_size)
        vs_parents.append(vs_each_parent(proband_bed, parent_bed, args.pval_cutoff, min_support=args.min_support))

    expanded_loci = vs_all_parents(vs_parents)
    output(expanded_loci, args.output)

main()

