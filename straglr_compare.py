import argparse
from pybedtools import BedTool
from collections import defaultdict
from scipy import stats
import numpy as np
from operator import itemgetter
import sys
import os

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

def vs_each_control(test_bed, control_bed, pval_cutoff, min_expansion=100, min_support=0):
    expanded_loci = {}
    new_loci = []
    common_loci = []
    for cols in test_bed.intersect(control_bed, loj=True, wao=True, f=0.9, F=0.9):
        if cols[-2] == '.':
            new_loci.append(cols)
        else:
            common_loci.append(cols)

    for cols in new_loci:
        #locus = tuple(cols[:5])
        locus = (cols[0], int(cols[1]), int(cols[2]), cols[3], float(cols[4]))
        ref_size = locus[-1]
        test_alleles = []
        supports = {}

        for i in range(5, 8, 2):
            if cols[i] == '-' or float(cols[i]) - ref_size < min_expansion:
                continue
            test_calls = cols[i+1].split(',')
            if len(test_calls) < min_support:
                continue
            test_alleles.append(float(cols[i]))
            supports[float(cols[i])] = len(test_calls)
        expanded_loci[locus] = test_alleles, ['-'], '-', supports

    for cols in common_loci:
        test_cols = cols[:9]
        locus = (test_cols[0], int(test_cols[1]), int(test_cols[2]), test_cols[3], float(test_cols[4]))
        ref_size = locus[-1]

        expanded_alleles = []
        supports = {}
        control_alleles = []
        pvals = []
        # loop thru test allele
        for i in range(6, 9, 2):
            if test_cols[i] == '-':
                continue
            test_calls = list(map(float, test_cols[i].split(',')))
            if len(test_calls) < min_support:
                continue

            test_allele = float(test_cols[i-1])
            if float(test_cols[i-1]) - ref_size < min_expansion:
                continue

            # loop thru controls
            not_expanded = False
            control_cols = cols[9:9+10]
            control_pvals = []

            # loop thru control alleles
            for j in range(-3, 0, 2):
                if control_cols[j] == '-':
                    continue
               
                control_allele = float(control_cols[j-1])
                if not control_allele in control_alleles:
                    control_alleles.append(control_allele)
                control_calls = list(map(float, control_cols[j].split(',')))

                result = stats.ttest_ind(np.array(control_calls), np.array(test_calls), alternative='less') 

                if result.pvalue > pval_cutoff or test_allele - control_allele + 1 < min_expansion:
                    not_expanded = True
                else:
                    control_pvals.append(np.format_float_scientific(result.pvalue, precision=2))

            if not not_expanded:
                #print('expanded', locus, test_allele, test_calls, control_cols, control_pvals, type(control_pvals[0]), test_allele, control_alleles, test_allele-max(control_alleles))
                expanded_alleles.append(test_allele)
                supports[test_allele] = len(test_calls)
                pvals.append(','.join(control_pvals))
            #else:
                #print('not_expanded', locus, test_allele, test_calls, control_cols)
            
        if expanded_alleles:
            expanded_loci[locus] = expanded_alleles, [','.join(list(map(str, control_alleles)))], ';'.join(pvals), supports

    return expanded_loci

def vs_all_controls(vs_controls):
    if len(vs_controls) == 1:
        return vs_controls[0]

    expanded_loci = {}
    for locus in set.intersection(*map(set, vs_controls)):
        expanded_alleles = set.intersection(*map(set, [vc[locus][0] for vc in vs_controls]))
        if not expanded_alleles:
            continue

        supports = []
        for allele in expanded_alleles:
            for vc in vs_controls:
                if allele in vc[locus][3]:
                    supports.append(vc[locus][3][allele])
                    break

        control_alleles = []
        pvals = []
        for vc in vs_controls:
            if vc[locus][1]:
                control_alleles.append(vc[locus][1][0])
                pvals.append(vc[locus][2])
            else:
                control_alleles.append('-')
                pvals.append('-')
        #print('bb', locus, supports, control_alleles, pvals)
        expanded_loci[locus] = list(expanded_alleles), control_alleles, pvals, ','.join(map(str, list(supports)))

    return expanded_loci

def output(expanded_loci, out_file):
    header = ('chrom', 'start', 'end', 'repeat', 'ref_size', 'test_allele', 'test_allele_support', 'control_alleles', 'pvals')
    with open(out_file, 'w') as out:
        out.write('#{}\n'.format('\t'.join(header)))
        for locus in sorted(expanded_loci.keys(), key=itemgetter(0,1,2)):
            expanded_alleles = ','.join(list(map(str, expanded_loci[locus][0])))
            control_alleles = ';'.join(expanded_loci[locus][1])
            if control_alleles == '-;-':
                control_alleles = '-'
            pvals = ';'.join(expanded_loci[locus][2])
            if pvals == '-;-':
                pvals = '-'
            cols = list(map(str, locus)) + [expanded_alleles, expanded_loci[locus][3], control_alleles, pvals]
            out.write('{}\n'.format('\t'.join(cols)))

def parse_list(list_file):
    return None

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("test", type=str, help="Straglr results of test_sample")
    parser.add_argument("controls", type=str, nargs='+', help="controls")
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

    control_results = []
    control_bams = []
    print(args.controls, type(args.controls))
    if type(args.controls) == str:
        parse_list(args.controls)
    elif type(args.controls) == list:
        print(os.path.splitext(args.controls[0])[1])
        if os.path.splitext(args.controls[0])[1] == '.tsv':
            control_results = args.controls
        elif os.path.splitext(args.controls[0])[1] == '.bam':
            control_bams = args.controls

    if not control_results and not control_bams:
        sys.exit()

    expanded_loci = {}
    test_bed = parse_straglr_tsv(args.test, use_size=args.use_size, skip_chroms=args.skip_chroms)
    if control_results:
        vs_controls = []
        for control_result in control_results:
            print('cc', control_result)
            control_bed = parse_straglr_tsv(control_result, use_size=args.use_size)
            vs_controls.append(vs_each_control(test_bed, control_bed, args.pval_cutoff, min_support=args.min_support))
    
        if vs_controls:
            expanded_loci = vs_all_controls(vs_controls)

    if expanded_loci:
        output(expanded_loci, args.output)
    ''' 
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
    '''
main()
