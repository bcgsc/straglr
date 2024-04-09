import numpy as np
from collections import Counter
from .cluster import Cluster
from operator import itemgetter
from .vcf import VCF

class Variant:
    """
    0: chrom
    1: start
    2: end
    3: alleles
    4: repeat
    5: coverage
    6: genotypes
    7: genotype_summary
    8: ref_motif
    9: ref_seq
    """
    tsv_headers = ['chrom',
                   'start',
                   'end',
                   'target_repeat',
                   'locus',
                   'coverage',
                   'genotype',
                   ]

    bed_headers = ['chrom', 'start', 'end', 'repeat_unit']

    @classmethod
    def set_genotype_config(cls, method=None, min_reads=None, max_num_clusters=3, eps=None):
        cls.genotype_config = {'min_reads': 4, 'max_num_clusters': max_num_clusters}

        # minimum number of reads per cluster
        if min_reads is not None:
            cls.genotype_config['min_reads'] = min_reads

        cls.clustering = Cluster(cls.genotype_config['min_reads'], cls.genotype_config['max_num_clusters'])

    @classmethod
    def genotype(cls, variant, report_in_size=False):
        # cluster - always use sizes
        sizes = sorted([a[4] for a in variant[3] if a[-1] == 'full'])
        clusters = cls.clustering.cluster(sizes)

        # genotype labels: mean of either copy numbers(default) or size
        for cluster in clusters:
            if report_in_size:
                alleles = cluster
            else:
                alleles = [allele[3] for allele in variant[3] if allele[4] in cluster and allele[-1] == 'full']
            variant[6].append(round(np.mean(alleles), 1))

        # assign genotype to each allele
        for allele in variant[3]:
            assigned = False
            for i in range(len(clusters)):
                if allele[4] in clusters[i]:
                    allele.append(variant[6][i])
                    assigned = True
                    break

            # 'NA' assigned if read is an outlier in clustering or 'partial'
            if not assigned:
                allele.append('NA')

        partials = [r for r in variant[3] if r[-2] == 'partial']
        if partials:
            # get biggest allele size that can be clustered
            if report_in_size:
                partials_sorted = sorted(partials, key=itemgetter(4), reverse=True)
                biggest_partial_size = partials_sorted[0][4]
            else:
                partials_sorted = sorted(partials, key=itemgetter(3), reverse=True)
                biggest_partial_size = partials_sorted[0][3]

            bigger_alleles = [a for a in variant[6] if a > biggest_partial_size]
            if bigger_alleles:
                alleles_sorted = sorted(variant[6], reverse=True)
                max_sizes = []
                min_sizes = []
                i = 4 if report_in_size else 3
                for allele in alleles_sorted:
                    allele_sizes = [a[i] for a in variant[3] if a[-1] == allele]
                    max_sizes.append(max(allele_sizes))
                    min_sizes.append(min(allele_sizes))

                for p in partials_sorted:
                    allele_assigned = None
                    size = p[4] if report_in_size else p[3]
                    for i in range(len(alleles_sorted)-1):
                        if size > max_sizes[i+1]:
                            allele_assigned = alleles_sorted[i]
                            break
                    if allele_assigned is not None:
                        p[-1] = allele_assigned
                    else:
                        for i in range(len(alleles_sorted)):
                            if size >= min_sizes[i]:
                                allele_assigned = alleles_sorted[i]
                                break
                        if allele_assigned is not None:
                            p[-1] = allele_assigned

            else:
                if variant[6]:
                    max_gt_size = sorted(variant[6], reverse=True)[0]
                else:
                    max_gt_size = 0 
                # find all partials with sizes bigger than current biggest size
                biggers = [p for p in partials if p[4] > max_gt_size]

                # if there are enough partials bigger than minimum, create allele
                if len(biggers) >= cls.genotype_config['min_reads']:
                    biggest_partial_size = max([r[4] for r in biggers])
                    gt = '>{}'.format(biggest_partial_size)
                    variant[6].append(gt)
                    for p in partials:
                        p[-1] = gt
                        p[3] = round(float(p[4]) / len(variant[4]), 1) 
            
    @classmethod
    def get_genotype(cls, variant):
        allele_counts = Counter([allele[-1] for allele in variant[3] if allele[-1] != '-' and allele[-1] != 'NA'])
        # for partial
        gt = [(a, allele_counts[a]) for a in allele_counts.keys() if type(a) is str and a[0] == '>']
        for allele in sorted([a for a in allele_counts.keys() if type(a) is not str], reverse=True):
            gt.append((allele, allele_counts[allele]))
        
        return gt

    @classmethod
    def summarize_genotype(cls, variant):
        gt = cls.get_genotype(variant)
        out = []
        for allele, support in gt:
            out.append('{}({})'.format(allele, support))
        variant[7] = ';'.join(out)

    @classmethod
    def to_tsv(cls, variant):
        cols = [variant[0],
                variant[1],
                variant[2],
                variant[4],
                '{}:{}-{}'.format(variant[0], variant[1], variant[2]),
                variant[5],
                variant[7],
                ]
        return list(map(str, cols))
    
    @classmethod
    def extract_alt_motifs(cls, variant, gt):
        alts = {}
        for g in gt:
            alts[g[0]] = {}
        consensus = variant[4]
    
        for allele in variant[3]:
            if allele[-1] == 'NA' or 'failed' in allele[-2]:
                continue
            motif = allele[2]
            if motif != consensus:
                if not motif in alts[allele[-1]]:
                    alts[allele[-1]][motif] = 0
                alts[allele[-1]][motif] += 1

        am = []
        for g in gt:
            if alts[g[0]]:
                counts = []
                for motif, count in alts[g[0]].items():
                    counts.append((motif, count))
                motifs = []
                for motif, count in sorted(counts, key=itemgetter(1,0), reverse=True):
                    motifs.append('{}({})'.format(motif, count))
                am.append(','.join(motifs))
            else:
                am.append('.')
        return '/'.join(am)

    @classmethod
    def to_vcf(cls, variant, vid='.'):
        gt = cls.get_genotype(variant)
        print('vv', variant)
        cols = [variant[0],
                variant[1],
                vid,
                variant[9]]
        alt_motifs = cls.extract_alt_motifs(variant, gt)
        cols.append(VCF.create_variant_format(variant))
        cols.append(VCF.extract_variant_gt(variant, gt, alt_motifs))
        return '\t'.join(list(map(str, cols)))

    @classmethod
    def above_min_expansion(cls, variant, min_expansion, min_reads):
        ref_size = int(variant[2]) - int(variant[1]) + 1

        if variant[5]:
            n = 0
            alleles = []
            for allele in variant[6]:
                if type(allele) is str and allele[0] == '>':
                    alleles.append((allele, float(allele[1:])))
                elif type(allele) is not str:
                    alleles.append((allele, allele))
            
            for allele, size in sorted(alleles, key=itemgetter(1), reverse=True):
                reads = [a for a in variant[3] if a[-1] == allele and a[4] - ref_size >= min_expansion]
                n += len(reads)

            if n >= min_reads:
                return True
            else:
                return False
        else:
            return False

    @classmethod
    def update_coords(cls, variant):
        genome_starts = [a[5] for a in variant[3] if a[-1] == 'full']
        genome_ends = [a[6] for a in variant[3] if a[-1] == 'full']
        if genome_starts and genome_ends:
            variant[1] = int(np.median(genome_starts))
            variant[2] = int(np.median(genome_ends))

    @classmethod
    def summarize_alleles(cls, alleles):
        reads = []
        sizes = []
        cns = []
        starts = []
        for allele in alleles:
            reads.append(allele[0])
            sizes.append(str(allele[4]))
            cns.append(str(allele[3]))
            starts.append(str(allele[1]))

        return ','.join(reads), ','.join(cns), ','.join(sizes), ','.join(starts)

class Allele:
    """
    0: read
    1: rstart
    2: repeat
    3: copy_number
    4: size
    5: genome_start
    6: genome_end
    7: strand
    8: repeat_seq
    9: label
    10: genotype
    """
    tsv_headers = ['read_name',
                   'actual_repeat',
                   'copy_number',
                   'size',
                   'read_start',
                   'strand',
                   'allele',
                   'read_status',
                   ]

    summary_headers = ['reads',
                       'copy_numbers',
                       'sizes',
                       'read_starts',
                       ]

    @classmethod
    def to_tsv(cls, cols):
        # __init__ input order to output order
        cols_ordered = [cols[0],
                        cols[2],
                        cols[3],
                        cols[4],
                        cols[1],
                        cols[7],
                        cols[-1],
                        cols[-2]
                        ]
        return list(map(str, cols_ordered))
