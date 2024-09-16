import numpy as np
from collections import Counter, defaultdict
from .cluster import Cluster
from operator import itemgetter
from .vcf import VCF
from .utils import is_same_repeat

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
    10: ref_cn
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
    def genotype(cls, variant, sex=None, report_in_size=False):
        # cluster - always use sizes
        sizes = sorted([a[4] for a in variant[3] if a[-1] == 'full'])
        max_num_clusters = 1 if sex is not None and sex.lower() == 'm' and variant[0] in ('chrX', 'X') else None
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
        ''' for vcf '''
        alts = {}
        for g in gt:
            alts[g] = {}
        ref = variant[8]

        motif_alias = {}
        for allele in variant[3]:
            if allele[-1] == 'NA' or 'failed' in allele[-2]:
                continue
            motif = allele[2]
            if motif != ref:
                if not alts[allele[-1]]:
                    alts[allele[-1]][motif] = 1
                else:
                    same = False
                    for mm in alts[allele[-1]]:
                        if is_same_repeat((mm, motif)):
                            same = True
                            alts[allele[-1]][mm] += 1
                            if mm != motif:
                                motif_alias[motif] = mm
                            break
                    # new, not alias of previous
                    if not same:
                        alts[allele[-1]][motif] = 1

        am = []
        for g in gt:
            if alts[g]:
                counts = []
                for motif, count in alts[g].items():
                    counts.append((motif, count))
                motifs = []
                for motif, count in sorted(counts, key=itemgetter(1,0), reverse=True):
                    motifs.append('{}({})'.format(motif, count))
                am.append(';'.join(motifs))
            else:
                am.append('.')
        return '/'.join(am), motif_alias

    @classmethod
    def get_alts(cls, variant, genotype_in_size, motif_alias, symbolic=False):
        ''' for vcf '''
        choices = defaultdict(list)
        for a in variant[3]:
            if a[9] is None or a[-2] != 'full':
                continue
            if a[9] > 0:
                if genotype_in_size:
                    diff = abs(a[-1] - a[4])
                else:
                    diff = abs(a[-1] - a[3])
                motif = a[2]
                if a[2] in motif_alias:
                    motif = motif_alias[a[2]]
                choices[a[9]].append((a[0], a[8], diff, motif))
        
        alts = []
        for gt in sorted(choices.keys()):
            # pick most common motif
            motifs = [a[3] for a in choices[gt]]
            motif = Counter(motifs).most_common(1)[0][0]

            choices_sorted = sorted([a for a in choices[gt] if a[3] == motif], key=itemgetter(2,0))
            rep = choices_sorted[0]
            if not symbolic:
                # sequence
                alt = rep[1]
            else:
                # copy number (symbolic)
                alt = '<TR{}>'.format(int(len(rep[1]) / len(rep[-1])))
            alts.append(alt)
        if alts:
            return ','.join(alts)
        else:
            return '.'

    @classmethod
    def convert_gt(cls, gt, variant, genotype_in_size):
        ''' for vcf '''
        col = 3 if genotype_in_size else 4
        gt2 = []
        for g in gt:
            alleles = [a[col] for a in variant[3] if a[-1] == g and a[-1] != '-' and a[-1] != 'NA' and a[-2] == 'full']
            gt2.append(round(np.mean(alleles), 1))
        return gt2
    
    @classmethod
    def get_allele_ranges(cls, gt, variant):
        ''' for vcf '''
        size_ranges = []
        cn_ranges = []
        for g in gt:
            alleles = [a for a in variant[3] if a[-1] == g and a[-1] != '-' and a[-1] != 'NA' and a[-2] == 'full']
            sizes = [a[4] for a in alleles]
            cns = [a[3] for a in alleles]
            size_ranges.append('{}-{}'.format(min(sizes), max(sizes)))
            cn_ranges.append('{}-{}'.format(min(cns), max(cns)))
        return size_ranges, cn_ranges

    @classmethod
    def find_fails(cls, variants):
        failed = [v for v in variants if not v[6]]
        fails = {}
        for variant in failed:
            failed_reasons = [a[10].split('(')[1].rstrip(')') for a in variant[3] if a[10] and 'failed' in a[10]]
            locus = tuple(map(str, variant[:3]))
            # no failed reason, but no genotype, most probable reason: No cluster was formed
            failed_reason = 'clustering_failed'
            if failed_reasons:
                failed_reason = Counter(failed_reasons).most_common(1)[0][0]
            fails[locus] = failed_reason

        return fails

    @classmethod
    def extract_gts(cls, variant, genotype_in_size):
        gt = cls.get_genotype(variant)

    @classmethod
    def to_vcf(cls, variant, genotype_in_size, vid='.', sex=None, fail=None, locus_id=None, symbolic=False):
        gt = cls.get_genotype(variant)
        
        gt_larger = [g for g in gt if type(g[0]) is str]
        gt_normal = [g for g in gt if type(g[0]) is not str]

        gt_sorted = sorted(gt_normal, key=itemgetter(0))
        gt1 = [g[0] for g in gt_sorted]
        if gt_larger:
            gt1 = [g[0] for g in gt_larger] + gt1
        gt2 = cls.convert_gt(gt1, variant, genotype_in_size)
        if gt_larger:
            for i in range(len(gt_larger)):
                gt2[i] = '>' + str(gt2[i])
        
        supports = [g[1] for g in gt_sorted]

        size_ranges, cn_ranges = cls.get_allele_ranges(gt1, variant)

        (gt_size, gt_cn) = (gt1, gt2) if genotype_in_size else (gt2, gt1)

        # found alt motifs first to find alias
        alt_motifs, motif_alias = cls.extract_alt_motifs(variant, gt1)
        has_alt_motifs = True if alt_motifs not in ('.', './.') else False
        if has_alt_motifs:
            symbolic = False
        alts = cls.get_alts(variant, genotype_in_size, motif_alias, symbolic=symbolic)
        qual = '.'
        filter = 'PASS'
        if fail is not None:
            filter = VCF.extract_filter(fail)
            filter = filter[0]

        cols = [variant[0],
                variant[1],
                vid,
                variant[9],
                alts,
                qual,
                filter]
        cols.append(VCF.extract_variant_info(variant, locus_id))
        cols.extend(VCF.extract_variant_gt(variant, gt_size, gt_cn, supports, size_ranges, cn_ranges, alt_motifs, sex=sex))
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
    9: gt
    10: label
    11: genotype
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
