import numpy as np
from sklearn.cluster import KMeans, AgglomerativeClustering
from collections import defaultdict, Counter
from scipy import stats
from .cluster import Cluster

class Variant:
    """
    0: chrom
    1: start
    2: end
    3: alleles
    4: repeat
    5: genotypes
    6: genotype_summary
    """
    tsv_headers = ['chrom',
                   'start',
                   'end',
                   'repeat_unit',
                   'genotype',
                   ]

    @classmethod
    def set_genotype_config(cls, method=None, min_reads=None, max_num_clusters=3, eps=None):
        genotype_config = {'method': 'gmm',
                           'min_reads': 4,
                           'max_num_clusters': max_num_clusters,
                           'eps': 50}

        # default = 'gmm', alternative: 'dbscan'
        if method == 'dbscan':
            genotype_config['method'] = method

        # minimum number of reads per cluster
        if min_reads is not None:
            genotype_config['min_reads'] = min_reads

        # default 50 is used if not specified
        if eps is not None:
            genotype_config['eps'] = eps

        cls.clustering = Cluster(genotype_config['method'],
                                 genotype_config['min_reads'],
                                 genotype_config['max_num_clusters'],
                                 genotype_config['eps'])

    @classmethod
    def genotype(cls, variant, report_in_size=False):
        # cluster - always use sizes
        sizes = sorted([a[4] for a in variant[3]])
        clusters = cls.clustering.cluster(sizes)

        # genotype labels: mean of either copy numbers(default) or size
        for cluster in clusters:
            if report_in_size:
                alleles = cluster
            else:
                alleles = [allele[3] for allele in variant[3] if allele[4] in cluster]
            variant[5].append(round(np.mean(alleles), 1))

        # assign genotype to each allele
        for allele in variant[3]:
            assigned = False
            for i in range(len(clusters)):
                if allele[4] in clusters[i]:
                    allele.append(variant[5][i])
                    assigned = True
                    break

            # '-' assigned if read is an outlier in clustering
            if not assigned:
                allele.append('-')

    @classmethod 
    def summarize_genotype(cls, variant):
        allele_counts = Counter([allele[-1] for allele in variant[3]])
        out = []
        for allele in sorted([a for a in allele_counts.keys() if type(a) is not str], reverse=True) +\
                             [a for a in allele_counts.keys() if type(a) is str]:
            if allele == '-' and len(allele_counts.keys()) > 1:
                continue
            out.append('{}({})'.format(allele, allele_counts[allele]))
        variant[6] = ';'.join(out)

    @classmethod
    def to_tsv(cls, variant):
        sorted_genotypes = sorted(variant[5], reverse=True)
        cols = [variant[0],
                variant[1],
                variant[2],
                variant[4],
                variant[6],
                ]
        return list(map(str, cols))

    @classmethod
    def above_min_expansion(cls, variant, min_expansion, min_reads):
        ref_size = int(variant[2]) - int(variant[1]) + 1

        if variant[5]:
            max_allele_size = sorted(variant[5])[-1]

            alleles = [a for a in variant[3] if a[7] == max_allele_size and a[4] - ref_size >= min_expansion]
            #sizes = [a[4] for a in variant[3] if a[7] == max_allele_size]
            #print('cc', variant, min_reads, len(variant[3]), len(alleles), len(alleles)>min_reads, np.median(sizes), np.median(sizes) - ref_size >min_expansion)
            return len(alleles) > min_reads
            #return max_allele_size - ref_size >= min_expansion
        else:
            return False

    @classmethod
    def update_coords(cls, variant):
        genome_starts = [a[5] for a in variant[3]]
        genome_ends = [a[6] for a in variant[3]]
        if genome_starts and genome_ends:
            variant[1] = int(np.median(genome_starts))
            variant[2] = int(np.median(genome_ends))

class Allele:
    """
    0: read
    1: rstart
    2: repeat
    3: copy_number
    4: size
    5: genome_start
    6: genome_end
    7: genotype
    """
    tsv_headers = ['read',
                   'copy_number',
                   'size',
                   'read_start',
                   'allele',
                   ]

    @classmethod
    def to_tsv(cls, cols):
        # __init__ input order to output order
        cols_ordered = [cols[0],
                        cols[3],
                        cols[4],
                        cols[1],
                        cols[7],
                        ]
        return list(map(str, cols_ordered))
