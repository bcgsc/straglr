import numpy as np
from sklearn.cluster import KMeans, AgglomerativeClustering
from collections import defaultdict
from scipy import stats

class Variant:
    tsv_headers = ['chrom',
                   'start',
                   'end',
                   'repeat_unit',
                   'support',
                   'gt1',
                   'gt1_support',
                   'gt2',
                   'gt2_support',
                   ]

    def __init__(self, chrom, start, end, repeat_unit):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.repeat_unit = repeat_unit
        self.alleles = []
        self.genotypes = []

    def update_coords(self):
        if self.alleles:
            # use median instead of min() and max() because of outliers
            self.start = int(np.median([a.chrom_start for a in self.alleles if a.chrom_start is not None]))
            self.end = int(np.median([a.chrom_end for a in self.alleles if a.chrom_end is not None]))

    def genotype(self):
        def cluster_into_2(cn):
            X = np.array(cn).reshape(-1,1)

            # k-means
            #km = KMeans(n_clusters=2)
            #km = km.fit(X)
            #genotypes = [c[0] for c in km.cluster_centers_]
            #clusters = defaultdict(list)
            #for i in range(len(cn)):
                #clusters[genotypes[km.labels_[i]]].append(cn[i])

            # hclust
            hclust = AgglomerativeClustering(n_clusters=2, affinity='euclidean', linkage='ward')
            clusters = defaultdict(list)
            labels = hclust.fit_predict(X)
            for i in range(len(cn)):
                clusters[labels[i]].append(cn[i])

            return clusters

        def verify_2_clusters(c1, c2, threshold=0.05):
            stat, pval = stats.ttest_ind(np.array(c1), np.array(c2), equal_var=False)

            same = True
            if pval < threshold:
                same = False

            return stat, pval, same

        # extract all copy numbers
        if len(self.alleles) > 2:
            cns = [a.copy_number for a in self.alleles]
            clusters = cluster_into_2(cns)
            means = sorted(clusters.keys())

            same = True
            if len(means) >= 2:
                if len(clusters[means[0]]) > 1 and len(clusters[means[1]]) > 1:
                    stat, val, same = verify_2_clusters(clusters[means[0]], clusters[means[1]])

            if same:
                self.genotypes = [round(np.median(cns), 1)]
                for allele in self.alleles:
                    allele.genotype = round(self.genotypes[0], 1)
            else:
                self.genotypes = [round(np.median(clusters[means[0]]), 1), round(np.median(clusters[means[1]]), 1)]
                cn_to_genotype_index = {}
                for i in range(len(means)):
                    for cn in clusters[means[i]]:
                        cn_to_genotype_index[cn] = i
                for allele in self.alleles:
                    allele.genotype = self.genotypes[cn_to_genotype_index[allele.copy_number]]

    def as_tsv(self):
        sorted_genotypes = sorted(self.genotypes, reverse=True)
        cols = [self.chrom,
                self.start,
                self.end,
                self.repeat_unit,
                len(self.alleles),
                sorted_genotypes[0],
                len([a for a in self.alleles if a.genotype == sorted_genotypes[0]]),
                sorted_genotypes[1] if len(sorted_genotypes) == 2 else '-',
                len([a for a in self.alleles if a.genotype == sorted_genotypes[1]]) if len(sorted_genotypes) == 2 else '-',
                ]

        return map(str, cols)

class Allele:
    tsv_headers = ['read',
                   'copy_number',
                   'size',
                   #'repeat_unit',
                   'read_start',
                   'genotype',
                   'chrom_start',
                   'chrom_end']

    def __init__(self, read, read_pos, repeat_unit, copy_number, seq, chrom_start, chrom_end):
        self.read = read
        self.read_pos = read_pos
        self.repeat_unit = repeat_unit
        self.copy_number = copy_number
        self.seq = seq
        self.chrom_start = chrom_start
        self.chrom_end = chrom_end
        self.genotype = None

    def as_tsv(self):
        cols = [self.read,
                self.copy_number,
                len(self.seq),
                #self.repeat_unit,
                self.read_pos,
                self.genotype,
                self.chrom_start,
                self.chrom_end]

        return map(str, cols)