import numpy as np
from numpy.random import rand
from sklearn.mixture import GaussianMixture
from scipy import stats
from operator import itemgetter
from collections import defaultdict
import argparse

class Cluster:
    def __init__(self, min_pts, max_num_clusters=2, max_check_size=5000):
        self.min_pts = min_pts
        self.max_num_clusters = max_num_clusters
        self.max_check_size = max_check_size

    def cluster(self, data):
        if len(data) == 1:
            return [[data[0]]]
        clusters = self.gmm(data, self.max_num_clusters + 4)

        self.merge_clusters(clusters, d=10, max_n=self.max_num_clusters)

        min_pts = len(data) * self.min_pts if self.min_pts < 1 and self.min_pts > 0 else self.min_pts 
        clusters = [c for c in clusters if not (len(c) < min_pts and np.median(c) < self.max_check_size)]
        return clusters

    def gmm(self, x, max_num_clusters):
        X = np.array([[i] for i in x])
        m = max_num_clusters if max_num_clusters < len(x) else len(x)
        N = np.arange(1, m + 1)
        models = [None for i in range(len(N))]
        AICs = []
        BICs = []
        if len(x) == 1:
            return []
        if len(set(x)) == 1:
            return [x]
        for i in range(len(N)):
            models[i] = GaussianMixture(N[i], init_params="k-means++", random_state=42).fit(X)
        AIC = [m.aic(X) for m in models]
        BIC = [m.bic(X) for m in models]
        M_best = models[np.argmin(AIC)]
        labels = M_best.predict(X)

        nclusters = M_best.n_components
        clusters_dict = defaultdict(list)
        for label, i in zip(labels, X):
            clusters_dict[label].append(i[0])
        clusters = [sorted(clusters_dict[i]) for i in clusters_dict.keys() if clusters_dict[i]]
        clusters.sort(key = lambda c:c[0])

        return clusters
    
    def merge_clusters(self, clusters, d=0, close=0.1, max_n=None):
        def remove_bad():
            rms = []
            for i in range(len(clusters)-1,-1,-1):
                cluster = clusters[i]
            #for cluster in clusters[::-1]:
                mean = np.mean(cluster)
                median = np.median(cluster)
                diff = abs(mean - median)
                if diff > close * mean or diff > close * median:
                    if len(cluster) <= 5:
                        print('bad', cluster)
                        rms.append(i)
            for i in rms:
                del clusters[i]

        def force_merge():
            n_merges = len(clusters) - max_n
            seps = []
            for i in range(len(clusters) - 1):
                j = i + 1
                sep = clusters[j][0] - clusters[i][-1]
                seps.append((i, sep))

            seps_sorted = sorted(seps, key=itemgetter(1))
            merges = [(s[0], s[0]+1) for s in seps_sorted[:n_merges]]
            merges.sort(key=itemgetter(0), reverse=True)
            return merges
        
        def merge_close():
            merges = []
            for i in range(len(clusters) - 1):
                j = i + 1
                sep = clusters[j][0] - clusters[i][-1]
                if sep <= d:
                    merges.append((i, j))
            return merges
        
        def merge_singles():
            singles = [i for i in range(len(clusters)) if len(clusters[i]) == 1]
            merges = []
            for i in singles:
                seps = []
                if i > 0:
                    seps.append((i - 1, clusters[i - 1][-1], clusters[i][0] - clusters[i - 1][-1]))
                if i < len(clusters) - 1:
                    seps.append((i + 1, clusters[i + 1][0], clusters[i + 1][0] - clusters[i][0]))
                seps.sort(key=itemgetter(2))
                j, val, sep = seps[0]
                if sep <= close * val and sep <= close * clusters[i][0]:
                    merges.append(sorted((i, j)))
            merges = sorted(list(set([tuple(m) for m in merges])), key=itemgetter(0))
            return merges
             
        def merge(merges):
            merges.sort(key=itemgetter(0), reverse=True)

            for m in merges:
                new_cluster = clusters[m[0]] + clusters[m[1]]
                del clusters[m[1]]
                del clusters[m[0]]
                clusters.append(new_cluster)
                clusters.sort(key=itemgetter(0))
        
        # remove "bad" clusters
        remove_bad()

        merges = []
        if d > 0:
            merges = merge_close()
            if merges:
                merge(merges)

        merges = []
        if close > 0:
            merges = merge_singles()
            if merges:
                merge(merges)

        # remove singletons
        for c in clusters[::-1]:
            if len(c) < 2:
                clusters.remove(c)

        merges = []
        if max_n is not None:
            merges = force_merge()
            if merges:
                merge(merges)

def parse_tsv(tsv):
    loci = defaultdict(list)
    with open(tsv, 'r') as ff:
        for line in ff:
            if line[0] == '#':
                continue
            cols = line.rstrip().split()
            locus = tuple(cols[:4] + [cols[6]])
            status = cols[14]
            if status != 'full':
                continue
            size = int(cols[10])
            gt = cols[13]
            loci[locus].append((size, gt))
    return loci

def format_array(size_gt):
    clusters = []
    last_gt = None
    cluster = []
    gts = sorted(set([float(s[1]) for s in size_gt if s[1] != 'NA']))
    for gt in gts:
        sizes = [int(s[0]) for s in size_gt if s[1] != 'NA' and float(s[1]) == gt]
        clusters.append(sorted(sizes))
    return clusters

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("tsv", type=str, help="straglr tsv output")
    parser.add_argument("--min_cluster_size", type=float, help="minimum number or fraction of supporting reads for allele clustering. Default:0.1", default=0.1)
    parser.add_argument("--max_num_clusters", type=int, help="maximum number of clusters to try. Default:2", default=2)
    parser.add_argument("--max_check_size", type=int, help="threshold of cluster (allele) size (median/mean) to skip filtering by --min_cluster_size. Default:5000", default=5000)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    clu = Cluster(args.min_cluster_size, args.max_num_clusters, args.max_check_size)
    loci = parse_tsv(args.tsv)
    for locus, sizes in loci.items():
        locus_str = '{}:{}-{}'.format(locus[0], locus[1], locus[2])
        old_clusters = format_array(sizes)
        old_means = ['{:.1f}'.format(np.mean(c)) for c in old_clusters]
        old_medians = ['{:.1f}'.format(np.median(c)) for c in old_clusters]
        new_clusters = clu.cluster([s[0] for s in sizes])
        new_means = ['{:.1f}'.format(np.mean(c)) for c in new_clusters]
        new_medians = ['{:.1f}'.format(np.median(c)) for c in new_clusters]
        print(locus_str, len(old_clusters), old_means, old_medians, old_clusters, len(new_clusters), new_means, new_medians, new_clusters, old_clusters==new_clusters)

if __name__ == '__main__':
    main()
