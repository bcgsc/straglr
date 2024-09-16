import numpy as np
from numpy.random import rand
from sklearn.mixture import GaussianMixture
from scipy import stats
from operator import itemgetter
from collections import defaultdict

class Cluster:
    def __init__(self, min_pts, max_num_clusters=2):
        self.min_pts = min_pts
        self.max_num_clusters = max_num_clusters

    def cluster(self, data):
        if len(data) == 1:
            return [[data[0]]]
        clusters = self.gmm(data, self.min_pts, self.max_num_clusters + 4)

        self.merge_clusters(clusters, d=10, max_n=self.max_num_clusters)

        clusters = [c for c in clusters if len(c) >= self.min_pts]
        return clusters

    def gmm(self, x, min_pts, max_num_clusters):
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
            models[i] = GaussianMixture(N[i], init_params="k-means++", random_state=0).fit(X)
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
    
    def merge_clusters(self, clusters, d=0, max_n=None):
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
        
        def merge(merges):
            merges.sort(key=itemgetter(0), reverse=True)

            for m in merges:
                new_cluster = clusters[m[0]] + clusters[m[1]]
                del clusters[m[1]]
                del clusters[m[0]]
                clusters.append(new_cluster)
                clusters.sort(key=itemgetter(0))
        
        merges = []
        if d > 0:
            merges = merge_close()
            if merges:
                merge(merges)

        #clusters = [c for c in clusters if len(c) >= self.min_pts]
        # remove singletons
        for c in clusters[::-1]:
            if len(c) < 2:
                clusters.remove(c)

        merges = []
        if max_n is not None:
            merges = force_merge()
            if merges:
                merge(merges)
