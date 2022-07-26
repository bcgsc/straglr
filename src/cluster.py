import numpy as np
from numpy.random import rand
from sklearn.mixture import GaussianMixture
from scipy import stats

class Cluster:
    def __init__(self, min_pts, max_num_clusters):
        self.min_pts = min_pts
        self.max_num_clusters = max_num_clusters
        
    def cluster(self, data):
        try:
            return self.gmm(data)
        except:
            return []
        
    def gmm(self, x):
        X = np.array([[i] for i in x])
        N = np.arange(1, self.max_num_clusters + 1)
        models = [None for i in range(len(N))]
        AICs = []
        BICs = []
        if len(set(x)) == 1:
            return [x]
        for i in range(len(N)):
            models[i] = GaussianMixture(N[i], init_params="k-means++").fit(X)
        AIC = [m.aic(X) for m in models]
        BIC = [m.bic(X) for m in models]
        M_best = models[np.argmin(AIC)]
        labels = M_best.predict(X)
        
        nclusters = M_best.n_components
        clusters = [[] for i in range(M_best.n_components)]
        for label, i in zip(labels, X):
            clusters[label].append(i[0])
            
        # sort clusters
        for cluster in clusters:
            cluster.sort()
        clusters.sort(key = lambda c:c[0])
        for cluster in clusters:
            m = np.mean(cluster)
            s = np.std(cluster)

        merged_clusters = self.merge_gmm_clusters(clusters)

        return [c for c in merged_clusters if len(c) >= self.min_pts]

    def merge_gmm_clusters(self, clusters, threshold=50, rel_threshold=0.05):
        if len(clusters) <= 1:
            return clusters

        bounds = {}
        means = {}
        for i in range(len(clusters)):
            m = np.mean(clusters[i])
            s = np.std(clusters[i])
            bounds[i] = m - 2*s, m + 2*s
            means[i] = m
        for i in range(len(clusters)-1):
            j=i+1

        i = 0
        merged = []
        while i < len(clusters) - 1:
            merged.append([i])
            j = i + 1
            while j < len(clusters):
                i = j
                if means[j] - means[merged[-1][-1]] < max(threshold, rel_threshold * min(means[j], means[merged[-1][-1]])):
                    merged[-1].append(j)
                    j += 1
                else:
                    if j == len(clusters) - 1:
                        merged.append([j])
                    break

        new_clusters = []
        for i in merged:
            cc = []
            for j in i:
                cc += clusters[j]
            new_clusters.append(cc)

        return new_clusters
