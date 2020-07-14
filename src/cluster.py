import numpy as np
from numpy.random import rand
from sklearn.mixture import GaussianMixture
from .DBSCAN_1D import DBSCAN_1D
from scipy import stats

class Cluster:
    def __init__(self, method, min_pts, max_num_clusters, eps=50):
        self.method = method
        self.min_pts = min_pts
        self.max_num_clusters = max_num_clusters
        self.eps = eps
        
    def cluster(self, data):
        try:
            if self.method == 'gmm':
                return self.gmm(data)

            elif self.method == 'dbscan':
                # make duplicates by adding 0.1 (assuming inputs are sizes in integers)
                seen = {}
                data_cluster = []
                for d in data:
                    if not d in seen.keys():
                        data_cluster.append(d)
                        seen[d] = d
                    else:
                        seen[d] += 0.1
                        data_cluster.append(seen[d])

                # set eps
                #percent1 = float(self.min_pts) * 100 / len(data)
                #percent2 = 100.0 - percent1
                #iqr = stats.iqr(data, rng=(percent1, percent2))
                ##iqr = stats.iqr(data, rng=(10,90))
                #if iqr < 100:
                    #eps = 20
                #elif iqr < 200:
                    #eps = 50
                #else:
                    #eps = 100
                #eps = max(20, stats.iqr(data)/4)
                #eps = 100

                iqr = stats.iqr(data, rng=(10, 90))
                eps = min(500, max(20, iqr/4))

                # reverting them back to integers
                cluster_ints = []
                for cluster in self.dbscan(data_cluster, eps):
                    cluster_int = []
                    for d in cluster:
                        cluster_int.append(int(d))
                    cluster_ints.append(cluster_int)
                return cluster_ints
                #return self.dbscan(data)
                #dbscan = DBSCAN_1D(data, self.eps, self.min_pts)
                #return dbscan.C
        except:
            return []
        
    def gmm(self, x):
        X = np.array([[i] for i in x])
        #N = np.arange(1, 4)
        N = np.arange(1, self.max_num_clusters + 1)
        models = [None for i in range(len(N))]
        AICs = []
        BICs = []
        for i in range(len(N)):
            models[i] = GaussianMixture(N[i]).fit(X)
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
            #print('nn', i, clusters[i], m, s, bounds[i])
        for i in range(len(clusters)-1):
            j=i+1
            #print('oo', i, j,  bounds[i], bounds[j], clusters[i], clusters[j], bounds[i][1] > bounds[j][0])

        i = 0
        merged = []
        while i < len(clusters) - 1:
            merged.append([i])
            j = i + 1
            while j < len(clusters):
                i = j
                if means[j] - means[merged[-1][-1]] < max(threshold, rel_threshold * min(means[j], means[merged[-1][-1]])):
                #if means[j] - means[merged[-1][-1]] < threshold:
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
        
    def dbscan(self, x, eps=None):
        if eps is None:
            eps = self.eps
        dbscan_1d = DBSCAN_1D(x, eps, 2)
        
        return [c for c in dbscan_1d.C if len(c) >= self.min_pts]
        
