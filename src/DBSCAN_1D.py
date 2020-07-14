import numpy as np
from numpy.random import rand
from collections import Counter

#MY_DBSCAN : a custom DBSCAN algorithm implementation for 1D values only
class DBSCAN_1D:
    # 1D distance
    def distance_1D(self, num1, num2):
        return abs(num1 - num2)
        
    def regionQuery(self, P):
        neighbourPts = []
        for point in self.D:
            if point not in self.visited:
                if self.distance_1D(P,point)<self.eps:
                    neighbourPts.append(point)

        return neighbourPts

    def inAnyCluster(self, point):
        for cluster in self.C:
            if point in cluster:
                return True
        return False
            
    def expandCluster(self, P, neighbourPts):
        # first append the current point to this new cluster
        # self.C[self.c_n].append(P)
        # for each of the points in the neighbourhood

        for point in neighbourPts:
            if point not in self.visited:
                self.visited.append(point)
                neighbourPts_2 = self.regionQuery(point)
                
                #if len(neighbourPts_2) >= self.MinPts:
                # adds all the neighbours to the list of neighbours
                # this includes previous points already in the list potentially, but we don't care
                # as those will be filtered by the visited list
                neighbourPts += neighbourPts_2 
                # adds the point to the cluster if not in any cluster yet
                if not self.inAnyCluster(point):
                    self.C[self.c_n].append(point)
                    
    def add_duplicates(self, D):
        counts = Counter(D)
        for point, count in counts.iteritems():
            if count > 1:
                for cluster in self.C:
                    if point in cluster:
                        for i in range(2, count + 1):
                            cluster.append(point)
                            
                for n in self.noise:
                    if point in n:
                        for i in range(2, count + 1):
                            n.append(point)
                            
    def sort_clusters(self):
        for cluster in self.C:
            cluster.sort()
        self.C.sort(key = lambda c:c[0])    
    
    def __init__(self, D, eps=100, MinPts=2):
        self.D = D
        self.eps = eps
        self.MinPts = MinPts
        self.noise = []
        self.visited = []
        self.C = []
        self.c_n = -1
        
        # run through all the points in the data
        for point in D:
            self.visited.append(point) #marking point as visited
            
            # gets all the neighbouring points within the distance defined by eps
            neighbourPts = self.regionQuery(point)

            if not self.inAnyCluster(point):
                self.C.append([])
                self.c_n+=1
                self.C[self.c_n].append(point)
                # see if we can expand the cluster further by adding points 
                # until there is a gap of eps
                self.expandCluster(point, neighbourPts)
            # point was completely expanded and cluster is complete
            
            # if the length of the cluster is not long enough discard it
            if len(self.C[self.c_n]) < (self.MinPts):
                self.noise.append(self.C[self.c_n])
                del(self.C[self.c_n])
                self.c_n-=1

        # add back duplicate numbers
        #self.add_duplicates(D)
        
        # sort lists
        self.sort_clusters()        