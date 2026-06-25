import pandas as pd
import numpy as np
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm
import scipy.signal
import scipy.stats as stats
from scipy.spatial.distance import pdist
from scipy.spatial.distance import cdist
from matplotlib import pyplot as plt

import sys
import os

class cluster:

    def __init__(self, data, distance_metric, gravity=1):
        #load the data, assume pd dataframe
        self.name=''
        self.index=0
        self.data=data
        self.center=self.data.mean(axis=0)
        self.distance_metric=distance_metric
        self.gravity=gravity
        #reserve an aray to mark draft clusters
        self.sidenote=pd.Series(data=np.zeros(data.shape[0]), index=data.index)
        self.steps=0
        self.eminense=0
        self.eminence_metric=["density","mean_distance_inside","steps","mean_nearest", "stdev"]

    def gravity_center(self, scoop):
        #this function defines the relative contribution of old vs new cluster elements
        #in the position of the formal element, the new starting point for the next iteration
        #high gravity makes the clustering conservative (default 1 makes it as conservative as kmeans)
        #lower gravity tends to chaise long traces of objects in multidimensional space
        if(self.gravity==1):
            center = self.data.mean(axis=0)
        else:
            #here the contribution of the cluster members from the previous steps are lighter by the gravity factor
            center*=self.gravity
            center=(scoop.mean(axis=0)+self.center)/2
        return center

    #make a single trip with the hypersphere, given the starting point and the radius
    #estimate and return the eminence of the draft cluster
    def walk(self, start, radius):
        #print("walking from:",start)
        #print("radius:",radius)
        #print("internal:", self.data.index.values)
        center = self.data.loc[start].values.reshape(1, -1)
        self.steps=0
        covered=np.array([])
        while(True):
            self.steps+=1
            #find distance to all data points
            dist=cdist(self.data.values, center, metric=self.distance_metric)
            #select points within the radius
            scoop=self.data[dist<=radius]
            #mark all selected points
            covered=np.append(covered,scoop.index.values)
            #find the gravity center and the eminence of the provisional cluster
            new_center=self.gravity_center(scoop).values.reshape(1, -1)
            trip=cdist(center, new_center, metric=self.distance_metric)[0][0]
            if(trip == 0):
                break
            else:
                center=new_center
        #all objects touched by the hypersphere on the trip are marked on a side
        #and eventually returned as a draft cluster
        covered=np.unique(covered)
        self.data=self.data.loc[covered]
        return(self)

    #weight up how good the draft cluster is by one of the many possible metrics and return the
    #estimated eminence
    def estimate_eminence(self, metric='density'):
        if metric=='density':
            #estimate cluster radius as the distance from center to the most distant element
            #then estimate density and the number of elements within this radius
            self.center=self.data.mean(axis=0)
            center2d = self.center.values.reshape(1, -1)
            dist = cdist(self.data.values, center2d, metric=self.distance_metric)
            rmax=dist.max()
            density=rmax/self.data.shape[0]
            self.eminense=density
            return(self.eminense)
        elif metric=='mean_distance_inside':
            #calculate the distance from each point to the center, take average, turn to 1/n
            self.center=self.data.mean(axis=0)
            center2d = self.center.values.reshape(1, -1)
            dist = cdist(self.data.values, center2d, metric=self.distance_metric)
            mean_distance_inside=dist.mean()
            #just to make it the bigger, the better in relation to the most different pair of samples
            self.eminense=dist.max()-mean_distance_inside
            return(self.eminense)
        elif metric=='steps':
            #simple metric of how many steps did it take to cover all elements in a walk
            self.eminense=self.steps
            return(self.eminense)
        elif metric=='mean_nearest':
            #average distance between nearest neighbors
            self.eminense=0
            for i in range(self.data.shape[0]):
                for j in range(i+1,self.data.shape[0]):
                    d=euclidean_distances(self.data[i],self.data[j])
                    if(d>self.eminense):
                        self.eminense=d
            return(self.eminense)
        elif metric=='stdev':
            #use stdev as a metric
            self.eminense=self.data.std(axis=0)
            return(self.eminense)

#perform ForEl clustering
def ForEl(data, low_margin=0.1, high_margin=0.1, num_increments=100, min_cluster=2, metric='mean_distance_inside', distance_metric='euclidean'):
    global cluster
    best_eminense=0
    best_cluster=None
    cluster_index=0
    cluster_list=[]
    #test all possible values of a radius in a range between the minimal distance and the maximal distance
    #within the data set with a given number of increments (default: 100). Cut the margins from the lowers and highest values
    distances=pdist(data,distance_metric)
    min_val = distances[distances != 0].min()
    max_val = distances.max()
    #minimal radius - min distance, max radius - max distance; try every possible radius value in this range with a given increment
    increment=(max_val-min_val)/num_increments
    #reduce the range by cutting some margins from the extreeme values, default cut 10% from each side
    high=max_val-high_margin*(max_val-min_val)
    low=min_val+low_margin*(max_val-min_val)
    #make a copy of the original data to operate on
    op_data=data.copy()
    #find the best cluster, extract, iterate until all objects in data are classified
    while (1):
        print("iteration:",cluster_index+1)
        all_starts = op_data.index
        best_eminense=0
        best_cluster=None
        for start in tqdm(all_starts):
            r=low
            while(r<=high):
                clu = cluster(op_data, distance_metric)
                clu.walk(start, r)
                clu.estimate_eminence(metric)
                if(clu.eminense > best_eminense):
                    best_eminense=clu.eminense
                    best_cluster=clu
                r+=increment
        #assign the best cluster its index
        best_cluster.index=cluster_index
        #increment the index for the next iteration
        cluster_index+=1
        #add the best cluster to the list
        cluster_list.append(best_cluster)
        #delete the best cluster from the data set
        op_data = op_data.loc[np.setxor1d(best_cluster.data.index.values,op_data.index.values)]
        #check the exit conditions
        if(op_data.shape[0]<min_cluster):
            #leftover items are too few
            break
        #additional exit conditions may include cluster eminence deteriorating below the cutoff value
    #assemble the resulting dataframe for output
    cluster_index=0
    data['cluster']=-1
    for cluster in cluster_list:
        #for each cluster, assign cluster index to all corresponding members in the initial dataframe
        data.loc[cluster.data.index,'cluster'] = cluster_index
        cluster_index+=1
    #add the remaining unclustered operational data as singletons
    for singleton in op_data.index:
        data.loc[singleton, 'cluster'] = cluster_index
        cluster_index+=1
    return data

#df=pd.read_csv('liver_fpkm.csv', index_col='gene_id').fillna(0)

# Load the Iris dataset
from sklearn.datasets import load_iris
iris = load_iris()
# Create the DataFrame
irisdf = pd.DataFrame(data=iris.data, columns=iris.feature_names)
irisdf['species'] = pd.Categorical.from_codes(iris.target, iris.target_names)
cluster_labels=irisdf.pop('species')
#irisdf = irisdf.iloc[:50]
#perform the ForEl clustering
new_clustered_df=ForEl(irisdf,metric='mean_distance_inside')
print(new_clustered_df.all)

