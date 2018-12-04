'''
This file is to analyse the data generated from GetDataForAnalysis
'''


##########################################################
'''
Load the files
'''
import os
os.chdir("/home/leo/Documents/Database/Pipeline/Ready_2_2_0_1")
os.listdir()

import json
with open('training_data_for_DataAnalysis', 'r') as f:
    training = json.load(f)
with open('testing_data_for_DataAnalysis', 'r') as f:
    testing = json.load(f)
    
############################################################
    
           
'''
First define some functions to process the data
'''
def To_seq(aa_sequence):
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    
    TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
                    ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
                    ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]
    
    seq_obj = None

    seq_single_letter = ''
    for aa in aa_sequence:
        for TS in TripleSingle:
            if TS[0] == aa:
                seq_single_letter += TS[1]
    seq_obj = Seq(seq_single_letter, IUPAC.protein)
    
    return seq_obj

#Define a function to calculate the distance matrix using smilarity derived distance
def Distance_matrix(seq_list):
    # Define the scoring rules
    from Bio import Align
    from Bio.SubsMat.MatrixInfo import blosum62
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = blosum62
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1
    aligner.mode = 'global'
 
    import numpy as np
    distance_matrix = np.zeros((len(seq_list), len(seq_list)))
    similarity_matrix = np.zeros((len(seq_list), len(seq_list)))
    l = len(seq_list[0])
    max_similarity = -10000
    min_similarity = 10000
    for i in range(len(seq_list)):
        for j in range(len(seq_list)):
            if seq_list[i] != seq_list[j]:
                similarity_matrix[i, j] =1 - aligner.score(seq_list[i], seq_list[j])/(11 * l)
                if similarity_matrix[i, j] >= max_similarity:
                    max_similarity = similarity_matrix[i, j]
                if similarity_matrix[i, j] <= min_similarity:
                    min_similarity = similarity_matrix[i, j]
            else:
                similarity_matrix[i, j] = 6e6
                
    dim = np.shape(distance_matrix)
    for i in range(dim[0]):
        for j in range(dim[1]):
            if similarity_matrix[i, j] == 6e6:
                distance_matrix[i, j] = 0
            else:
                distance_matrix[i, j] = (max_similarity - similarity_matrix[i,j])
    normalizer = max_similarity - min_similarity
    return distance_matrix/normalizer
'''
Hcluster:
    Perform hierarchical clustering
Inputs: 
    distance_matrix: a symmetric matrix in the np.arrary form, with the main disagnal 0
Ouputs:
    Z: 
        in the form of an array, which could be used to show the 
        dendrogram.
             [[0, 1, 0.5, 2],
              [2, 3, 0.5, 2],
              [4, 5, 1, 4]]
        each row is in the form of [idx1, idx2, dist, sample_count]
        In the above example, there are for samples with index 0, 1, 2, 3.
        4 means the cluster formed in row number 4 - len(samples) = 0
        5 means the cluster formed in row number 5 - len(samples) = 1.
        Therefore, the last row means combine the clusters formed in the above 
        two rows, and the sample count is 4, all the samples.
  
'''

def Hcluster(distance_matrix):
    
    from scipy.cluster.hierarchy import linkage, optimal_leaf_ordering
    from scipy.spatial.distance import squareform
    
    condenced_distance_matrix = squareform(distance_matrix)
    unordered_Z = linkage(condenced_distance_matrix, method = 'complete')
    Z = optimal_leaf_ordering(unordered_Z, condenced_distance_matrix)
    
    return Z
       
'''
Show_elbow:
    To find the best cutting score
Inputs:
    tree: the same form of the output of Hcluster
Outputs:
    a graph showing the relationship between the number of clusters and the 
    cutting score.
'''

def Show_elbow(tree):
    
    from scipy.cluster.hierarchy import cut_tree
    
    n_clusters = []
    heights = [x * 0.01 for x in range(101)]
    for h in heights:
        cut = cut_tree(tree, height = h)
        n = -1
        for i in range(len(cut)):
            if cut[i][0] >= n:
                n = cut[i][0]
        n_clusters.append([n, h])
        
    from matplotlib import pyplot as plt
    
    
    plt.title('clusters Vs cut distance')
    plt.xlabel('cut distance')
    plt.ylabel('n_clusters')
    plt.plot(heights, n_clusters)
    plt.show()
    plt.close()
    
    
    
    return  n_clusters

def Draw_dendrogram(Z, n_clusters = 15):
    from matplotlib import pyplot as plt
    from scipy.cluster.hierarchy import dendrogram
    
    plt.figure(figsize=(10,10))
    dendrogram(
        Z,
        truncate_mode='lastp',
        p=n_clusters,
        leaf_rotation=90.,
        leaf_font_size=12.,
        show_contracted=True,
        show_leaf_counts=True,
    )
    plt.show()
 
'''
Get_cut_clusters_ids:
    To cut the hcluster with given number of clusters, and find the indices in each cluster
Inputs:
    h_cluster: the return of Hcluster
    n_clusters: the number of clusters that the samples are divided into
Outputs:
    clusters: a list in the form of [[0, 1, 2], [3, 4, 6, 7, 8]]
              [0, 1, 2] means samples with indices 0, 1, 2 are in the same cluster.
'''    
def Get_cut_clusters_ids(h_cluster, n_clusters):
    
    from scipy.cluster.hierarchy import cut_tree
    cut1 = cut_tree(Z=h_cluster, n_clusters=n_clusters)
    n = 0
    for i in range(len(cut1)):
        if cut1[i][0] >= n:
            n = cut1[i][0]

    clusters = []
    for i in range(n+1):
        subcluster = []
        for j in range(len(cut1)):
            if cut1[j][0] == i:
                subcluster.append(j)
        clusters.append(subcluster)
    return clusters

###########################################################################
    
'''
Analyse the data
'''
class DataAnalysis(object):

    def __init__(self, training, testing):
        self.training_data = training
        self.testing_data = testing


    '''
    The input data should be in the same form as the loaded training data, or the same 
        as the output of the GetDataForAnalysis
    '''
    def Calculate_distance_matrix(self):

        seq_data__training_Ab = []
        seq_data__training_Ag = [] 
        
        for parepi in self.training_data:
            seq_data__training_Ab.append(To_seq(parepi[0]))
            seq_data__training_Ag.append(To_seq(parepi[1]))
            
        self.training_Ab_distance_martix = Distance_matrix(seq_data__training_Ab)
        self.training_Ag_distance_martix = Distance_matrix(seq_data__training_Ag)

        seq_data__testing_Ab = []
        seq_data__testing_Ag = [] 
        
        for parepi in self.testing_data:
            seq_data__testing_Ab.append(To_seq( parepi[0]))
            seq_data__testing_Ag.append(To_seq(parepi[1]))

        self.testing_Ab_distance_martix = Distance_matrix(seq_data__testing_Ab)
        self.testing_Ag_distance_martix = Distance_matrix(seq_data__testing_Ag)
        
    def Get_hierarchical_clusters(self):
        
        self.training_Ab_hclusters = Hcluster(self.training_Ab_distance_martix)
        self.training_Ag_hclusters = Hcluster(self.training_Ag_distance_martix)
        
    def Get_cluster_ids(self, Ab_n_clusters = None, Ag_n_clusters = None):
        if Ab_n_clusters == None or Ag_n_clusters == None:
            print ('You have to assign the number of clusters you want to cut')
        else:
            self.Ab_cluster_ids = Get_cut_clusters_ids(self.training_Ab_hclusters, Ab_n_clusters)
            self.Ag_cluster_ids = Get_cut_clusters_ids(self.training_Ag_hclusters, Ag_n_clusters)
            
            self.Ab_cluster_ids.sort(key = lambda x:len(x))
            self.Ag_cluster_ids.sort(key = lambda x:len(x))
            
    '''
    self.distance_matrix_Ag_in_Ab is in the form of dictionary with the keys the tuples 
    of the reference cluster ids in the form of (0, 1, 2, 3). All the other dictionaries 
    with the same kind of keys.
    '''
            
    def Cluster_within_cluter(self):
        
        Ab_cluster_ids = self.Ab_cluster_ids
        Ag_cluster_ids = self.Ag_cluster_ids
        
        distance_matrix_Ag_in_Ab = {}
        distance_matrix_Ab_in_Ag = {}
        
        for Ab_cluster in Ab_cluster_ids:
            Ag_in_Ab_cluster_distance_matrix = np.zeros((len(Ab_cluster), len(Ab_cluster)))
            for i in range(len(Ab_cluster)):
                for j in range(len(Ab_cluster)):
                    Ag_in_Ab_cluster_distance_matrix[i, j] = self.training_Ag_distance_martix[Ab_cluster[i], Ab_cluster[j]]
            distance_matrix_Ag_in_Ab[tuple(Ab_cluster)] = Ag_in_Ab_cluster_distance_matrix
            
        for Ag_cluster in Ag_cluster_ids:
            Ab_in_Ag_cluster_distance_matrix = np.zeros((len(Ag_cluster), len(Ag_cluster)))
            for i in range(len(Ag_cluster)):
                for j in range(len(Ag_cluster)):
                    Ab_in_Ag_cluster_distance_matrix[i, j] = self.training_Ab_distance_martix[Ag_cluster[i], Ag_cluster[j]]
            distance_matrix_Ab_in_Ag[tuple(Ag_cluster)] = Ab_in_Ag_cluster_distance_matrix
            
        self.distance_matrix_Ag_in_Ab = distance_matrix_Ag_in_Ab
        self.distance_matrix_Ab_in_Ag = distance_matrix_Ab_in_Ag
        
        hcluster_Ag_in_Ab = {}
        hcluster_Ab_in_Ag = {}
        
        for key in distance_matrix_Ag_in_Ab:
            hcluster_Ag_in_Ab[key] = Hcluster(distance_matrix_Ag_in_Ab[key])
        
        for key in distance_matrix_Ab_in_Ag:
            hcluster_Ab_in_Ag[key] = Hcluster(distance_matrix_Ab_in_Ag[key])
            
        self.hcluster_Ag_in_Ab = hcluster_Ag_in_Ab
        self.hcluster_Ab_in_Ag = hcluster_Ab_in_Ag
        
    def Show_elbow_wthin_clusters(self):
        
        hcluster_Ag_in_Ab = self.hcluster_Ag_in_Ab
        hcluster_Ab_in_Ag = self.hcluster_Ab_in_Ag
        
        n_clusters_Ag_in_Ab = {}

        n_clusters_Ab_in_Ag = {}
        
        key_list_Ab = list(hcluster_Ab_in_Ag.keys())
        for i in range(len(key_list_Ab)):
            n_clusters_Ab_in_Ag[key_list_Ab[i]] = Show_elbow(hcluster_Ab_in_Ag[key_list_Ab[i]])
            print(str(i)+'   hcluster_Ab_in_Ag')
            
        key_list_Ag = list(hcluster_Ag_in_Ab.keys())
        for i in range(len(key_list_Ag)):
            n_clusters_Ag_in_Ab[key_list_Ag[i]] = Show_elbow(hcluster_Ag_in_Ab[key_list_Ag[i]])
            print(str(i)+'   hcluster_Ag_in_Ab')
            
        
                    

import numpy as np
tryout =  DataAnalysis(training, testing)                
tryout.Calculate_distance_matrix()            
tryout.testing_Ab_distance_martix[:5, :5] 
tryout.Get_hierarchical_clusters() 
Ab_clusters = Show_elbow(tryout.training_Ab_hclusters)
Ag_clusters = Show_elbow(tryout.training_Ag_hclusters)
Ab_clusters
Ag_clusters
tryout.Get_cluster_ids(Ab_n_clusters=23, Ag_n_clusters=21)     
len(tryout.Ab_cluster_ids)
len(tryout.Ag_cluster_ids)
len(tryout.Ab_cluster_ids[0])
tryout.Cluster_within_cluter()
type(tryout.distance_matrix_Ab_in_Ag)
len(tryout.distance_matrix_Ab_in_Ag)
type(tryout.hcluster_Ab_in_Ag)
len(tryout.hcluster_Ab_in_Ag)
tryout.Show_elbow_wthin_clusters()

Z_Ab = Hcluster(tryout.training_Ab_distance_martix)
Draw_dendrogram(Z_Ab, n_clusters=23)
clusters = Get_cut_clusters_ids(Z_Ab, n_clusters=23)
for i in clusters:
    print(len(i))
from scipy.spatial.distance import squareform        
from scipy.cluster.hierarchy import cophenet
    
Z_complete = Hcluster(tryout.testing_Ab_distance_martix)  
         
condensed_dm_Ab = squareform(tryout.testing_Ab_distance_martix)
condensed_dm_Ag = squareform(tryout.training_Ag_distance_martix)  

cophenet_dm = cophenet(Z_complete) 

corr_coef = np.corrcoef(condensed_dm_Ab, cophenet_dm)[0, 1]    
corr_coef            
            
            
a = [[1, 2], [3]]
b = a.sort(key = lambda x:len(x)) 
b = a 
b         
TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
                ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
                ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]           




seq = []
for parepi in testing:
    seq.append(To_seq(parepi[0]))            

len(seq)
seq            


            
            
            
            
            
            
            
            
            