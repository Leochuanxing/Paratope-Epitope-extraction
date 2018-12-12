'''
This file is to analyse the data generated from GetDataForAnalysis
'''
import random
import numpy as np
import os
import json
import math
import copy
##########################################################
'''
Load the files
'''

os.chdir("/home/leo/Documents/Database/Pipeline/Ready_2_2_0_1")
os.listdir()


with open('training_2_2_1_1', 'r') as f:
    training = json.load(f)
with open('testing_2_2_1_1', 'r') as f:
    testing = json.load(f)
with open('negative_centers_5', 'r') as f:
    negative_samples = json.load(f)

len(training)    
len(testing)
testing[:18]
training[:6]
negative_samples[:6]
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
    '''
    Distance_matrix:
        to calculate the distance matrix
    Input:
        seq_list: a list with elements seq_object
    Output:
        distance, a distance matrix given by the method definde in  the function
    '''
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
    l = len(seq_list[0])
#    max_score = -1000
#    min_score = 1000
    for i in range(len(seq_list)):
        for j in range(len(seq_list)):
            if seq_list[i] != seq_list[j]:
                distance_matrix[i, j] =1 - (4*l + aligner.score(seq_list[i], seq_list[j]))/(15*l)#Global
#                 distance_matrix[i, j] =1 -  aligner.score(seq_list[i], seq_list[j])/(15*l)#Local
#                distance_matrix[i, j] = - aligner.score(seq_list[i], seq_list[j])#Unnormalized
            else:
                distance_matrix[i, j] = 0
                
   # fliped_dm = np.zeros(np.shape(distance_matrix))
    return distance_matrix
# check the diatance mat
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

#from scipy.spatial.distance import squareform  
#condenced = squareform(tryout.all_Ab_distance_martix)
#len(condenced)
#np.shape(tryout.all_Ag_distance_martix)
#len(tryout.Ab_cluster_ids)
#len(tryout.Ag_cluster_ids)     
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
    
    n_clusters_heights =[]
    n_clusters = []
    heights = [x * 0.01 for x in range(101)]
    for h in heights:
        cut = cut_tree(tree, height = h)
        n = -1
        for i in range(len(cut)):
            if cut[i][0] >= n:
                n = cut[i][0]
        n_clusters.append(n)
        n_clusters_heights.append([n, h])
        
#    from matplotlib import pyplot as plt
#    
#    
#    plt.title('clusters Vs cut distance')
#    plt.xlabel('cut distance')
#    plt.ylabel('n_clusters')
#    plt.plot(heights, n_clusters)
#    plt.show()
#    plt.close()   
    
    return  n_clusters_heights


'''
Find_the_cut_height:
    for a given cluster, which is the returned value of Show_elbow, find the number of clusters
    corresponding to the steepest drop, and the number of clusters should be in the range of 
    [lower_n_cluster, upper_n_cluster]. If the total number of clusters is smaller than the 
    lower_n_cluster, then the number of clusters with cut height 0.01 will be returned
Inuput:
    clusters: the return value of Show_elbow
    lower_n_cluster: the lower limit of the returned cluster number
    upper_n_cluster: the upper limit of the returned cluster number
Output:
    cut_cluster: The number of clusters after the steepest drop within the given range.
    
'''  
def Find_the_cut_n_clusters(n_clusters, lower_n_cluster, upper_n_cluster):
    
    n_cluster = []
    for i in n_clusters:
        n_cluster.append(i[0])
    diff = np.diff(n_cluster)
    
    steepest = 10000
    cut_cluster = -1
    
    if n_cluster[0] <= lower_n_cluster:
        cut_cluster = n_clusters[1][0]
    else:
        for j in range(len(diff)):
            if n_cluster[j+1] >= lower_n_cluster and n_cluster[j+1] <= upper_n_cluster:
                if diff[j] <= steepest:
                    steepest = diff[j]
                    cut_cluster = n_clusters[j+1][0]
    return cut_cluster
            

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
#from scipy.cluster.hierarchy import cut_tree
#from scipy.cluster.hierarchy import linkage

#z = linkage(condenced, method='complete')
#cut1 = cut_tree(z, n_clusters=50)
#len(cut1)
#cut1[10:26]
'''
The returned value is the cluster index
Input: 
    clusters: a lists of clusters in the form like [[1, 2], [3, 4]], where the numbers 
    are the samples
    element: an integer, gives the index of the sample
    top_n: int, gives the number of clusters the element is assigned to.
    distance_matrix: gives the distance between paires of samples
    
    All the index of the samples should refer to the original data list
    
Output: a list of indices of the the clusters. For example, if top_n = 2, the returned 
     values maybe in the form of [1, 2]. This means the element should be assigned to cluster
     clusters[1], and clusters[2].
     
'''

def Assign_to_cluster(clusters, element, distance_matrix, mode = 'complete', top_n = 1):
    
    distances_complete = []
    distances_single = []
    distances_average = []
    
    for i in range(len(clusters)):
        max_dist = -1000
        min_dist = 1000
        average = 0
        for aa in clusters[i]:
            if distance_matrix[aa, element] >= max_dist:
                max_dist = distance_matrix[aa, element]
            if distance_matrix[aa, element] <= min_dist:
                min_dist = distance_matrix[aa, element]
            average += distance_matrix[aa, element]
            
        distances_complete.append([i, max_dist])
        distances_single.append([i, min_dist])
        distances_average.append([i, average/len(clusters[i])])
        
    distances_complete.sort(key = lambda x : x[1])
    distances_single.sort(key = lambda x : x[1])
    distances_average.sort(key = lambda x : x[1] )
    
    
    assigned_cluster = []
    if mode == 'complete':
        for i in range(top_n):
            assigned_cluster.append(distances_complete[i][0])
    if mode == 'single':
        for i in range(top_n):
            assigned_cluster.append(distances_single[i][0])
    if mode == 'average':
        for i in range(top_n):
            assigned_cluster.append(distances_average[i][0])
            
    return assigned_cluster

#a = [[1, 2], [2, 1]]
#a.sort(key = lambda x:x[1])
#a
'''
Prediction:
     For a given x, predict y, if the prediction is correct, the returned value is 1, otherwise, the returned value is 0
Inputs:
    parameter, a dictionary, should contain the following values
    x_trained cluster, this is in the form of [[1,2], [3,4]], where 1,2,3,4 corresponding to the original data
    y_trained_clustter_dict: a disctionary with keys the tuples of the elements in x_trained_clusters
    x_testing, an integer, corresponding to the numbers of the original data
    y_testing, an integer, corresponding to the numbers of the original data, usually x_testing = y_testing
    top_x_cluster, an integer, gives the nearest 'top_x_cluster' x_testing is assigned to
    top_y_cluser, the number of predictions given
    mode: a string, takes values of either 'complete' or 'single'. If it is complete, the maximum distance is taken as the 
        distance between the cluster and the element. If it is 'single', the smallest distance is takeen as the distance 
        between the cluster and the element.
Output:
    successful_prediction, an integer. If the prediction is successful, the returned value is 1. Otherwise it is 0.
    
     The clusters in y_trained_cluster_dict must be arranged in decreasing order
'''
def Prediction(parameter):
    # take out the parameters
    x_trained_cluster = parameter['x_trained_cluster']
    y_trained_cluster_dict = parameter['y_trained_cluster']
    x_testing = parameter['x_testing']
    y_testing = parameter['y_testing']
    top_x_n_cluster = parameter['top_x_n_cluster']
    top_y_n_cluster = parameter['top_y_n_cluster']
    x_distance_matrix = parameter['x_distance_matrix']
    y_distance_matrix = parameter['y_distance_matrix']
    mode = parameter['mode']
    # 
    successful_prediction = 0
    x_assign_to_cluster = Assign_to_cluster(x_trained_cluster, x_testing, x_distance_matrix,  mode, top_x_n_cluster)
    for i in x_assign_to_cluster:
        y_assign_to_cluster = Assign_to_cluster(y_trained_cluster_dict[tuple(x_trained_cluster[i])], y_testing,
                                                                       y_distance_matrix, mode, 1)
        if y_assign_to_cluster[0] <= top_y_n_cluster - 1:
            successful_prediction = 1
            break
        print(y_assign_to_cluster)
            
    return successful_prediction

def Meaningfull_prediction_support(parameter):
    train_l = parameter['len_training_data']
    test_l = parameter['test_testing_data']
    Ag_cluster_ids = parameter['Ag_cluster_ids']
    Ab_cluster_ids = parameter['Ab_cluster_ids']
    n_prediction = parameter['n_prediction']
    all_Ab_distance_martix = parameter['all_Ab_distance_martix']
    all_Ag_distance_martix = parameter['all_Ag_distance_martix']
    support_Ab = parameter['support_Ab']
    support_Ag = parameter['support_Ag']
    mode = parameter['mode']
    x_method = parameter['x_method']
    y_method = parameter['y_method']
    
    if mode is 'Ag_to_Ab':
        success = 0
        for i in range(train_l, train_l+test_l):
            Ag_prediction = Assign_to_cluster(Ag_cluster_ids, i,
                                              all_Ag_distance_martix, x_method, top_n = 1)
            Ab_prediction = Assign_to_cluster(support_Ab, i, 
                                              all_Ab_distance_martix, y_method, n_prediction )
            if Ag_prediction[0] in Ab_prediction:
                success += 1                
        
    if mode is 'Ab_to_Ag':
        success = 0
        for i in range(train_l, train_l+test_l):
            Ab_prediction = Assign_to_cluster(Ab_cluster_ids, i,
                                              all_Ab_distance_martix, x_method, top_n = 1)
            Ag_prediction = Assign_to_cluster(support_Ag, i, 
                                              all_Ag_distance_martix, y_method, n_prediction )
            if Ab_prediction[0] in Ag_prediction:
                success += 1    
    
    return round(success/test_l, 2)
    
##############################################################################################
###############################################################################################
    
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

        seq_data__all_Ab = seq_data__training_Ab[:]
        seq_data__all_Ag = seq_data__training_Ag[:]
        
        for parepi in self.testing_data:
            seq_data__all_Ab.append(To_seq( parepi[0]))
            seq_data__all_Ag.append(To_seq(parepi[1]))
        # the following matrix will be used in the prediction
        self.all_Ab_distance_martix = Distance_matrix(seq_data__all_Ab)
        self.all_Ag_distance_martix = Distance_matrix(seq_data__all_Ag)
        
    def Get_hierarchical_clusters(self):
        
        self.training_Ab_hclusters = Hcluster(self.training_Ab_distance_martix)
        self.training_Ag_hclusters = Hcluster(self.training_Ag_distance_martix)
        
    def Show_elbow(self):
        self.Ab_height_n_clusters = Show_elbow(self.training_Ab_hclusters)
        self.Ag_height_n_clusters = Show_elbow(self.training_Ag_hclusters)

        
    def Cut_n_clusters(self, lower_n_cluster, upper_n_cluster):
        
        self.cut_Ab_n_clusters = Find_the_cut_n_clusters(self.Ab_height_n_clusters, 
                                                         lower_n_cluster=lower_n_cluster, upper_n_cluster= upper_n_cluster)

        self.cut_Ag_n_clusters = Find_the_cut_n_clusters(self.Ag_height_n_clusters, 
                                                         lower_n_cluster=lower_n_cluster, upper_n_cluster= upper_n_cluster)        


        
    def Get_cluster_ids(self):
        self.Ab_cluster_ids = Get_cut_clusters_ids(self.training_Ab_hclusters, self.cut_Ab_n_clusters)
        self.Ag_cluster_ids = Get_cut_clusters_ids(self.training_Ag_hclusters, self.cut_Ag_n_clusters)
        
        self.Ab_cluster_ids.sort(key = lambda x:len(x))
        self.Ag_cluster_ids.sort(key = lambda x:len(x))
            
    '''
    self.distance_matrix_Ag_in_Ab is in the form of dictionary with the keys the tuples 
    of the reference cluster ids in the form of (0, 1, 2, 3). All the other dictionaries 
    with the same kind of keys.
    '''
            
    def Cluster_within_cluter(self, min_cluster_number, max_cluster_number):
        
        Ab_cluster_ids = copy.deepcopy(self.Ab_cluster_ids)
        Ag_cluster_ids = copy.deepcopy(self.Ag_cluster_ids)
        # select the clusters with one element to avoid error in the following clustering process
        Ab_cluster_id_single = []
        Ag_cluster_id_single = []
        
        for Ab in Ab_cluster_ids:
            if len(Ab) <= 1:
               Ab_cluster_id_single.append(Ab)
        for Ag in Ag_cluster_ids:
            if len(Ag) <= 1:
                Ag_cluster_id_single.append(Ag)
                
        for Ab in Ab_cluster_id_single:
            Ab_cluster_ids.remove(Ab)
        for Ag in Ag_cluster_id_single:
            Ag_cluster_ids.remove(Ag)
       
        # calculate distance matrix within each cluster
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
            

        # hierarchical clustering within each cluster
        hcluster_Ag_in_Ab = {}
        hcluster_Ab_in_Ag = {}
        
        for key in distance_matrix_Ag_in_Ab:
            hcluster_Ag_in_Ab[key] = Hcluster(distance_matrix_Ag_in_Ab[key])
        
        for key in distance_matrix_Ab_in_Ag:
            hcluster_Ab_in_Ag[key] = Hcluster(distance_matrix_Ab_in_Ag[key])
            
        # calculate the cluster numbers with differenct cut height
        n_clusters_heights_Ag_in_Ab = {}
        n_clusters_heights_Ab_in_Ag = {}
        
        key_list_Ab = list(hcluster_Ab_in_Ag.keys())
        for i in range(len(key_list_Ab)):
            n_clusters_heights_Ab_in_Ag[key_list_Ab[i]] = Show_elbow(hcluster_Ab_in_Ag[key_list_Ab[i]])

        
        key_list_Ag = list(hcluster_Ag_in_Ab.keys())
        for i in range(len(key_list_Ag)):
            n_clusters_heights_Ag_in_Ab[key_list_Ag[i]] = Show_elbow(hcluster_Ag_in_Ab[key_list_Ag[i]])

            
        # find proper number of cluster to cluser
        Cut_n_clusters_Ab_in_Ag = {}
        for key in n_clusters_heights_Ab_in_Ag:
            Cut_n_clusters_Ab_in_Ag[key] = Find_the_cut_n_clusters(n_clusters_heights_Ab_in_Ag[key],
                                   lower_n_cluster= min_cluster_number, upper_n_cluster=max_cluster_number)
        
        
        Cut_n_clusters_Ag_in_Ab = {}
        for key in n_clusters_heights_Ag_in_Ab:
            Cut_n_clusters_Ag_in_Ab[key] = Find_the_cut_n_clusters(n_clusters_heights_Ag_in_Ab[key], 
                                   lower_n_cluster=min_cluster_number, upper_n_cluster=max_cluster_number)

        
        # get the with in cluster ids
        
        Ab_in_Ag_cluster_ids = {}
        for key in hcluster_Ab_in_Ag:
            indirect_ids = Get_cut_clusters_ids(hcluster_Ab_in_Ag[key], Cut_n_clusters_Ab_in_Ag[key])
            direct_group_ids = []
            for group in indirect_ids:#in the form of [[1, 2], [3,4,5]]
                direct_ids = []
                for idx in group:
                    direct_ids.append(key[idx])
                direct_group_ids.append(direct_ids)
            direct_group_ids.sort(key = lambda x:len(x), reverse = True)
            Ab_in_Ag_cluster_ids[key] = direct_group_ids
        self.Ab_in_Ag_cluster_ids = Ab_in_Ag_cluster_ids
            
        Ag_in_Ab_cluster_ids = {}
        for key in hcluster_Ag_in_Ab:
            indirect_ids = Get_cut_clusters_ids(hcluster_Ag_in_Ab[key], Cut_n_clusters_Ag_in_Ab[key])
            direct_group_ids = []
            for group in indirect_ids:#in the form of [[1, 2], [3,4,5]]
                direct_ids = []
                for idx in group:
                    direct_ids.append(key[idx])
                direct_group_ids.append(direct_ids)
            direct_group_ids.sort(key = lambda x:len(x), reverse = True)
            Ag_in_Ab_cluster_ids[key] = direct_group_ids
        self.Ag_in_Ab_cluster_ids = Ag_in_Ab_cluster_ids
        
        # add the single ones
        if Ab_cluster_id_single != []:
            for single in Ab_cluster_id_single:
                self.Ag_in_Ab_cluster_ids[tuple(single)] = [single]
                
        if Ag_cluster_id_single != []:
            for single in Ag_cluster_id_single:
                self.Ab_in_Ag_cluster_ids[tuple(single)] = [single]
                    

        
    def Ag_to_Ab_prediction_correct_rate(self, top_x_n_cluster, top_y_n_cluster):
        parameter = {}
        parameter['x_trained_cluster'] = self.Ag_cluster_ids
        parameter['y_trained_cluster'] = self.Ab_in_Ag_cluster_ids
        parameter['top_x_n_cluster'] = top_x_n_cluster
        parameter['top_y_n_cluster'] = top_y_n_cluster
        parameter['x_distance_matrix'] = self.all_Ag_distance_martix
        parameter['y_distance_matrix'] = self.all_Ab_distance_martix
        parameter['mode'] = 'single'
        n = 0
        for i in range(len(self.training_data), len(self.training_data) + len(self.testing_data)):
            parameter['x_testing'] = i
            parameter['y_testing'] = i
            n += Prediction(parameter)
            
        return round(n/len(self.testing_data), 2) 
            
    def Ab_to_Ag_prediction_correct_rate(self, top_x_n_cluster, top_y_n_cluster):
        parameter = {}
        parameter['x_trained_cluster'] = self.Ab_cluster_ids
        parameter['y_trained_cluster'] = self.Ag_in_Ab_cluster_ids
        parameter['top_x_n_cluster'] = top_x_n_cluster
        parameter['top_y_n_cluster'] = top_y_n_cluster
        parameter['x_distance_matrix'] = self.all_Ab_distance_martix
        parameter['y_distance_matrix'] = self.all_Ag_distance_martix
        parameter['mode'] = 'single'
        n = 0
        for i in range(len(self.training_data), len(self.training_data) + len(self.testing_data)):
            parameter['x_testing'] = i
            parameter['y_testing'] = i
            n += Prediction(parameter)
            
        return round(n/len(self.testing_data), 2) 
    
    def Inter_intra_cluster_distance():
        pass   

    

    
    def Selected_support_predcition(self, top_n_clusters_support=3, n_prediction = 3, mode = 'Ag_to_Ab',
                                    x_method = 'complete', y_method = 'average'):

        support_Ag_dict = {}
        for key, value in self.Ag_in_Ab_cluster_ids.items():
            support_Ag_dict[key]  = []
            for i in range(min(top_n_clusters_support, len(value))):
                support_Ag_dict[key].extend(value[i])
                
        support_Ab_dict = {}
        for key, value in self.Ab_in_Ag_cluster_ids.items():
            support_Ab_dict[key]  = []
            for i in range(min(top_n_clusters_support, len(value))):
                support_Ab_dict[key].extend(value[i])
                
        support_Ag = [0]*len(self.Ab_cluster_ids)
        for i in range(len(support_Ag)):
            for key, value in support_Ag_dict.items():
                if key == tuple(self.Ab_cluster_ids[i]):
                    support_Ag[i] = value
                    
        support_Ab = [0]*len(self.Ag_cluster_ids)
        for i in range(len(support_Ab)):
            for key, value in support_Ab_dict.items():
                if key == tuple(self.Ag_cluster_ids[i]):
                    support_Ab[i] = value
                    
        parameter = {}
        parameter['len_training_data'] = len(self.training_data)
        parameter['test_testing_data'] = len(self.testing_data)
        parameter['Ag_cluster_ids'] = self.Ag_cluster_ids
        parameter['Ab_cluster_ids'] = self.Ab_cluster_ids
        parameter['n_prediction'] = n_prediction
        parameter['all_Ab_distance_martix'] = self.all_Ab_distance_martix
        parameter['all_Ag_distance_martix'] = self.all_Ag_distance_martix
        parameter['support_Ab'] = support_Ab
        parameter['support_Ag'] = support_Ag
        parameter['mode'] = mode
        parameter['x_method'] = x_method
        parameter['y_method'] = y_method        
        
        return Meaningfull_prediction_support(parameter)
    
    def Generate_center_ids_for_RBFN(self):
        RBFN_centers_Ab_in_Ag = []
        for keys, value in self.Ab_in_Ag_cluster_ids.items():
            Ab = value[0][0]
            RBFN_centers_Ab_in_Ag.append(self.training_data[Ab])
        self.RBFN_centers_Ab_in_Ag = RBFN_centers_Ab_in_Ag

        RBFN_centers_Ag_in_Ab = []
        for keys, value in self.Ag_in_Ab_cluster_ids.items():
            Ag = value[0][0]
            RBFN_centers_Ag_in_Ab.append(self.training_data[Ag])
        self.RBFN_centers_Ag_in_Ab = RBFN_centers_Ag_in_Ab                
            
##################################################################################
##################################################################################


tryout =  DataAnalysis(training, testing)                
tryout.Calculate_distance_matrix()            
tryout.Get_hierarchical_clusters() 
tryout.Show_elbow()
tryout.Cut_n_clusters(lower_n_cluster=20, upper_n_cluster=100)
tryout.Get_cluster_ids()
len(tryout.Ab_cluster_ids)
len(tryout.Ag_cluster_ids)
tryout.Cluster_within_cluter(min_cluster_number=3, max_cluster_number=10)
keys1 = list(tryout.Ab_in_Ag_cluster_ids.keys())
keys2 = list(tryout.Ag_in_Ab_cluster_ids.keys())
len(keys1)
len(keys2)
len(tryout.Ab_cluster_ids)
tryout.Selected_support_predcition(top_n_clusters_support=15, 
                                   n_prediction= math.floor(0.1 * len(tryout.Ag_cluster_ids)),
                                   mode='Ag_to_Ab', x_method='average', y_method='single')
tryout.Selected_support_predcition(top_n_clusters_support=15, 
                                   n_prediction= math.floor(0.1 * len(tryout.Ab_cluster_ids)),
                                   mode='Ab_to_Ag', x_method='average', y_method='single')
tryout.Generate_center_ids_for_RBFN()
len(tryout.RBFN_centers_Ab_in_Ag)
len(tryout.RBFN_centers_Ag_in_Ab)

with open('Ab_Ag_centers_from_NewDataAnalysis', 'w') as f:
    json.dump(tryout.Ab_Ag_centers, f)
 
#################################################################################
###############################################################################
'''
Try to use the paires with lower contact number, but it does not work well.
'''
#len(negative_samples)
#sliced_testing = testing[:30]
#sliced_testing
#negative_centers = DataAnalysis(negative_samples, sliced_testing)  
#negative_centers.Calculate_distance_matrix()            
#negative_centers.Get_hierarchical_clusters() 
#negative_centers.Show_elbow()
#negative_centers.Cut_n_clusters(lower_n_cluster=10, upper_n_cluster=20)
#negative_centers.Get_cluster_ids()
#len(negative_centers.Ab_cluster_ids)
#len(negative_centers.Ag_cluster_ids)
#negative_centers.Cluster_within_cluter(min_cluster_number=3, max_cluster_number=10)  
#negative_centers.Selected_support_predcition(top_n_clusters_support=15, 
#                                   n_prediction= math.floor(0.1 * len(negative_centers.Ag_cluster_ids)),
#                                   mode='Ag_to_Ab', x_method='average', y_method='single')
#negative_centers.Generate_center_ids_for_RBFN()
#len(negative_centers.RBFN_centers_Ab_in_Ag)
#len(negative_centers.RBFN_centers_Ag_in_Ab)    
    
#################################################################################
####################################################################################    
    
    
#tryout.Ab_in_Ag_cluster_ids[keys1[40]]
#tryout.Ag_in_Ab_cluster_ids[keys2[40]]
# 
# simple check
#def simple_check(dictionary, cluster):
#    for key, value in dictionary.items():
#        list_key = list(key)
#        list_key.sort
#        ids = []
#        for i in value:
#            ids.extend(i)
#        ids.sort()
#        if ids == list_key and ids in cluster:
#            print('Good')
#        else:
#            print('Not good')

#simple_check(tryout.Ab_in_Ag_cluster_ids, tryout.Ag_cluster_ids)
#simple_check(tryout.Ag_in_Ab_cluster_ids, tryout.Ab_cluster_ids)
#tryout.Ab_in_Ag_cluster_ids[keys1[0]]
#tryout.Ag_cluster_ids[0]
#tryout.Ab_in_Ag_cluster_ids[keys1[-1]]
#tryout.Ag_cluster_ids[-1]
#tryout.Ag_in_Ab_cluster_ids[keys2[80]]
#n = 0
#for i in keys1:
#    n += len(i)
#n
#m = 0
#for i in keys2:
#    m += len(i)
#m   

#keys1[:2]
#keys2[:2]
#tryout.Ab_cluster_ids[90]
#tryout.Ag_cluster_ids[:2]
#len(tryout.Ab_cluster_ids)

#support_Ag_dict = {}
#for key, value in tryout.Ag_in_Ab_cluster_ids.items():
#    support_Ag_dict[key]  = []
#    for i in range(min(3, len(value))):
#        support_Ag_dict[key].extend(value[i])
#support_Ag = [0]*len(tryout.Ab_cluster_ids)
#for i in range(len(support_Ag)):
#    for key, value in support_Ag_dict.items():
#        if key == tuple(tryout.Ab_cluster_ids[i]):
#            support_Ag[i] = value
#
#
#keys = list(support_Ag_dict.keys())
#tryout.Ag_in_Ab_cluster_ids[keys[8]]
#support_Ag_dict[keys[8]]
#support_Ag[8]
#tryout.Ab_cluster_ids[8]













#tryout.Ag_to_Ab_prediction_correct_rate(top_x_n_cluster=1, top_y_n_cluster=2)
#tryout.Ab_to_Ag_prediction_correct_rate(top_x_n_cluster=1, top_y_n_cluster=2)






######################################################################################
#####################################################################################
'''
Create a baseline by shuffling the data.
'''

def Shuffle(training_data):#permuataion method Just mass up with the training data, see how the prediction performs
    import copy
    data = copy.deepcopy(training_data)
    Ab_data = []
    for parepi in data:
        Ab_data.append(parepi[0])
        
    random.shuffle(Ab_data)
    
    for i in range(len(data)):
        data[i][0] = Ab_data[i]
        
    return data

def Random_sample(training_data, length):
    TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
                    ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
                    ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]
    
    random_training = []
    for i in range(len(training_data)):
        rsample = []
        for j in range(length):
            rsample.append(random.sample(TripleSingle, 1)[0][0])
        random_training.append([rsample, training_data[i][1], training_data[i][2], training_data[i][3]])       
    return random_training 

r_training = Random_sample(training, length = 2)        
r_training[:6]
tryout.training_data[:6]

shuffled = Shuffle(training)
training[:6]
shuffled[:6]


Ag_to_Ab = []
Ab_to_Ag = []
for i in range(6):
    r_training = Shuffle(training) 
    
    baseline = DataAnalysis(r_training, testing)                
    baseline.Calculate_distance_matrix()            
    baseline.Get_hierarchical_clusters() 
    baseline.Show_elbow()
    baseline.Cut_n_clusters(lower_n_cluster=20, upper_n_cluster=100)
    baseline.Get_cluster_ids()
    baseline.Cluster_within_cluter(min_cluster_number=3, max_cluster_number=10)
#    baseline = list(tryout.Ab_in_Ag_cluster_ids.keys())
#    baseline = list(tryout.Ag_in_Ab_cluster_ids.keys())
    #baseline.Ag_to_Ab_prediction_correct_rate(top_x_n_cluster=1, top_y_n_cluster=2)
    
    Ag_to_Ab.append(baseline.Selected_support_predcition(top_n_clusters_support=15,
                                   n_prediction= math.floor(0.1 * len(tryout.Ag_cluster_ids)),
                                   mode='Ag_to_Ab', x_method='average', y_method='single'))
    Ab_to_Ag.append(baseline.Selected_support_predcition(top_n_clusters_support=15, 
                                  n_prediction= math.floor(0.1 * len(tryout.Ab_cluster_ids)),
                                   mode='Ab_to_Ag', x_method='average', y_method='single'))

np.average(Ag_to_Ab)
np.average(Ab_to_Ag)
########################################################################################
#################################################################################
'''
Use the above baseline to generate negative centers is a relatively good choice.
'''
baseline.Generate_center_ids_for_RBFN()
len(baseline.RBFN_centers_Ab_in_Ag)
len(baseline.RBFN_centers_Ag_in_Ab)
baseline.RBFN_centers_Ab_in_Ag[:6]
negative_centers_Ab_in_Ag = []
for parepi in baseline.RBFN_centers_Ab_in_Ag:
    negative_centers_Ab_in_Ag.append([parepi[0], parepi[1], 0, parepi[3]])
negative_centers_Ab_in_Ag[:6]
len(negative_centers_Ab_in_Ag)
####################################################################################
###################################################################################
'''
Generate the center for RBFN and save it.
'''
RBFN_centers_Ab_in_Ag = []
for i in tryout.RBFN_centers_Ab_in_Ag:
    RBFN_centers_Ab_in_Ag.append(i)
for j in negative_centers_Ab_in_Ag:
    RBFN_centers_Ab_in_Ag.append(j)
len(RBFN_centers_Ab_in_Ag)
RBFN_centers_Ab_in_Ag

with open('RBFN_centers_Ab_in_Ag', 'w') as f:
    json.dump(RBFN_centers_Ab_in_Ag, f)

