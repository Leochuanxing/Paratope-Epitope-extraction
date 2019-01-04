
##########################################################
'''
Import modules
'''
import random
import numpy as np
import os
import json
import math
import copy
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
        for j in range(i, len(seq_list)):
            if seq_list[i] != seq_list[j]:
                distance_matrix[i, j] =1 - (4*l + aligner.score(seq_list[i], seq_list[j]))/(15*l)#Global
#                 distance_matrix[i, j] =1 -  aligner.score(seq_list[i], seq_list[j])/(15*l)#Local
#                distance_matrix[i, j] = - aligner.score(seq_list[i], seq_list[j])#Unnormalized
                distance_matrix[j, i] = distance_matrix[i, j]
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
    unordered_Z = linkage(condenced_distance_matrix, method = 'single')
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

    

    
    def Selected_support_predcition(self, top_n_clusters_support=3, n_prediction = 3, mode = 'Ag_to_Ab',
                                    x_method = 'complete', y_method = 'average'):


                    
        parameter = {}
        parameter['len_training_data'] = len(self.training_data)
        parameter['test_testing_data'] = len(self.testing_data)
        parameter['Ag_cluster_ids'] = self.Ag_cluster_ids
        parameter['Ab_cluster_ids'] = self.Ab_cluster_ids
        parameter['n_prediction'] = n_prediction
        parameter['all_Ab_distance_martix'] = self.all_Ab_distance_martix
        parameter['all_Ag_distance_martix'] = self.all_Ag_distance_martix
        parameter['support_Ab'] = self.Ag_cluster_ids
        parameter['support_Ag'] = self.Ab_cluster_ids
        parameter['mode'] = mode
        parameter['x_method'] = x_method
        parameter['y_method'] = y_method        
        
        return Meaningfull_prediction_support(parameter)
    
               
            
##################################################################################
##################################################################################





 
#################################################################################
###############################################################################
 
  





######################################################################################
#####################################################################################
'''
Create a baseline by shuffling the data.
'''
#
#def Shuffle(training_data):#permuataion method Just mass up with the training data, see how the prediction performs
#    import copy
#    data = copy.deepcopy(training_data)
#    Ab_data = []
#    for parepi in data:
#        Ab_data.append(parepi[0])
#        
#    random.shuffle(Ab_data)
#    
#    for i in range(len(data)):
#        data[i][0] = Ab_data[i]
#        
#    return data
#
#def Random_sample(training_data, length):
#    TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
#                    ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
#                    ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]   
#    random_training = []
#    for i in range(len(training_data)):
#        rsample = []
#        for j in range(length):
#            rsample.append(random.sample(TripleSingle, 1)[0][0])
#        random_training.append([rsample, training_data[i][1], training_data[i][2], training_data[i][3]])       
#    return random_training 
#
#r_training = Random_sample(training, length = 2)        
#
#
#
#
#Ag_to_Ab = []
#Ab_to_Ag = []
#for i in range(2):
#    r_training = Random_sample(training, 3) 
#    
#    baseline = DataAnalysis(r_training, testing)                
#    baseline.Calculate_distance_matrix()            
#    baseline.Get_hierarchical_clusters() 
#    baseline.Show_elbow()
#    baseline.Cut_n_clusters(lower_n_cluster=20, upper_n_cluster=5000)
#    baseline.Get_cluster_ids()
#    
#    Ag_to_Ab.append(baseline.Selected_support_predcition(top_n_clusters_support=15,
#                                   n_prediction= math.floor(0.1 * len(tryout.Ag_cluster_ids)),
#                                   mode='Ag_to_Ab', x_method='single', y_method='single'))
#    Ab_to_Ag.append(baseline.Selected_support_predcition(top_n_clusters_support=15, 
#                                  n_prediction= math.floor(0.1 * len(tryout.Ab_cluster_ids)),
#                                   mode='Ab_to_Ag', x_method='single', y_method='single'))
#
#np.average(Ag_to_Ab)
#np.average(Ab_to_Ag)
#Ab_to_Ag
#results = []
#with open('Results', 'w') as f:
#    json.dump(results, f)
#def Results_collect(lower_contact_limit,x_method, y_method, top_n_clusters_support, 
#                    mode1, correct_rate1, baseline1, mode2,correct_rate2, baseline2):
#    with open('Results', 'r') as f:
#        results = json.load(f)
#    results.append([lower_contact_limit,  x_method, y_method, top_n_clusters_support, 
#                    mode1, correct_rate1, baseline1, mode2, correct_rate2, baseline2])
#    with open('Results', 'w') as f:
#        json.dump(results, f)
#    return print('Results saved')
#
#Results_collect(lower_contact_limit=9, x_method='single', y_method='single', top_n_clusters_support=15,
#                mode1='Ag_to_Ab', correct_rate1=0.72, baseline1 = 0.11, mode2='Ab_to_Ag', correct_rate2=0.78, baseline2=0.77)
#with open('Results', 'r') as f:
#    results = json.load(f)
#results
########################################################################################
#################################################################################
'''
Define a main function, which gives what I want in one run
'''
def main():
    '''
    This file is to analyse the data generated from GetDataForAnalysis
    '''
    ##########################################################
    
    os.chdir("/home/leo/Documents/Database/Pipeline/Ready_2_2_1_1")
    
    
    with open('training_2_2_1_1_all_jump', 'r') as f:
        training = json.load(f)
    with open('testing_2_2_1_1_all_jump', 'r') as f:
        testing = json.load(f)
    with open('negative_samples_2_2_1_1', 'r') as f:
        negative_samples = json.load(f)

####################################################################################
    tryout =  DataAnalysis(training, testing)                
    tryout.Calculate_distance_matrix()            
    tryout.Get_hierarchical_clusters() 
    tryout.Show_elbow()
    tryout.Cut_n_clusters(lower_n_cluster=20, upper_n_cluster=1000000)
    tryout.Get_cluster_ids()
    len(tryout.Ab_cluster_ids)
    len(tryout.Ag_cluster_ids)
    
    Ag_to_Ab_rate = tryout.Selected_support_predcition(top_n_clusters_support=15, 
                                       n_prediction= math.floor(0.1 * len(tryout.Ag_cluster_ids)),
                                       mode='Ag_to_Ab', x_method='single', y_method='single')
    Ab_to_Ag_rate = tryout.Selected_support_predcition(top_n_clusters_support=15,
                                       n_prediction= math.floor(0.1 * len(tryout.Ab_cluster_ids)),
                                       mode='Ab_to_Ag', x_method='single', y_method='single')
    
    baseline = DataAnalysis(negative_samples, testing)                
    baseline.Calculate_distance_matrix()            
    baseline.Get_hierarchical_clusters() 
    baseline.Show_elbow()
    baseline.Cut_n_clusters(lower_n_cluster=20, upper_n_cluster=5000)
    baseline.Get_cluster_ids()
    
    baseline_Ag_to_Ab_rate = baseline.Selected_support_predcition(top_n_clusters_support=15,
                                   n_prediction= math.floor(0.1 * len(baseline.Ag_cluster_ids)),
                                   mode='Ag_to_Ab', x_method='single', y_method='single')
    baseline_Ab_to_Ag_rate = baseline.Selected_support_predcition(top_n_clusters_support=15, 
                                  n_prediction= math.floor(0.1 * len(baseline.Ab_cluster_ids)),
                                   mode='Ab_to_Ag', x_method='single', y_method='single')

    return ['Hierarchical clustering', 'Ag_to_Ab', len(tryout.Ag_cluster_ids), math.floor(0.1 * len(tryout.Ag_cluster_ids)),\
            Ag_to_Ab_rate, 'Ab_to_Ag', len(tryout.Ab_cluster_ids), math.floor(0.1 * len(tryout.Ab_cluster_ids)), Ab_to_Ag_rate],\
            ['Hierarchical clustering baseline', 'Ag_to_Ab', len(baseline.Ag_cluster_ids), math.floor(0.1 * len(baseline.Ag_cluster_ids)),\
            baseline_Ag_to_Ab_rate, 'Ab_to_Ag', len(baseline.Ab_cluster_ids), math.floor(0.1 * len(baseline.Ab_cluster_ids)), baseline_Ab_to_Ag_rate]
 ################################################################################################

################################################################################################           
if __name__ == '__main__':
    Hcluster_results, baseline_results = main()
    results={}
    results['Hcluster'] = Hcluster_results
    results['baseline'] = baseline_results
    
    os.chdir("/home/leo/Documents/Database/Pipeline/Ready_2_2_1_1")
    with open ('Hcluster_results', 'w') as f:
        json.dump(results, f)

#os.chdir("/home/leo/Documents/Database/Pipeline/Ready_2_2_1_1")
#with open ('Hcluster_results', 'r') as f:
#    results = json.load(f)
#results
