# Do some preparations and load the data
import os
os.chdir("C:\\Users\\leo\\Documents\\Research\\Database\\Pipeline\\Analysis2-1-1-1")
os.listdir()
import json
with open('Homo_2_1_1_1', 'r') as f:
    Homo_2_1_1_1 = json.load(f)
with open('Mouse_2_1_1_1', 'r') as f:
    Mouse_2_1_1_1 = json.load(f)
########################################################
'''
In order to make this analysis more easily to be generalized to other frame length and 
other cases like 0-free or 1-free, we should make this analysis as generalized as possible
'''
# First we would like to cluster reference frames
# second, according to the clustered reference frames, allocate the complementary aas 
# into different clusters.
#Third, cluster the the complementary aas in each cluster into different clusters 
# Fourth, according to the size of the clusters fo the complementary aas in each 
# class to make predictions
# FiFth, can we show that the complemantary aas allocated to each cluster of the reference 
# frame can form a reasonable clustering just by referring to the sequence of the complementary 
# aas. To what extent is this consistant with the allocation, we can do a smimple calculation 
# of the between cluster distance and inter cluster distance.
####################################################################
#lets do the first step, we have to use a lot of function in module AlignmentConstraint
# change the sequence into Seq object
# define a simple Triple_single_converter
'''
inputs: amino_acids, a list either in the form of ['ASP', 'ASP'], or in the form of ['A', 'L']
return: converted, a list, triple-letter or single letter notation corresponding to the amino_acids
'''
def T_S_converter(amino_acids):
    TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
   ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
   ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]
    # Find the form of the amino_acids
    if len(amino_acids[0]) == 3:
        ts = 0 #ts means triple or single
    else:
        ts = 1
    # convert
    converted = []
    for aa in amino_acids:
        for table in TripleSingle:
            if aa == table[ts]:
                converted.append(table[1-ts])
    return converted

'''
inputs: sequence, a list in the form of ['ALA', 'ARG', 'TYR']
return: Seq_ob, a Seq object, in the form of Seq(ARY, IUPACprotein())
'''
# Import needed packages
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
def To_Seq(sequence):
    TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
       ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
       ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]
    single_letter_seq = ''
    for i in sequence:
        for j in TripleSingle:
            if j[0] == i:
                single_letter_seq += j[1] 
    Seq_ob = Seq(single_letter_seq, IUPAC.protein)
    return Seq_ob

# Set the aligner
from Bio import Align
from Bio.SubsMat.MatrixInfo import blosum62
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = blosum62
aligner.open_gap_score = -5
aligner.extend_gap_score = -1
# creat similarity matrix 
'''
inputs: 
'''
import numpy as np
def Similarity_matrix(Seq_ob):
    similarity_matrix = np.zeros([len(Seq_ob), len(Seq_ob)])
    for i in range(len(Seq_ob)):
        for j in range(len(Seq_ob)):
            similarity_matrix[i,j] = aligner.score(Seq_ob[i], Seq_ob[j])
    return similarity_matrix

# Calculate the cluster dist
# here we just copy and paste alot from AlignmentConstraint
# define a function to caltulate the distances between clusters
'''
inputs: clusters, a list in the form of [[0, 1, 2], [3,5], ...]
        dist_matrix, a distance matrix, returns from Dist_matrix
return: cluster_dist, a list in the form of [[[i,j],dist], ....]
        i, j, means the i th and the j th cluster, dist is the distance between 
        those two clusters, here we use the least similarities, and the results 
        given are the paires with the highest similarity
'''
def Cluster_dist(clusters, similarity_matrix):
    # Creat an empty container to contain the results
    cluster_dist = []
    for i in range(len(clusters)-1):
        # Creat a variable to contain the possible maximum score and the corresponding
        # position temporarily
        score = -10000
        pos = -1
        for j in range(i+1, len(clusters)):
            # calculate the distance between cluster i and cluster j
            # here we use the smallest similarity score
            # creat a variable to contain this smalles score
            similarity = 10000           
            for k in clusters[i]:
                for l in clusters[j]:
                    if similarity_matrix[k,l] < similarity:
                        similarity = similarity_matrix[k,l]
            # find the clusters with the highest similarity                        
            if similarity > score:
                pos = j
                score = similarity
        cluster_dist.append([[i, pos], score])
    return cluster_dist
# the original Merge function has to be modified
'''
Inputs: clusters, a list in the form of [[0, 1, 2], [3,5], ...], corresponding 
                to cluster_dist. 0, 1, 2,...are the positions of the items need to be 
                clustered, they can be applied directly on the similarity_matrix
        cluster_number, an integer, the number of clusters returned should be less or equal 
                 to this number, and the nearest to this number
        highest_score, a number, gives the nearest upper limit of the similarity score
               cluster_number and highest_score can be used to control 
               the clustering process simutaneously.

return, merged_clusters, a list of the form [[0, 1, 2], [3,5], ...], where 0, 1, 2
        ... are corresponding to the items need to be clustered
        highest, a float, gives the largest similarity that to be merged in this 
        merge step.
'''

def Merge(clusters, similarity_matrix,  cluster_number = None, highest_score = None):
    # find the paires with the highest similarity
    cluster_dist = Cluster_dist(clusters, similarity_matrix)
    # Generate the list need to be merged with the maximum similarity, for example
    #[[[0, 45], 37.0], [[1, 2], 50.0], [[2, 389], 34.0], [[3, 4], 50.0]]
    # we select outh the emlemts with the highest score and merge
    # select out the elements with the highest score
    highest = -10000
    for dist in cluster_dist:
        if dist[1] > highest:
            highest = dist[1]
    # make two conditions
    condition1 = cluster_number != None and len(clusters) > cluster_number 
    condition2 = highest_score != None and highest >= highest_score 
    if condition1 or condition2:                   
        bag = []
        need_merge = []
        for i in cluster_dist:    
            if i[1] == highest:
                if i[0][0] not in bag and i[0][1] not in bag:
                    need_merge.append(i[0])
                bag.extend(i[0])
                
        # merge the indices in the need_merge
        for merge in need_merge:
            clusters[merge[0]].extend(clusters[merge[1]])
            clusters[merge[1]] = 'D'
        while 'D' in clusters:
            clusters.remove('D')                                 
        clusters.sort()            
        if len(clusters) == 1:
            return clusters, highest
        else:
            return Merge(clusters, similarity_matrix, cluster_number, highest_score ) 
    else:
        return clusters, highest



'''
inputs: sequence_list, a list [[ALA, ARG], [ALA, GLU], ], or in single letter form
        cluster_number, the nearest upper limit of clusters 
        highest_score, the nearst uper limit of similarity score
        #Only one of highest_score or cluster_number can be used to control the merging process
returns: clusters, a list, in the form of ['GLWY', 'AALV',...]
         highest, the nearest upper limit of the inter-cluster similarily. If the inter-cluster
         similarity is no less than this value, two clusters will be merged
'''
def MainMerge(sequence_list, cluster_number = None, highest_score = None):
    # convert to Seq
    if len(sequence_list[0][0]) == 3:
        Seq_list = []
        for seq in sequence_list:
            temp = ''
            for i in T_S_converter(seq):
                temp += i
            Seq_list.append(Seq(temp, IUPAC.protein))
    else:
        Seq_list = []
        for seq in sequence_list:
            temp = ''
            for i in seq:
                temp += i
            Seq_list.append(Seq(temp, IUPAC.protein)) 
           
    similarity_matrix = Similarity_matrix(Seq_list)
    starters = []
    for i in range(len(Seq_list)):
        starters.append([i])
    clusters_ind, highest = Merge(starters, similarity_matrix, cluster_number, highest_score)
    
    # Convert the returned clusters in term of sigle-letter sequences
    clusters = []
    for clus in clusters_ind:
        temp = []
        for ind in clus:
            temp.append(str(Seq_list[ind]))
        clusters.append(temp)            
    
    return clusters, highest
        
#Need to check this function
# Now it looks good after hours of debugging   

#######################################################################
# Now, a simple implementation
# Merge the data
len(Homo_2_1_1_1) 
len(Mouse_2_1_1_1)
import copy; data = copy.deepcopy(Homo_2_1_1_1)
data.extend(Mouse_2_1_1_1)
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Hold about 10% out for cross validation
import random
hold_out_data = random.sample(data, 150)
len(hold_out_data)
# delete the hold_out_data from the data
for i in hold_out_data:
    data.remove(i)
len(data)
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# Here the length of the complementary aas are different, in order to make the 
# alignment make sense more, we add the same amimo acid to the single amino acid
# to make the length 2
for i in range(len(data)):
    if len(data[i][0]) == 1:
        data[i][0].extend(data[i][0])
testing_data = copy.deepcopy(data[:30])
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\        
# Two ways to cluster, one is to set the upper limit of the cluster similarity
# the other one is to set the upper limit of the number of clusters, in order to 
# find a proper upper limit of the cluster distance, we have to do a simulation
# the simulation is that, we generate all possible biaas, and calcualte the largest 
# similarity when we clustered them into the number of clusters we wanted. 
# we use this number as the upper limit of the cluster similarity.
# Lets do a simulation first
'''
inputs: num_clusters, an integer, nearest upper bound of the number of clusters
        highest, a number(integer or float), gives the nearest upper bound of similarity score
        # Only one of num_clusters and highest can be used to control the simulation.
        aa_length, give the length of the all possible aas need to be simulated, here the
        aa_length should not be very long, because the computation cost is exponential
returns: cluster, the clusters under the control of num_clusters or highest
          score, the simulated score
'''

def Simulation_cluster(aa_length, num_clusters = None, highest = None):
    # The 20 amino acids
    twenty = ['S', 'N', 'D', 'T', 'V', 'R', 'F', 'W', 'Q', 'P',
     'K', 'H', 'E', 'Y', 'M', 'C', 'G', 'I', 'A', 'L'] 
    # generate all possible aas of gien length
    starter = copy.deepcopy(twenty)      
    while len(starter[0]) <aa_length:
        temp = []
        for i in twenty:
            for j in range(len(starter)):
                temp.append(starter[j]+i)
        starter = copy.deepcopy(temp)
    # Do a sample in case the starter is too big
    if len(starter) >1000:
        import random; sample = random.sample(starter, 1000)
    else:
        sample = starter
        
    cluster, score =  MainMerge(sample, num_clusters, highest ) 
    
    return cluster, score
# Check the function
cluster, score = Simulation_cluster(aa_length = 2, num_clusters = None, highest = 6.0)
len(cluster)    
score
cluster[:10]
# check if the function works this time
#n = 0
#for i in cluster[1]:
#    for j in cluster[1]:
#        if i != j and aligner.score(i,j) > n:
#            n = aligner.score(i,j) 
# it looks ok!
# Convert the returned clusters in term of sigle-letter sequences


# it looks good
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

####################################################################

class Complementary_cluster(object):
    '''
    data, a list, gives the data need to be analyzed, in the form of [[['ASP', 'ASP'], ['SER'], 7, 17],..]
        it is the return of the module GetDataForAnalysis
    ref_chain, a string, takes the value of either 'Ab' or 'Ag', gives information of the 
        name of the reference chain
    top_cluster, give  the number of the top clusters that should be 
          used to calculate the percentage, those cluster take up the total number of complementary
          sequences
    feed_ref, is a list of ref sequence, used as an independent variable to make predictions 
          about the complementary amino acids, in the form of [[ALA, TYR], [SER, PRO],..]
    top_prediction, an integer, give the number of top predictions, our prediction is given in
        term of clusters, it tells you which cluster it the complementary sequence of the feed_ref
        belong to.
    feed_complement, a list of complementary aas, they are used as the independent variable
        to make predictions about the ref sequences, it shares the same top_prediction with 
        feed_ref. And the prediction is given in terms of clusters of the feed_ref   
    cluster_para, a dictionary in the form of {'1': [20, 7], ...} where '1' is the key
        corresponding to the keys of Grouped_dict or group_no_clu. [20, None], gives the 
        parameter of clustering control parameters, cluster_number and highest_score.
    
    '''
    def __init__(self, data, ref_chain = 'Ag', top_cluster = 0.25, feed_ref = None, 
                 top_prediction = 5, feed_complement = None, data_single_letter = None,
                 clustered_ref = None, ref_cluster_score =None, key_dict = None, 
                 Grouped_dict = None, group_no_clu = None, frequency = None, cluster_para = None,
                 clustered_complementary = None, complementary_cluster_score = None, 
                 ):
        self.data = data
        if ref_chain == 'Ag':
            self.ref_ind = 1
        else:
            self.ref_ind = 0
        self.top_cluster = top_cluster
        self.feed_ref = feed_ref
        self.top_prediction = top_prediction
        self.feed_complement = feed_complement
        self.cluster_para = cluster_para  #
    # before we do anything, lets convert the notation of the amino acids into single letter form
    # the return is exactly the same as the data, except that it is in sigle-letter notation
    def Convert_data(self):
        data_single_letter = []
        for core in self.data:
            data_single_letter.append([T_S_converter(core[0]), T_S_converter(core[1]), core[2], core[3]])
        self.data_single_letter = data_single_letter
    
        # cluster the reference frame
    '''
    Inputs: the constraints, highest_scores = None, num_cluster = None
    Return: the clustered keys in the form of [[A, T, G],...]        
    '''
    def Cluster_ref(self, highest_score = None, cluster_number = None):
        ref_list = [ref[self.ref_ind] for ref in self.data]
        # divide up the combined 
        clustered_ref, ref_cluster_score = MainMerge(ref_list, cluster_number, highest_score)
        self.clustered_ref = clustered_ref
        self.ref_cluster_score = ref_cluster_score 
        

    # separate them into groups and store them in a dictionary, with keys the reference 
    # sequence.
    '''
    inputs: data_single_letter, a list in the form of [[['A', 'T'], ['W'], 7, 17]...]
            ref_ind, an integer, specify the position of the keys
    return: groups, a dictionary, with keys specificied by the ref_ind, and the data
            will be grouped according to the specified keys
    '''
    #[['S'], ['D', 'D', 'D', 'D'],...] the clusterd ref looks like this
    def Groupup_to_ref_cluster(self):
        # creat a key dictionary for the clusted ref
        key_dict = {}
        for i in range(len(self.clustered_ref)):
            key_dict[str(i)] = self.clustered_ref[i]
        self.key_dict = key_dict
        # load data to to dictionary with clustered ref keys
        import copy
        single_letter_data = copy.deepcopy(self.data_single_letter)
        Grouped_dict ={}
        for i in self.key_dict:
            Grouped_dict[i] = []
            for aas in self.key_dict[i]:
                for core in single_letter_data:
                    temp = ''
                    for j in core[self.ref_ind]:
                        temp += j
                    if temp == aas:
                        Grouped_dict[i].append(core[1-self.ref_ind])
                        single_letter_data.remove(core)
                        break
        self.Grouped_dict = Grouped_dict
    # Group directly without clustering the ref, this is particularly for ref with 
    # ref length 1
    def Groupup_without_ref_cluster(self):
        group_no_clu = {}
        for core in self.data_single_letter:
            if core[self.ref_ind][0] not in group_no_clu:
                group_no_clu[core[self.ref_ind][0]] = core[1-self.ref_ind]
            else:
                group_no_clu[core[self.ref_ind][0]].append(core[1-self.ref_ind])
        self.group_no_clu = group_no_clu
    #Calculate the frequency 
    def Frequency(self, dictionary):
        frequency = []
        for i in dictionary:
            frequency.append([i, len(dictionary[i])])
        self.frequency = frequency
    # Cluster the complementary aas
    def Cluster_complementary(self, dictionary):
        clustered_complementary = {}
        for i in dictionary:
            clustered_complementary[i], complementary_cluster_score = MainMerge(dictionary[i],self.cluster_para[i][0], 
                                                self.cluster_para[i][1] )
        self.clustered_complementary = clustered_complementary
        self.complementary_cluster_score = complementary_cluster_score
            
            
        
            
            

                        

analysis = Complementary_cluster(data, ref_chain='Ag')
analysis.Convert_data()
analysis.data_single_letter
analysis.ref_ind 
#analysis.Cluster_ref(highest_score=None, cluster_number=13)
analysis.Groupup_without_ref_cluster()   
len(analysis.group_no_clu)
#analysis.clustered_ref[:3]
#analysis.Groupup_to_ref_cluster()
#analysis.Grouped_dict.keys()
#analysis.key_dict
analysis.Frequency(analysis.group_no_clu)
analysis.frequency.sort(key = lambda x:x[1], reverse = True)
analysis.frequency
cluster, score = MainMerge(analysis.group_no_clu['Q'], 12, None)
length = []
for i in cluster:
    length.append(len(i))
length.sort(reverse = True)
length
