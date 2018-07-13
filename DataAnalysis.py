# Do some preparations and load the data
import copy
import os
os.chdir("C:\\Users\\leo\\Documents\\Research\\Database\\Pipeline\\Analysis2_1")
os.listdir()
import json

#with open('data_2_1', 'r') as f:
#    Mouse_2_1= json.load(f)
#with open('data_2_1', 'r') as f:
#    Homo_2_1= json.load(f)
#with open('Mouse_2_2_1_2', 'r') as f:
#    Mouse_2_2_1_2 = json.load(f)

with open('training', 'r') as f:
    training = json.load(f)
with open('testing', 'r') as f:
    testing = json.load(f)
#    
len(training)
len(testing)
testing[:10]
training[:10]
#######################################################
# Divide the data up into training data and testing data
# Here the length of the complementary aas are different, in order to make the 
# alignment make sense more, we add the same amimo acid to the single amino acid
# to make the length 2
              
#len(Homo_2_1) 
#len(Mouse_2_1)
#data = copy.deepcopy(Mouse_2_1)
#data.extend(Homo_2_1)
#len(data)
#data[:10]

#for i in range(len(data)):
#    if len(data[i][0]) == 1:
#        data[i][0].extend(data[i][0])
        
#data2_2 = []
#for i in range(len(data)):
#    if len(data[i][0]) == 2:
#        data2_2.append(data[i])
#len(data2_2)
##        
## Hold about 10% out for cross validation
#import random
#hold_out_data = random.sample(data, 130)
#len(hold_out_data)
#hold_out_data[:10]

## delete the hold_out_data from the data
#for i in hold_out_data:
#    data.remove(i)
#len(data)
#data[:10]
##
##
# save all the above data
#with open('training2_1', 'w') as f:
#    json.dump(data, f)
#with open('testing2_1', 'w') as f:
#    json.dump(hold_out_data, f)
    
#with open('function_checking_data', 'w') as f:
#    json.dump(testing_data, f)

##\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

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
inputs: sequence_list, a list [['ALA', 'ARG'], ['ALA', 'GLU'], ], or in single letter form
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
          

#######################################################################

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
#cluster, score = Simulation_cluster(aa_length = 2, num_clusters = None, highest = 6.0)
#len(cluster)    
#score
#cluster[:10]
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
    def __init__(self, data, testing_data, ref_chain = 'Ag', feed_complement = None, data_single_letter = None,
                 clustered_ref = None, ref_cluster_score =None, key_dict = None, Grouped_dict = None, 
                  frequency = None, cluster_para = None,clustered_complementary = None, 
                 complementary_cluster_score = None, processed_testing_data = None, reverse_correct_rate =None,
                 intra_inter_cluster_similarity = None):
        self.data = data
        self.testing_data = testing_data
        if ref_chain == 'Ag':
            self.ref_ind = 1
        else:
            self.ref_ind = 0
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
    def Cluster_ref(self, cluster_number = None,  highest_score = None):
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
    #[['S'], ['D', 'D', 'D', 'D'],...] the clusterd ref looks like this.
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
    # ref length 1, This function need to be deleted
#    def Groupup_without_ref_cluster(self):
#        group_no_clu = {}
#        for core in self.data_single_letter:
#            if core[self.ref_ind][0] not in group_no_clu:
#                group_no_clu[core[self.ref_ind][0]] = core[1-self.ref_ind]
#            else:
#                group_no_clu[core[self.ref_ind][0]].append(core[1-self.ref_ind])
#        self.group_no_clu = group_no_clu
    #Calculate the frequency 
    def Frequency(self, dictionary):
        frequency = []
        for i in dictionary:
            frequency.append([i, len(dictionary[i])])
        self.frequency = frequency
    # Cluster the complementary aas
    def Cluster_complementary(self):
        clustered_complementary = {}
        complementary_cluster_score = {}
        for i in self.Grouped_dict:
            clustered_complementary[i], complementary_cluster_score[i] = MainMerge(self.Grouped_dict[i],self.cluster_para[i][0], 
                                                self.cluster_para[i][1] )
            # sort according to the length before return
            clustered_complementary[i].sort(key = lambda x:len(x), reverse = True)
        self.clustered_complementary = clustered_complementary
        self.complementary_cluster_score = complementary_cluster_score
    #Find the simulation cluster number and the merging scores
    # The simulation is randomly drawn from the combination of coresponding aa_length
    # To further improving the results more persuasive, the data should be drwan from
    # the protein where those cores come from 
    '''leave this out righ now''' 
    # Define a function to find the percentage the num_cluster takes up of 
    # all the core contacts
    def Top_cluster_percentage(self, num_cluster):#bug here, watch out!####
        # creat an empty container to contain all the top 'num_cluster' percentage
        top_n_cluster_percentage = []
        for i in self.clustered_complementary:
            self.clustered_complementary[i].sort(key = lambda x:len(x), reverse = True)
        for i in self.frequency:
            t = 0
            for j in range(num_cluster):
                t += len(self.clustered_complementary[i[0]][j])
            top_n_cluster_percentage.append([i[0], round(t/i[1], 2)])
            top_n_cluster_percentage.sort(key = lambda x:x[1], reverse = True)
        self.top_n_cluster_percentage = top_n_cluster_percentage
                    
    # Processing the testing data, majorly to transform the data into single-letter form
    # and conbine the reference into a string
    def Processing_testing_data(self):
        processed_testing_data = []
        for core in self.testing_data:
            processed_testing_data.append([T_S_converter(core[0]), T_S_converter(core[1]), core[2], core[3]])
        self.processed_testing_data = processed_testing_data
    # define a function to make predictions about the paratope on basis of the given 
    # core of the epitope, we use the tope 'top_prediction'  number of predictions of 
    # the cluster to make prediction

    # define a function to to make prediction        
    def Prediction(self): 
        # load data from self object
        testing_data = self.processed_testing_data
        key_dict = self.key_dict
        trained = self.clustered_complementary
        # sort the trained according to the size of the clusters
        for i in trained:
            trained[i].sort(key = lambda x:len(x), reverse = True)
        # the results of the prediction will be stored in 
        all_predictions = []
        #Extracting information to match up with the keys, and store this information in temp
        for testing_core in testing_data:
            temp = ''
            for i in testing_core[self.ref_ind]:
                temp += i
            #match up the keys by calculating the distances between the clustered ref and testing
            # reference aas, here we use the average similarity, because the keys store in key_dict has information
            # about frequencies and this frequency is consistent with the number of complementary seqs
            # in the dictionary.
            key_match = ''
            key_match_score = -1000
            for key in key_dict:
                # calculate the average similarity between key_dict[key] and temp
                t = 0
                for seq in key_dict[key]:
                    t += aligner.score(Seq(temp, IUPAC.protein), Seq(seq, IUPAC.protein))
                mean = t/len(key_dict[key])
                # update the key_match_score and the key_match
                if mean > key_match_score:
                    key_match_score = mean
                    key_match = key # this value should be returned as a prediction
            # now we have found the matched key for testing_core and store it in key_match
            # It is time to make prediction according to the clusters in trained[key_match]
            # We do it by calculating the average distance between the complementary seq and 
            # the seq in the trained cluster, to see whether the cluster with the 
            # largest similarity score falls in the top "num_prediction" clusters, if 
            # it is true, our prediction is correct. Otherwise, it is not.
            # Let us calculate the average distance
            prediction = -1
            similarity = -1000
            # prepare the testing seq will be used in calculating the average similarity
            temp = ''
            for i in testing_core[1-self.ref_ind]:
                temp += i
            for cluster_id  in range(len(trained[key_match])):#[['AT', 'ST'],['VV']...]
                #Calculate the average similarity between clusters[cluster_id] and temp
                t = 0
                for trained_seq in trained[key_match][cluster_id]:
                    t += aligner.score(Seq(temp, IUPAC.protein), Seq(trained_seq, IUPAC.protein))
                mean = t/len(trained[key_match][cluster_id])
                # update prediction and similarity
                if mean > similarity:
                    similarity = round(mean, 3)
                    prediction = cluster_id + 1
            # our prediction together together with the key_match will be stored
            # with the testing data
            all_predictions.append([testing_core[0], testing_core[1], key_match, prediction, similarity])
        self.all_predictions = all_predictions # End of this function
        
        # define a function to claculate the successful rate of prediction
    def Correct_rate(self, num_prediction):
        # any prediction value in function Prediction larger than num_prediction is wrong
        n = 0
        for prediction in self.all_predictions:
            if prediction[3] <= num_prediction:
                n +=1
        correct_rate = round(n/len(self.all_predictions), 2)
        self.correct_rate = correct_rate #End of this function

     # define a function to do reverse prediction   
    def Reverse_prediction(self, top_prediction, top_n_cluster):
        testing = self.processed_testing_data
        key_dict = self.key_dict
        # creat a new set of trained data for reverse prediction and store them in trained 
        trained = {}
        for i in self.clustered_complementary:
            # Sort the self.clustered_complementary[i], because we want the top top_n_cluster
            self.clustered_complementary[i].sort(key = lambda x:len(x), reverse = True)
            combine_cluster = []
            if len(self.clustered_complementary[i]) >= top_n_cluster:
                for cluster in self.clustered_complementary[i][:top_n_cluster]:
                    combine_cluster.extend(cluster)
            else:
                for cluster in self.clustered_complementary[i]:
                    combine_cluster.extend(cluster)                
            trained[i] = combine_cluster
##        The predicted results will be stored in predictions
#        predictions = []
        num_correct = 0
        #calculate the distance 
        for core in testing:
            #extract the complementary seq information
            temp = ''
            for i in core[1-self.ref_ind]:# pay attention to the index
                temp += i
          
            # use the smallest similarity as distance
            scores = []
            s = 1000
            for ref in trained:
                for seq in trained[ref]:
                    t = aligner.score(Seq(temp, IUPAC.protein), Seq(seq, IUPAC.protein))
                    if t < s:
                        s = t
                scores.append([ref, s]) 
            
            # sort the scores and take the top top_prediction
            scores.sort(key = lambda x:x[1], reverse = True)        
            # extract the correct label
            testing_ref_seq = ''
            for i in core[self.ref_ind]:
                testing_ref_seq += i
#            # tell whether the prediction is correct or not
#            n = 0
#            for pred in scores[:top_prediction]:# this is a bug, don't have to be in, close enough is enough
#                if temp in key_dict[pred[0]]:
#                    n = 1
#                    num_correct += 1
#                    break
            # Prediction shouldn't use in or not, should use whether it is the nearest
            # calculate the nearest similarity between temp and all the clustered keys
            # Here we use the least similarity as the similarity  score between the 
            # temp and the clustered keys
            ref_similarity = []
            for key_ref in key_dict:
#                s = 1000
#                for training_ref_seq in key_dict[key_ref]:
#                    t = aligner.score(Seq(testing_ref_seq, IUPAC.protein), Seq(training_ref_seq, IUPAC.protein))
#                    if t < s:
#                        s = t
                # lets use average similarity to see what will happen
                t = 0
                for training_ref_seq in key_dict[key_ref]: 
                    t += aligner.score(Seq(testing_ref_seq, IUPAC.protein), Seq(training_ref_seq, IUPAC.protein))
                s = round(t/len(key_dict[key_ref]), 3)
                # load
                ref_similarity.append([key_ref, s])
            ref_similarity.sort(key = lambda x:x[1], reverse = True)
            # now lets compare this results with the predicted value
            predicted_key_ref = []
            for pred in scores[:top_prediction]:
                predicted_key_ref.append(pred[0])
            # calculate the number of correct predictions
            
            if ref_similarity[0][0] in predicted_key_ref:
                num_correct += 1
            
                    
#            predictions.append([scores[:top_prediction], core[self.ref_ind], n])
#            self.correct_rever_pred = num_correct# gives the percentage
        self.reverse_correct_rate = round(num_correct/len(testing), 2)
#        return predictions 
    def Intra_Inter_cluster_similarity(self, top_n_cluster):
        # creat an empty container to contain the returned results
        intra_inter_cluster_similarity = []
        clustered_complementary = self.clustered_complementary
        # Sort the clustered_complementary again, to make sure there is no bug
        for ref in clustered_complementary:
            clustered_complementary[ref].sort(key = lambda x:len(x), reverse = True)
        # Creat an empty conatainer 'cluster_container' to contain all the clusters
        # used to calculate the intra_inter_cluster_similarity
        cluster_container = []
        for ref_key in clustered_complementary:
            # use the top n_cluster to represent the characteristic of the ref
            temp = []
            cluster_range = min(top_n_cluster, len(clustered_complementary[ref_key]))
            for n in range(cluster_range):# n_cluster is parameter
                temp.extend(clustered_complementary[ref_key][n])
            if len(temp) >= 2:# make sure we are able to calculate the intra cluster similarities
                cluster_container.append([ref_key, temp])
        # Now the data are prepared, lets get to the work of calculating similarities
        for ref_key_cluster in cluster_container:
            # calculate intra similarity, we use the mean similarity
            t = 0
            n = 0
            for seq1 in range(len(ref_key_cluster[1])-1):      
                for seq2 in range(seq1, len(ref_key_cluster[1])):
                    n += 1
                    t += aligner.score(Seq(ref_key_cluster[1][seq1], IUPAC.protein), Seq(ref_key_cluster[1][seq2], IUPAC.protein))
            intra_similarity = round(t/n, 3)
            # For inter similarity, we use the mean similarity as well.
            inter_similarity = []
            for ref_key_cluster_1 in cluster_container:
                if ref_key_cluster_1 != ref_key_cluster:
                    s = 0
                    n = 0
                    for i in ref_key_cluster_1[1]:
                        for j in ref_key_cluster[1]:
                            n += 1
                            s += aligner.score(Seq(i, IUPAC.protein), Seq(j, IUPAC.protein))
                    inter= round(s/n, 3)
                    inter_similarity.append(inter)# inter_similarity contains all the inter_cluster similarities
            # load the data
            inter_mean = round(sum(inter_similarity)/ len(inter_similarity), 3)
            intra_inter_cluster_similarity.append([ref_key_cluster[0],intra_similarity, inter_mean])
        # Pass the data to self
        self.intra_inter_cluster_similarity = intra_inter_cluster_similarity
        
            

##############################################################################  

#############################################################################

analysis = Complementary_cluster(training, testing_data = testing, ref_chain='Ag')
analysis.Convert_data()
analysis.data_single_letter
analysis.Cluster_ref(highest_score=None, cluster_number=20)
analysis.Groupup_to_ref_cluster()
analysis.key_dict.keys()
analysis.key_dict['2']
len(analysis.key_dict['5'])
len(analysis.clustered_ref)
len(analysis.Grouped_dict)
analysis.Grouped_dict['5']
#n = 0
#for i in analysis.key_dict:
#    n += len(analysis.key_dict[i])
#    print(len(analysis.key_dict[i]))
#n
#####
analysis.Frequency(analysis.Grouped_dict)
analysis.frequency.sort(key = lambda x:x[1], reverse = True)
analysis.frequency
cluster_parameter_dict = {}
for i in analysis.frequency:
    cluster_parameter_dict[i[0]] = [10, None]
len(cluster_parameter_dict ) 
analysis.cluster_para = cluster_parameter_dict 
analysis.Cluster_complementary() 
for i in analysis.clustered_complementary:
    print(len(analysis.clustered_complementary[i]))
analysis.complementary_cluster_score
analysis.clustered_complementary.keys()
analysis.clustered_complementary['2']
len(analysis.clustered_complementary['4'])
analysis.Top_cluster_percentage(2)
analysis.top_n_cluster_percentage
analysis.clustered_complementary['4'].sort(key = lambda x:len(x), reverse = True) 
for i in analysis.clustered_complementary['4']:
    print(len(i))
analysis.Processing_testing_data()
analysis.processed_testing_data
analysis.Prediction()
analysis.all_predictions 
analysis.Correct_rate(2)
analysis.correct_rate  
predictions = analysis.Reverse_prediction(4, 1)
analysis.reverse_correct_rate    
analysis.Intra_Inter_cluster_similarity(1)
analysis.intra_inter_cluster_similarity  
#################################################################################       
   
#################################################################################
       

   


