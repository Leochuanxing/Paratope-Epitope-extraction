import os
os.getcwd()
os.chdir("C:\\Users\\leo\\Documents\\Research\\Database\\Biopython")

import json
with open('seq_homo', 'r') as f:
    seq_homo = json.load(f)
with open('id_dict_4A', 'r') as f:
    id_dict_4A = json.load(f)
#################################################################################
# Creat functions to do conversion
'''
Inputs: sequence, a list of sequence in terms of triple-letter notation
Returns: single_letter_seq, a string sequence in terms of single-letter notation
         corresponding to the input sequence.
'''

def triple_to_single(sequence):
    TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
       ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
       ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]
    single_letter_seq = ''
    for i in sequence:
        for j in TripleSingle:
            if j[0] == i:
                single_letter_seq += j[1]           
    return single_letter_seq

# Define a tiny fuction to swich sequence to Seq object
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
def To_Seq(seq):
    Seq_obj = Seq(seq, IUPAC.protein)  
    return Seq_obj
'''
Inputs: iddict, a list if ids in the form of [[H, L, A], [D, C, G],...]
        seq, a dictionary of sequence, corresponding to the same pdb as iddict
Returns: CDR_seq, a list in the form of [['CDRL', L, [SEQ], [SEQ], [SEQ]], ['CDRL', B, [], [], []], ...]
         where the SEQ are in single-letter format, and all the SEQ are from antibody chain that contact 
         with the antigens
'''
def CDR_Seq_Ob(iddict, seq):
    CDR_seq = []
    for i in iddict:
        if i[2] != '':
            #we extract the triple_letter sequence from the seq_homo
            if i[0] != '':
                CDRH = [triple_to_single(seq[i[0]][25:36]),triple_to_single(seq[i[0]][46:65]),
                        triple_to_single(seq[i[0]][90:110])]
                CDR_seq.append(['CDRH', i[0], To_Seq(CDRH[0]),To_Seq(CDRH[1]), To_Seq(CDRH[2])])
            if i[1] != '':
                CDRL = [triple_to_single(seq[i[1]][23:36]), triple_to_single(seq[i[1]][45:56]),
                        triple_to_single(seq[i[1]][88:97])]
                CDR_seq.append(['CDRL', i[1], To_Seq(CDRL[0]),To_Seq(CDRL[1]), To_Seq(CDRL[2])])
    return CDR_seq

'''
Inputs: iddict_all, a dictionary in the form of {'1adq': [['H','L','A']...],...}
         seq_all, a dictionary in the form of {'1adq':{'H':[ALA, ARG,....]....}}
Returns: Seq_dic, a dictionary in the form of {'1adq': [['CDRL', L, [SEQ], [SEQ], [SEQ]],...],...}
          and each SEQ is an Seq object.
'''
def Seq_all(iddict_all, seq_all):
    Seq_dic = {}
    for i in iddict_all:
        Seq_dic[i] = CDR_Seq_Ob(iddict_all[i], seq_all[i])
    return Seq_dic

# Check if it works
Seq_dic = Seq_all(id_dict_4A, seq_homo)
len(Seq_dic)
Seq_dic.keys()
Seq_dic['1adq'] 
Seq_dic['5u3k']   
id_dict_4A['1n0x']   
####################################################################################
# Find a proper upper limit of alignment scores
# Define an aligner
from Bio import Align
aligner = Align.PairwiseAligner()
aligner.open_gap_score = 0
aligner.extend_gap_score = 0
aligner.match = 1
aligner.mismatch = 0
'''
Inputs: list1 and list2 are in the form of  ['4ydl','H',Seq1, Seq2, Seq3]
Return: score, the alignment score of the two Seq objects
'''
def score(list1, list2):
    score = 0
    for i in range(2, 5):
        score += aligner.score(list1[i], list2[i])
    return score
# Test this function, the results should be 9
a = [ 1, 1, To_Seq('ATWY'), To_Seq('WRSV'), To_Seq('DEYF')]
b = [1 , 1, To_Seq('AT'), To_Seq('WRVS'), To_Seq('DEYTF')]
score(a,b)
# Separate all the Seq objects into to groups, CDRL and CDRH
'''
Inputs: Seq_dic, a dictionary, of the same form of the returns of Seq_all function
Returns: HL_seq, a dictionary with keys 'CDRH', and 'CDRL', the values are in the form 
         of [['1adq', 'L', Seq1, Seq2, Seq3], ...], Seq1, Seq2, Seq3 are Seq objects.
         'L' is the chain id of either heavy chain or light chain indicated by the CDRL or 
         CDRH
'''
def HL_group(Seq_dic):
    HL_seq = {}
    HL_seq['CDRH'] = [] 
    HL_seq['CDRL'] = [] 
    
    for i in Seq_dic:
        for j in Seq_dic[i]:
            if j[0] == 'CDRH':
                HL_seq['CDRH'].append([i, j[1], j[2], j[3], j[4]])
            if j[0] == 'CDRL':
                HL_seq['CDRL'].append([i, j[1], j[2], j[3], j[4]])
    return HL_seq
# Check if it works
HL_Seq = HL_group(Seq_dic)
len(HL_Seq['CDRL'])
len(HL_Seq['CDRH'])     
HL_Seq['CDRL'][:5] 
HL_Seq['CDRH'][:5]

####################################################
## calculate all the alignment scores
#def All_scores(HL_Seq):
#    HL_all_scores ={}
#    HL_all_scores['CDRH'] = []
#    HL_all_scores['CDRL'] = []
#    length = len(HL_Seq['CDRH'])
#    for i in range(length-1):
#        for j in range(i+1, length):
#            HL_all_scores['CDRH'].append(score(HL_Seq['CDRH'][i], HL_Seq['CDRH'][j]))
#          
#    length = len(HL_Seq['CDRL'])
#    for i in range(length-1):
#        for j in range(i+1, length):
#            HL_all_scores['CDRL'].append(score(HL_Seq['CDRL'][i], HL_Seq['CDRL'][j]))
#    return HL_all_scores
##Check if it works
#HL_all_scores = All_scores(HL_Seq)
#len(HL_all_scores['CDRH'])
#len(HL_all_scores['CDRL'])
#HL_all_scores['CDRH'][:5]
## plot a histogram to show the distribution of scores
#def plot_hist(scores):
#    import matplotlib.pyplot as plt
#    from math import ceil
#    u = ceil(max(scores)) +1
#    l = ceil(min(scores)) - 1
#    bins = list(range(l, u))
#    plt.hist(scores, bins)
#    return plt.show
## Check if it works
#plot_hist(HL_all_scores['CDRH']) 
#plot_hist(HL_all_scores['CDRL']) 
#'''The total length of CDRH is 50, and the CDRL is 35， we do a simple interpretation
#  that less than 80 percent of CDRH seq are identical means the score is less than
#  40
#'''   
#HL_all_scores['CDRH'].count(50)
##########################################################################
# select from the HL_Seq elements with alignment scores under the threshold
#selected_HL_Seq = {}
#selected_HL_Seq['CDRH'] = []
#selected_HL_Seq['CDRL'] = []
#
#''' Don't simple select, need to do a cluster'''
#def select(Seq, selected_Seq, bound):
#    if Seq != []:
#        Seq.sort(key = lambda x:x[0])
#        Seq0 = Seq[0]
#        selected_Seq.append(Seq0)
#        Seq.remove(Seq0)
#        if Seq != []:
#            for i in Seq:
#                if score(Seq0, i) > bound:
#                    Seq.remove(i)
#        return select(Seq, selected_Seq, bound)
#    if Seq == []:
#        return selected_Seq
## Check if it works
#                      
#selected_H_Seq = select(HL_Seq['CDRH'], [], bound = 45)   
#len(selected_H_Seq) 
#selected_L_Seq = select(HL_Seq['CDRL'], [], bound = 31)
#len(selected_L_Seq)     
# Change the results to dictionary    
###########################################################
# lets do a dendrogra
# creat a distance matrix
import numpy as np
# sort the HL_Seq, according to the pdb id and chain id
H_Seq = HL_Seq['CDRH']
H_Seq.sort(key = lambda x: x[1])
H_Seq.sort(key = lambda x:x[0])
H_Seq[:10]




def Similarity_matrix(H_Seq):
    similarity_matrix = np.zeros([len(H_Seq), len(H_Seq)])
    for i in range(len(H_Seq)):
        for j in range(len(H_Seq)):
            similarity_matrix[i,j] = score(H_Seq[i], H_Seq[j])
    return similarity_matrix
# caltulate the results
similarity_matrix = Similarity_matrix(H_Seq) 
# take a look at the results           
similarity_matrix[:5, :5]
score(H_Seq[0], H_Seq[3])

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
        score = 0
        pos = -1
        for j in range(i+1, len(clusters)):
            # calculate the distance between cluster i and cluster j
            # here we use the smallest similarity score
            # creat a variable to contain this smalles score
            similarity = 51           
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

# check this function 
# here starters, means the begining, when each one is a cluster
# and the ids for each H_seq element are give in numbers
starters = []
for i in range(len(H_Seq)):
    starters.append([i])
cluster_dist = Cluster_dist(starters, similarity_matrix)
starters[:10]
cluster_dist[:10]
len(cluster_dist)
H_Seq[7]
H_Seq[23]

# define a sub function
'''
inputs: need_merge, a list, for example: [ [1,2],[2, 3], [5, 6], [1, 7]]
return: need_mrege, a modified list on basis of need_merge, this can be achieved 
        because list is an mutable object. The reults will be [[1,2,3,7], [5,6]]
'''
def Merge_sub(need_merge):
    for i in range(len(need_merge)-1):
        for j in range(i+1, len(need_merge)):
            intersection = [v for v in need_merge[i] if v in need_merge[j]]
            if intersection !=[]:
                need_merge[i].extend(need_merge[j])
                need_merge[i] = list(set(need_merge[i]))
                del need_merge[j]
                break
        if intersection != []:
            break
    if intersection ==[]:
        return need_merge
    return Merge_sub(need_merge)

# define a function to merge the result generated by Cluster_dist
'''
Inputs: cluster_dist, a list in the form of [[[i,j],dist], ....], the terms i,j 
        corresponding to the positions in the clusters below
        clusters, a list in the form of [[0, 1, 2], [3,5], ...], corresponding 
        to cluster_dist. 0, 1, 2,...are the positions of the items need to be 
        clustered, they can be applied directly on the similarity_matrix
return, merged_clusters, a list of the form [[0, 1, 2], [3,5], ...], where 0, 1, 2
        ... are corresponding to the items need to be clustered
        highest, a float, gives the largest similarity that to be merged in this 
        merge step.
'''
def Merge(clusters, similarity_matrix, threshold = 50):
    # find the paires with the highest similarity
    cluster_dist = Cluster_dist(clusters, similarity_matrix)
    # Generate the list need to be merged with the maximum similarity, for example
    #[[[0, 45], 37.0], [[1, 2], 50.0], [[2, 389], 34.0], [[3, 4], 50.0]]
    # we select outh the emlemts with the highest score and merge
    # select out the elements with the highest score
    highest = max([i[1] for i in cluster_dist])
    if highest >= threshold:
        merge = []
        for i in cluster_dist:
            if i[1] == highest:
                merge.append(i[0])
        # merge the positions
        Merge_sub(merge) 
        # merge according to the information given in merge
        for i in merge:
            for j in range(1, len(i)):
                clusters[i[0]].extend(clusters[i[j]])
                clusters[i[j]] = 'D' 
        # get rid of the 'D'
        while 'D' in clusters:
            clusters.remove('D')                
#        clusters = list(filter(lambda x: x != 'D', clusters))
        # sort the merged_clusters beform return
        for i in range(len(clusters)):
            clusters[i].sort()
        clusters.sort(key = lambda x:x[0])
        return Merge(clusters, similarity_matrix, threshold) 
    else:
        return clusters
    
# check if the function works                 
starters[:10]    
import copy; clusters = copy.deepcopy(starters)
merged = Merge(clusters, similarity_matrix, threshold = 49)
len(merged)
len(clusters)
clusters[:11]
merged[:11]
for i in clusters:
    n = True
    if i not in merged:
        n = False
        break
if n:
    print('Same')
else:
    print('Different')
        
cluster_dist = Cluster_dist(clusters, similarity_matrix)
highest = max([i[1] for i in cluster_dist])
highest
similarity_matrix[1, 2]
H_Seq[2]
#################################################################
# find the corresponding pdb id and the chain id for the clusters
'''
inputs: clusters, a list, for example [[0],[1, 2]]
        Seq_sorted, a list, used to calculate the similarity matrix, the order 
        should be kept a constant
return: id_clusters, a list, of clustered id, for example[['1adqH'], ['1bvkB',''1bvkE' ]]
'''
def Id_clusters(clusters, Seq_sorted):
    # Creat an empty container.
    id_clusters = []
        # Load data
        
    for i in range(len(clusters)):
        #Creat a temporary container
        temp = []
        for j in clusters[i]:
            temp.append(Seq_sorted[j][0]+Seq_sorted[j][1])
        id_clusters.append(temp)
    return id_clusters

# Check the function
id_clusters = Id_clusters(clusters, H_Seq)
id_clusters[:3]
len(id_clusters)
##################################################################
# find the one with the largest contact within each cluster
'''
inputs: id_clusters, a list, for example [['1adqH'], ['1bvkB', '1bvkE']]
        contact, a dictionary, in the form of {'1adq':[four-coordinates. four-coordinates]}
        which is the output of the AAC module.
return: ar_contact, a dictionary, in the form of {'1adqH':[four-coordinates. four-coordinates]}
        for each cluster chose the one with the largest contact number
'''
def AR_contact(id_clusters, contact):# AR means, alignment restraint contact
    # Creat an empty dictionary
    ar_contact ={}
    for cluster in id_clusters:
        # Creat a temporary dictionary to record the total contact number
        temp_dict = {}
        for ids in cluster:
            temp_dict[ids] = 0
            if ids[:4] in contact:
                for fcdn in contact[ids[:4]]:
                    if fcdn[0][2] == ids[4]:
                        temp_dict[ids] += fcdn[3]
        # Chose the one with the largest contact number from temp_dict
        # creat a temporary key hoder and a contact number holder
        temp_key =''
        n = 0
        for key in temp_dict:
            if temp_dict[key] > n:
                temp_key = key
                n = temp_dict[key]
        # Load data to ar_contact
        if temp_key != '':
            ar_contact[temp_key] = []
            for fcdn in contact[temp_key[:4]]:
                if fcdn[0][2] == temp_key[4]:
                    ar_contact[temp_key].append(fcdn)
        # clear the temporary containers
    del temp_key, n, temp_dict
        
    return ar_contact

# Check if it works
with open('contact_homo_4A', 'r') as f:
    contact = json.load(f)
ar_contact = AR_contact(id_clusters, contact)
len(ar_contact)
len(clusters)
keys = list(ar_contact.keys())
keys[:10]
