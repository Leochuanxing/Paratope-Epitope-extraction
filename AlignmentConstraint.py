import os
os.getcwd()
os.chdir("/home/leo/Documents/Database/Pipeline/Homo")

import json
with open('sequence', 'r') as f:
    sequence = json.load(f)
with open('good_matched_ids', 'r') as f:
    good_matched_ids = json.load(f)
with open('dud_AAC', 'r') as f:
    dud_AAC = json.load(f)
with open('contact', 'r') as f:
    contact = json.load(f)
# Get rid of the dud
for dud in dud_AAC:
    del good_matched_ids[dud]
good_matched_ids.keys()
good_matched_ids['4om1']
contact.keys()
contact['1adq']
sequence['1adq'].keys()
len(sequence['1adq']['H'])
#################################################################################
# Creat functions to do conversion
'''
triple_to_single:
    to convert a sequence of triple letter notation to single letter notation
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
'''
To_Seq:
    is to convert the aa sequence into aa sequence object
Input:
    seq, a sequence of single-letter notation
Output:
    Seq_obj, a sequence object
'''
def To_Seq(seq):
    for aa in seq:       
        Seq_obj = Seq(seq, IUPAC.protein)  
    return Seq_obj
'''
CDR_Seq_Ob:
    to convert the CDR sequence into sequence object.
Inputs: iddict, a list if ids in the form of [[H, L, A], [D, C, G],...]
        seq, a dictionary of sequence, corresponding to the same pdb as iddict
Returns: CDR_seq, a list in the form of [['CDRL', L, [SEQ], [SEQ], [SEQ]], ['CDRL', B, [], [], []], ...]
         where the SEQ are in single-letter format, and all the SEQ are from antibody chain that contact 
         with the antigens. L and B are the names of the antibody chains.
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
Seq_all:
    to convert all the CDRs of all the complexes into sequence objects.
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

'''Check if it works'''
Seq_dic = Seq_all(good_matched_ids, sequence)
len(Seq_dic)
keys = list(Seq_dic.keys())
Seq_dic[keys[8]]  
####################################################################################
# Find a proper upper limit of alignment scores
# Define an aligner

'''referring to he hamming distance'''

from Bio import Align
aligner = Align.PairwiseAligner()
aligner.open_gap_score = 0
aligner.extend_gap_score = 0
aligner.match = 1
aligner.mismatch = 0
'''
score:
    find the alignment score for a given pair of CDR sequences under the defined aligner.
    Each CDR contains 3 sequences.
Inputs: 
    list1 and list2 are in the form of  ['4ydl','H',Seq1, Seq2, Seq3]
Return: 
    score, the alignment score of the two Seq objects
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
# Separate all the Seq objects into to groups, CDRL and CDRH. The motivation to do 
# this is we want align heavy chains with heavy chains and light chains with light 
# chains separately.
'''
HL_group:
    Separate all the Seq objects into to groups, CDRL and CDRH. The motivation to do 
    this is we want align heavy chains with heavy chains and light chains with light 
    chains separately.
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
#HL_Seq = HL_group(Seq_dic)
#len(HL_Seq['CDRL'])
#len(HL_Seq['CDRH'])     
#HL_Seq['CDRL'][:5] 
#HL_Seq['CDRH'][:5]
####################################################
   
###########################################################
'''
Sort:
    to sort the values in 'HL_seq' according to the pdb id and the chain id
Input:
    HL_seq, the return of fucntion HL_group
Output:
    the sorted HL_seq
'''
# sort the HL_Seq, according to the pdb id and chain id
def Sort(HL_Seq):
    for Ab_chain in HL_Seq:
        HL_Seq[Ab_chain].sort(key = lambda x: x[1])
        HL_Seq[Ab_chain].sort(key = lambda x:x[0])
    return

'''
Similarity_matrix:
    to give the similarity scores between each pair of sequences for a given set 
    of sequences
Input:
    H_seq, a list of sequence object
Output:
    A symmetric matrix, gives the similarity scores between any pair of sequence obj
    in H_seq.
'''
def Similarity_matrix(H_Seq):
    similarity_matrix = np.zeros([len(H_Seq), len(H_Seq)])
    for i in range(len(H_Seq)):
        for j in range(len(H_Seq)):
            similarity_matrix[i,j] = score(H_Seq[i], H_Seq[j])
    return similarity_matrix
# caltulate the results
#similarity_matrix = Similarity_matrix(H_Seq) 
## take a look at the results           
#similarity_matrix[:10, :10]
#score(H_Seq[0], H_Seq[3])
################################################################
# define a function to caltulate the distances between clusters
'''
Cluster_dist:
        for a fixed cluster i, find cluster j with the largest similarity between i and j.
        the similarity between clusters is defined as the smallest similarities between elements
        in those two clusters.
inputs: clusters, a list in the form of [[0, 1, 2], [3,5], ...]
        [0, 1, 2] is the first cluster.
        dist_matrix, a distance matrix, returns from Dist_matrix
return: cluster_dist, a list in the form of [[[i,j],dist], ....]
        i, j, means the i th and the j th cluster, dist is the distance between 
        those two clusters, here we use the least similarities, and the results 
        given are the paires with the highest similarity
        [[i,j],dist], if we fix i, j is the cluster with the highest similarity 
        with i, and the similarity between i and j is 'dist'.  The similarity between
        i and j is calculated as the minimum similarity between all possible sequences in 
        those clusters
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
            # creat a variable to contain this smallest score
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


# define a function to merge the result generated by Cluster_dist
'''
Merge:
    to merge clusters with similarities larger than the highest_score
Inputs: clusters, a list in the form of [[0, 1, 2], [3,5], ...], The same as in
           Cluster_dist. 0, 1, 2,...are the positions of the items need to be 
                clustered, they can be applied directly on the similarity_matrix
        highest_score, a number, gives the nearest upper limit of the similarity score
               #Only one of cluster_number and highest_score can be used to control 
               the clustering process.

return, merged_clusters, a list of the form [[0, 1, 2], [3,5], ...], where 0, 1, 2
        ... are corresponding to the items need to be clustered
        highest, a float, gives the largest similarity that to be merged in this 
        merge step.
'''
def Merge(clusters, similarity_matrix,  highest_score ):
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
    
    # go to the mode of highest_score

    if highest >= highest_score:                    
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
            return clusters
        else:
            return Merge(clusters, similarity_matrix, highest_score ) 
    else:
        return clusters
# check if the function works                 
#starters[:10]    
#import copy; clusters = copy.deepcopy(starters)
#merged = Merge(clusters, similarity_matrix, threshold = 45)
#len(merged)
#merged[:6]
#       
#cluster_dist = Cluster_dist(clusters, similarity_matrix)
#highest = max([i[1] for i in cluster_dist])
#highest
#################################################################
# find the corresponding pdb id and the chain id for the clusters
'''
Id_clusters:
    to translate the clusters from the function Merge to their corresponding Ids.
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
#id_clusters = Id_clusters(clusters, H_Seq)
#id_clusters[:6]
#len(id_clusters)
##################################################################
# find the one with the largest contact within each cluster
'''
AC_contact:
    to find the one with the largest contact within each cluster
inputs: id_clusters, a list, for example [['1adqH'], ['1bvkB', '1bvkE']]
        contact, a dictionary, in the form of {'1adq':[four-coordinates. four-coordinates]}
        which is the output of the AAC module.
        contact, the output of AAC_2
        hl, if it takes value h, it means heavy chain, and l means light chain.
return: ac_contact, a dictionary, in the form of {'1adqH':[four-coordinates. four-coordinates]}
        which is the one with the largest contact number.
'''
def AC_contact(id_clusters, contact, hl = 'h'):# AC means, alignment constraint contact
    # Creat an empty dictionary
    ac_contact ={}
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
        # Load data to ac_contact
        if temp_key != '':
            ac_contact[temp_key+hl] = []
            for fcdn in contact[temp_key[:4]]:
                if fcdn[0][2] == temp_key[4]:
                    ac_contact[temp_key+hl].append(fcdn)
        # clear the temporary containers
    del temp_key, n, temp_dict
        
    return ac_contact

# Check if it works
#with open('contact', 'r') as f:
#    contact = json.load(f)
#ac_contact = AC_contact(id_clusters, contact)
#len(ac_contact)
#len(id_clusters)
#ac_contact.keys()
#keys = list(ac_contact.keys())
#ac_contact[keys[6]]
#################################################
# Creat a main function to combine the above functions together, and deals with
# both the heavy chains and light chains
import numpy as np
def main(good_matched_ids, sequence, contact, threshold = [45, 31]):
    # the returned results stored in ac_contact
    # creat an empty container
    ac_contact = {}
    # Convert the CDR sequences into Seq objects
    Seq_dic = Seq_all(good_matched_ids, sequence)
    # separate Seq_dic into two groups, and store them in HL_Seq, the keys are
    # 'CDRL' and 'CDRH'
    HL_Seq = HL_group(Seq_dic)
    # sort the values in HL_Seq
    Sort(HL_Seq)
    # load the threshhold parameters
    threshold_dict = {}
    threshold_dict['CDRH'] = threshold[0]
    threshold_dict['CDRL'] = threshold[1]
    # calculate the similarity_matrix, and store the results in similarity_matrix
    # with keys 'CDRL' and 'CDRH'.And cluster.
    similarity_matrix = {}
    for CDR in HL_Seq:        
        similarity_matrix[CDR] = Similarity_matrix(HL_Seq[CDR])
        # creat starters as the begins of clustering
        starters = []
        for i in range(len(HL_Seq[CDR])):
            starters.append([i])
        # merge the indices.
        merged = Merge(starters, similarity_matrix[CDR], threshold_dict[CDR])
        # map the merged indices into ids.
        id_clusters = Id_clusters(merged, HL_Seq[CDR])
        # Find the one with the largest contact in each cluster
        ac_contact_sub = AC_contact(id_clusters, contact, hl = CDR[3].lower())
        # merge into the final result
        for Ab_chain in ac_contact_sub:
            ac_contact[Ab_chain] = ac_contact_sub[Ab_chain]
    return ac_contact, similarity_matrix, id_clusters, HL_Seq, 
# check if it works
ac_contact, similarity_matrix, id_clusters, HL_Seq = main(
                         good_matched_ids, sequence, contact, threshold = [40, 28])
len(ac_contact)#486
keys = list(ac_contact.keys())
keys
ac_contact['5ucbLl']
'''
It is better we include the pdbid into the ac_contact
'''
###############################################
## a simple verification of the results
#id_clusters[:10]
#HL_Seq['CDRL'][:10]
##Find the corresponding indices for the id_clusters according to HL_Seq
##Creat a list of corresponding Ab_chains
#Ab_chains = [i[0]+i[1] for i in HL_Seq['CDRL']]
#Ab_chains[:5]
#index = []
#for cluster in id_clusters:
#    temp = []
#    for e in cluster:
#        temp.append(Ab_chains.index(e))
#    index.append(temp)
#index[:5]        
#similarity_matrix['CDRL'][:5, :5] 
#cluster_dist = Cluster_dist(index, similarity_matrix['CDRL'])
#dist = 0
#for i in cluster_dist:
#    if i[1] > dist:
#        dist = i[1]
#dist
#################################################
#with open('ac_contact', 'w') as f:
#    json.dump(ac_contact, f)
with open('ac_contact', 'r') as f:
    ac_contact = json.load(f)
len(ac_contact)
#########################################################################
'''
triple_to_single:
    to convert a sequence of triple letter notation to single letter notation
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
'''
To_Seq:
    is to convert the aa sequence into aa sequence object
Input:
    seq, a sequence of single-letter notation
Output:
    Seq_obj, a sequence object
'''
def To_Seq(seq):
    if len(seq[0]) == 1:
            seq_single_letter = []
        for aa in seq:  
            aa 
        Seq_obj = Seq(seq, IUPAC.protein)  
        return Seq_obj
##########################################################################
'''
Write the above functions into a class, named 'AlignmentConstraint'

'''
class AlignmentConstraint(object):
    def __init__(self, iddict, sequence, contact):
        self.iddict = iddict
        self.sequence = sequence
        self.contact = contact
        self.TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
       ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
       ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]
        
        pass
    
    def To_seq(self, aa_sequence):
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC
        
        seq_obj = None

        seq_single_letter = ''
        for aa in aa_sequence:
            for TS in self.TripleSingle:
                if TS[0] == aa:
                    seq_single_letter += TS[1]
        seq_obj = Seq(seq_single_letter, IUPAC.protein)
        
        return seq_obj
    '''
    SeqCDR:
        to convert the sequence into sequence object and be prepared for futher alignment
    Input:
        The self object
    Output:
        in the form of a list [['1abvhH', Seq_obj], ['1abvlL', Seq_obj]]
        ['1abvhH', Seq_obj] means the pdbid is 1abv heavy chain with chain name H, and 
        the Seq_obj is the sequence object of the concatenated CDR1, CDR2, CDR3 of chain H.
    '''   
    def SeqCDR(self):
        # Concatenate the sequences 
        seqCDRH = []
        seqCDRL = []
        for pdbid in self.iddict:
            for combination in self.iddict[pdbid]:
                Concatenated_CDRH = []
                Concatenated_CDRL = []
                if combination[2] != '':
                    if combination[0] != '':
                        CDRH = [pdbid+'h'+combination[0]]
                        CDRH.append(self.To_seq(self.sequence[pdbid][combination[0]][25:36]))
                        CDRH.append(self.To_seq(self.sequence[pdbid][combination[0]][46:65]))
                        CDRH.append(self.To_seq(self.sequence[pdbid][combination[0]][90:110]))
                        seqCDRH.append(CDRH)
                    if combination[1] != '':
                        CDRL = [pdbid+'l'+combination[1]]
                        CDRL.append(self.To_seq(self.sequence[pdbid][combination[1]][23:36]))
                        CDRL.append(self.To_seq(self.sequence[pdbid][combination[1]][45:56]))
                        CDRL.append(self.To_seq(self.sequence[pdbid][combination[1]][88:97]))
                        seqCDRL.append(CDRL)                        
        self.seqCDRH = seqCDRH
        self.seqCDRL = seqCDRL
        
    '''
    Align according to the given criteria and choose the none redundent one
    '''
    
    def Hamming_like_dist (u, v):
        # 1 creat the similarity matrix
        # 2 Do hyerachical clustering
        # 3 Chose the representative
        dist = 0 
        from Bio import Align
        aligner = Align.PairwiseAligner()
        aligner.open_gap_score = -1 # it is unlikly that two sequence of the same source could be different by a gap
        aligner.extend_gap_score = 0
        aligner.match = 1
        aligner.mismatch = 0 
        aligner.mode = 'local'
        
        dist += aligner.score(u[1], v[1])
        dist += aligner.score(u[2], v[2])
        dist += aligner.score(u[3], v[3])
        l1 = len(u[1]) + len(u[2]) + len(u[3])
        l2 = len(v[1]) + len(v[2]) + len(v[3])
        l = min(l1, l2)
        dist =  1 - dist /l
        return dist
    
    def Hcluster(self):

        from Bio import Align
        aligner = Align.PairwiseAligner()
        aligner.open_gap_score = -1 # it is unlikly that two sequence of the same source could be different by a gap
        aligner.extend_gap_score = -1
        aligner.match = 1
        aligner.mismatch = 0 
        aligner.mode = 'local'
        
        from scipy.spatial.distance import pdist
        from scipy.spatial.distance import squareform
        from scipy.cluster.hierarchy import linkage
        
        CDR = [self.seqCDRH, self.seqCDRL]
        
        for idx in range(2):
            
            toydata = CDR[idx]
            distmatrix = np.zeros((len(toydata), len(toydata)))
            for i in range(len(toydata)):
                for j in range(len(toydata)):
                    distmatrix[i,j] += aligner.score(toydata[i][1], toydata[j][1])
                    distmatrix[i,j] += aligner.score(toydata[i][2], toydata[j][2])
                    distmatrix[i,j] += aligner.score(toydata[i][3], toydata[j][3])
                    l1 = len(toydata[i][1]) + len(toydata[i][2]) + len(toydata[i][3])
                    l2 = len(toydata[j][1]) + len(toydata[j][2]) + len(toydata[j][3])
                    l = min(l1, l2)
                    distmatrix[i,j] =  1 - distmatrix[i,j]/l
                    
            if idx == 0:
                self.distmatrix_CDRH = distmatrix
                self.hcluster_CDRH = linkage(squareform(distmatrix), method = 'complete')
            else:
                self.distmatrix_CDRL = distmatrix
                self.hcluster_CDRL = linkage(squareform(distmatrix), method = 'complete')
        


                
        
 



alignment =  AlignmentConstraint(good_matched_ids, sequence, contact)
alignment.SeqCDR()
alignment.Hcluster()


from Bio import Align
aligner = Align.PairwiseAligner()
aligner.open_gap_score = -1
aligner.extend_gap_score = -1
aligner.match = 1
aligner.mismatch = 0 
aligner.mode = 'local'

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
s1 = Seq('AYRWED', IUPAC.protein)
s2 = Seq('RWEDA', IUPAC.protein)
aligner.score(s1, s2)
alignments = aligner.align(s1, s2)
# creat a distance matrix
toydata = alignment.seqCDRH[:]
import numpy as np
distmatrix = np.zeros((len(toydata), len(toydata)))
for i in range(len(toydata)):
    for j in range(len(toydata)):
        distmatrix[i,j] += aligner.score(toydata[i][1], toydata[j][1])
        distmatrix[i,j] += aligner.score(toydata[i][2], toydata[j][2])
        distmatrix[i,j] += aligner.score(toydata[i][3], toydata[j][3])
        l1 = len(toydata[i][1]) + len(toydata[i][2]) + len(toydata[i][3])
        l2 = len(toydata[j][1]) + len(toydata[j][2]) + len(toydata[j][3])
        l = min(l1, l2)
        distmatrix[i,j] =  1 - distmatrix[i,j]/l
np.shape(distmatrix)
distmatrix[:6, :6] 
distmatrix[543, 543]   
alignment.distmatrix_CDRH[:6, :6]   

        



                                            

'''
'''


'''
Hcluster: A function to perform hierarchical clustering
Input: distance_matrix, a symmetric matrix, with the main diagal 0
Output: tree, in the form of an array, which could be used to show the 
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

    

import pandas as pd
df = pd.DataFrame(distmatrix)
df.to_csv("distmatrix.csv")

from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
z = np.array([[0, 1, 0.5, 2],
              [2, 3, 0.5, 2],
              [4, 5, 1, 4]])
            
plt.figure(figsize = (10, 5))
plt.title('Test graph')
plt.xlabel('Sample index')
plt.ylabel('distance')
dendrogram(
        z,
        leaf_rotation = 90,
        leaf_font_size = 8     
        )
plt.show()
z= alignment.hcluster_CDRH
np.shape(alignment.distmatrix_CDRH)
