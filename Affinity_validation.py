import copy
import numpy as np
import json
import os
os.chdir('/home/leo/Documents/Database/Pipeline/Codes and testing data')
from AAC_2 import Coordinates
from AAC_2 import Get_contact, Chain_seq
from FrameConstraint import Add, Sub_seq, Get_consecutive

from RBFN_coverage import To_seq
from RBFN_coverage import Multiplication_distance



##########################################################

from Bio import Align
aligner = Align.PairwiseAligner()
aligner.open_gap_score = -1 # it is unlikly that two sequence of the same source could be different by a gap
aligner.extend_gap_score = -1
aligner.match = 1
aligner.mismatch = -1 
aligner.mode = 'global'
#####################################################################


###################################################################

'''
Criteria:
    First we pick the number of mutations positions, they should be continuous in the 1_free sense.
    
    Second, we find all the opposite consecutive antigen amino acids in the 1_free sense.
             Among all those consecutive antigen amino acids with length no more than 4, we 
             choose the one with the largest contact numbers with the given mutation positions.   
'''
###################################################################################
      
'''
Max_contact:
    returns the nple with the most contact
Input:
    consecutive:
        a list of list of consecutives with the same length. It shouldn't be empty
    selected_contact:
        a sub set of chian_contact, selected by the position
    Ab_Ag:
        A string, takes the value of either 'Ab' or 'Ag'
Output:
    max_nple:
        a list of consecutives with the maximum contact number
'''

def Max_contact(consecutive, selected_contact, Ab_Ag):
    max_contact = 0
    max_nple = []
    if Ab_Ag == 'Ab':
        chain_pos = 2
    elif Ab_Ag == 'Ag':
        chain_pos = 3
    for nple in consecutive:
        n_contact = 0
        for i in selected_contact:
            if i[chain_pos] in nple:
                n_contact += i[3]
        if n_contact >= max_contact:
            max_contact = n_contact
            max_nple = nple
    return max_nple
############################################################################
'''
Order_Ag_sequence:
    Change the order of the paires, so that they match
Input:
    Ab_sequence:
        a list of numbers, give the positions of the Ab. 
        #### We always suppose the Ab_sequence is ordered from small to large and fixed###
    Ag_sequence:
        a list of numbers, gives the positions of the Ag
        #### Ag_sequence is also arranged from small to large, but it may be reversed 
             according to the criteria.#####
    contact:
        the four coordinates contact information about the amino acids between Ab_sequence and
        Ag_sequence
Output:
   Ag_sequence:  
       a list of numbers with either the same order as Ag_sequence or the reversed order.
'''
def Order_Ag_sequence(Ab_sequence, Ag_sequence, contact):
    start_para = Ab_sequence[0]
    end_para = Ab_sequence[-1]
    start_epi = Ag_sequence[0]
    end_epi = Ag_sequence[-1]
    # Use the same method as what we do in the FrameConstraint
    keep = 0
    rever = 0
    for i in contact:
        if i[1]==start_para and i[2] == start_epi:
            keep += 1
        elif i[1] == end_para and i[2] == end_epi:
            keep += 1
        elif i[1] == start_para and i[2] == end_epi:
            rever += 1
        elif i[1] == end_para and i[2] == start_epi:
            rever += 1
    if rever > keep:
        Ag_sequence.sort(reverse=True)
    return Ag_sequence
####################################################################
'''
Select_contact_opposite:
    a function to select the contact from among all the contact
    and the all the possible opposite
Input:
    mutation_match_parameter:
        a dictionary, contains
        mutation_match_parameter['pdbid']
        mutation_match_parameter['mutation_chain']
        mutation_match_parameter['mutations']
        mutation_match_parameter['matched_ids']
        mutation_match_parameter['combined_ids']
        mutation_match_parameter['opposite_chain']
        mutation_match_parameter['Ab_Ag']: gives the information of whether the 
                                           mutation chain is an Ab chain or an Ag chain.
'''
def Select_contact_opposite(mutation_match_parameter):
    # take out the parameter
    pdbid = mutation_match_parameter['pdbid']
    chain = mutation_match_parameter['mutation_chain']
    mutations = mutation_match_parameter['mutations']
    matched_ids = mutation_match_parameter['matched_ids']
    combined_ids = mutation_match_parameter['combined_ids']
    opposite_chain = mutation_match_parameter['opposite_chain']
    Ab_Ag = mutation_match_parameter['Ab_Ag']
    # Extract the required data
    with open(pdbid+'.pdb', 'r') as file:
        cdn = Coordinates(file, combined_ids)
    with open(pdbid+'.pdb', 'r') as file:
        pdbid_sequence = Chain_seq(file, combined_ids)
    # store the paires in paires
    cutoff = 4
    possible_opposite = []
    if Ab_Ag == 'Ab':
        chain_pos = 2
        aa_pos = 1
    elif Ab_Ag == 'Ag':
        chain_pos = 3
        aa_pos = 2
    while possible_opposite == [] and cutoff <=8:
        contact = Get_contact(cdn, matched_ids, cutoff)
        # Carry out the above process:
        # take out all the contact containing the chain_name
        selected_contact = []
        possible_opposite = []
        for i in contact:
            if chain == i[0][chain_pos] and i[aa_pos] in mutations:                
                if Ab_Ag == 'Ag' and i[0][2] == opposite_chain:
                    selected_contact.append(i)
                    possible_opposite.append(i[3-aa_pos])
                if Ab_Ag == 'Ab':
                    selected_contact.append(i)
                    possible_opposite.append(i[3-aa_pos])                    
                
        # Increase the cutoff by 1. Make sure this is not a dead loop
        if possible_opposite == []:
            cutoff += 1
            
    return selected_contact, possible_opposite, pdbid_sequence


###################################################################
'''
Paire_select:
    a function to carry out the above criteria
Input:
    mutation_match_parameter, the same as the input of Select_contact_opposite
        
Output:
    mutation_match_parameter with one more key value
    mutation_match_parameter['pairs']:
        gives the paires of the paratope and epitope of the complex with pdb id pdbid.
        
        #### The paratope are/is corresponding to the given positions of the mutations
        #### The epitope are the ones with the longest consecutive sequences in the 1_free 
        #### sense, but not more than four amino acids
'''
def Paire_select(mutation_match_parameter):

    # Extract the required information from the pdb file
    selected_contact, possible_opposite, pdbid_sequence = Select_contact_opposite(mutation_match_parameter)
    # Take out the values from the parameters
#    pdbid = mutation_match_parameter['pdbid']
    chain = mutation_match_parameter['mutation_chain']
    mutations = mutation_match_parameter['mutations']
#    matched_ids = mutation_match_parameter['matched_ids']
#    combined_ids = mutation_match_parameter['combined_ids']
    opposite_chain = mutation_match_parameter['opposite_chain']
    Ab_Ag = mutation_match_parameter['Ab_Ag']
    # make sure we have got something and Get the longest possible consecutives 
    # The basic assumuption is that the more information we know about the surrounding 
    # environment, the more accurate our prediction. This is the reason why we choose
    # the longest possible consecutives.
    paires = []
    if possible_opposite != []:
        possible_opposite.sort()
        for length in [4, 3, 2, 1]:
            longest_possible_consecutives = Get_consecutive(possible_opposite, length, free_type=1)
            if longest_possible_consecutives != []:
                break
    # Change the order/directions
    for i in longest_possible_consecutives:
        choosen_opposite_pos = i
        # Correct the direction of the choosen_opposite_pos
        choosen_opposite_pos.sort()
        # define a small function to change the order of the paires
#        if len(mutations) >= 2 and len(choosen_opposite_pos)>=2:
        if len(choosen_opposite_pos)>=2:
            if Ab_Ag == 'Ab':
                choosen_opposite_pos = Order_Ag_sequence(mutations, choosen_opposite_pos, selected_contact)
            else:
                choosen_opposite_pos = Order_Ag_sequence( choosen_opposite_pos,mutations, selected_contact)            

        # Load the amino acids to the paires according to the choosen_epitope_pos          
        mutations_aa = []
        for i in mutations:
            mutations_aa.append(pdbid_sequence[chain][i])
        opposite_aa = []
        for i in choosen_opposite_pos:
            opposite_aa.append(pdbid_sequence[opposite_chain][i]) 
        # Make a deep copy to be safe       
        kept_opposite_aa = copy.deepcopy(opposite_aa)
        kept_mutations_aa = copy.deepcopy(mutations_aa)
        # Here we have to make sure the Ab amino acids is the first element of 
        # the paires and the Ag amino acids is the second element of the paires.
        if Ab_Ag == 'Ab':
            paires.append([kept_mutations_aa, kept_opposite_aa])
        elif Ab_Ag == 'Ag':
            paires.append([kept_opposite_aa, kept_mutations_aa]) 
    
    # Load the results
    mutation_match_parameter['paires'] = paires
            
    return mutation_match_parameter
######################################################################
'''
Original_mutation_sets:
    to get the original matched parepi sets and the mutationed parepi sets
Input:
    mutation_match_parameter:
        with one more keys than the above
        
        ###mutation_match_parameter['mutations_aa']= [aa,aa]###
        
        a list of amino acids corresponding tho the mutations

Output:
    mutation_match_parameter with one more key value
    mutation_match_parameter['original_mutation_sets']:
        in the form of [[sets of original parepi],[sets of mutated parepi]]
'''
def Original_mutation_sets(mutation_match_parameter):
    # take out the values
    mutation_aa = mutation_match_parameter['mutations_aa']
    Ab_Ag = mutation_match_parameter['Ab_Ag']
    paires = mutation_match_parameter['paires']
    # Do the work
    mutated_paires = []
    if Ab_Ag == 'Ab':
        for parepi in paires:
            mutated_paires.append([mutation_aa, parepi[1]])
    elif Ab_Ag == 'Ag':
        for parepi in paires:
            mutated_paires.append([parepi[0], mutation_aa])
    # Load the results to original_mutation_sets.
    original_mutation_sets = [paires, mutated_paires]
    mutation_match_parameter['original_mutation_sets'] = original_mutation_sets
    
    return mutation_match_parameter
           
#######################################################################
'''
Prediction:
    a function to calculate the predicted values of the testing_set
    
Input:
    testing_set:
        a list of the given Ab_aa and Ag_aa paires in the form of 
       [[[ASP, LYS], [ALA, ARG, GLU]],......]
       
    wd_results:
        gives the working directory of the results
        
    wd_negative_samples:
        gives the working directory of the negative_samples
Output:
    predictions:
        an arrey in the shape of (len(testing_data), a), gives the predicted values
        of the paires.
'''
def Predition(testing_set, wd_results, wd_negative_samples):
#testing_set = testing5
    header = str(len(testing_set[0][0]))+'_'+str(len(testing_set[0][1]))
    #header
    os.chdir(wd_results)
    with open(header+'_aa_train', 'r') as f:
        positive_training_set = json.load(f)        
    with open(header+'_test_results_contact', 'r') as f:
        results_contact = json.load(f)   
    with open(header +'_test_results', 'r') as f:
        results_binary = json.load(f)
    #results_binary.keys()
    #len(results_binary['centers'])        
    os.chdir(wd_negative_samples)
    with open(header+'_train_negative', 'r') as f:
        negative_training_set = json.load(f)    
    
    results = results_binary # We can alter to choose the set of results we want to use.
    #results_contact.keys()

    centers = results['centers']
    # Calculate the testing design matrix
    training_set = copy.deepcopy(positive_training_set)
    training_set.extend(negative_training_set)
    selected_centers = []
    for i in centers:
        selected_centers.append(training_set[i])
    
    testing_design_matrix = np.zeros((len(testing_set), len(selected_centers)))
    testing_distance_matrix = np.zeros_like(testing_design_matrix)
    for i in range(len(testing_set)):
        for j in range(len(selected_centers)):
            Ab_seq1 = To_seq(testing_set[i][0])
            Ag_seq1 = To_seq(testing_set[i][1])
#            print(selected_centers[j])
            Ab_seq2 = To_seq(selected_centers[j][0])
            Ag_seq2 = To_seq(selected_centers[j][1])
            
            testing_distance_matrix[i,j] = Multiplication_distance(Ab_seq1, Ab_seq2, Ag_seq1, Ag_seq2)
    testing_design_matrix = np.exp(- testing_distance_matrix)
    testing_design_matrix = np.hstack((testing_design_matrix, np.ones((len(testing_set),1))))

    coeff = results['all_coeff'][0]['coeff']
    coeff = np.array(coeff).reshape((-1,1))
    #Calculate the prediction results
    predictions = testing_design_matrix.dot(coeff)
    
    return predictions
  
##################################################################
'''
Compare:
    a function to compare the predicted values of the original data set and
    the mutated data set
Input:
     wd_results: The same as above
     wd_negative_samples:The same as above
     mutation_match_parameter['fold_change']: a float number, give the fold
     of the Kd changes
Output:
    compare:
        a list of float numbers, [original_pred_sum, mutation_pred_sum]
'''
def Compare(wd_results, wd_negative_samples, mutation_match_parameter):
    # Take out the values
#    fold_change = mutation_match_parameter['fold_change']
    original_mutation_sets = mutation_match_parameter['original_mutation_sets']
    # Get the original sets and the mutated sets
    original_sets = original_mutation_sets[0]
    mutated_sets = original_mutation_sets[1]
    # Make predictions for each set
    original_pred = Predition(original_sets, wd_results, wd_negative_samples)
    mutated_pred = Predition(mutated_sets, wd_results, wd_negative_samples)
    
    return [np.sum(original_pred), np.sum(mutated_pred)]

######################################################################
'''
Preprocessing

Input:
    one_element_from_selected_affinity
    good_matched_ids
    good_combined_ids
    
Output:
    list_1
    list_2
    fold_change_1
    fold_change_2
'''
def Preprocessing(one_element_from_selected_affinity, good_matched_ids, good_combined_ids):
    list_1 = []
    list_2 = []
    # Take out the pdbids
    pdbid_1 = one_element_from_selected_affinity['affinities'][0][0]
    pdbid_2 = one_element_from_selected_affinity['affinities'][0][1]
    # Take out the mutations
    mutations = one_element_from_selected_affinity['mutations']

    # Load the data to the mutaion_match_parameter dictionary
    for mutation in mutations:
        if mutation[0] == pdbid_1:
            # Make sure there is no duplicates
            for opposite in mutation[-1]:
                list_1.append(Load_mutation_match_parameter(mutation,\
                                              opposite, pdbid_1, good_matched_ids,good_combined_ids))
        if mutation[0] == pdbid_2:
            for opposite in mutation[-1]:
                list_2.append(Load_mutation_match_parameter(mutation,\
                                              opposite, pdbid_2, good_matched_ids,good_combined_ids))
                
    # Calculate the fold_changes
    affinity_1 = float(one_element_from_selected_affinity['affinities'][1][0])
    affinity_2 = float(one_element_from_selected_affinity['affinities'][1][1])
    
    fold_change_1 = affinity_1/affinity_2
    fold_change_2 = affinity_2/affinity_1
    
    return list_1, list_2, fold_change_1, fold_change_2

'''
This is a subfunction of Preprocessing
'''
def Load_mutation_match_parameter(mutation, opposite, pdbid, good_matched_ids,good_combined_ids):                
        temp_mutation_match_parameter = {}
        # Load the pdbid
        temp_mutation_match_parameter['pdbid'] = pdbid
        # Load the mutation chain
        temp_mutation_match_parameter['mutation_chain'] = mutation[1]
        # Load the mutations
        temp_mutation_match_parameter['mutations'] = mutation[2]
        # Load the mutations_aa
        temp_mutation_match_parameter['mutations_aa'] = mutation[3]
        # Load the opposite chain 
        temp_mutation_match_parameter['opposite_chain'] = opposite
        # Load the matched_ids
        matched_ids = good_matched_ids[pdbid]
        temp_mutation_match_parameter['matched_ids'] = matched_ids
        # Load teh combined_ids
        combinded_ids = good_combined_ids[pdbid]
        temp_mutation_match_parameter['combined_ids'] = combinded_ids
        # Deep copy and add to list_1
        return temp_mutation_match_parameter      

                    


            
            
#pdbid = mutation_match_parameter['pdbid']
#chain = mutation_match_parameter['mutation_chain']
#mutations = mutation_match_parameter['mutations']
#matched_ids = mutation_match_parameter['matched_ids']
#combined_ids = mutation_match_parameter['combined_ids']
#opposite_chain = mutation_match_parameter['opposite_chain']
#Ab_Ag = mutation_match_parameter['Ab_Ag']        
#mutation_match_parameter['mutations_aa']    

        
        
        
######################################################################
os.chdir('/home/leo/Documents/Database/Pipeline/Complexes with Affinity')
with open('good_matched_ids', 'r') as f:
    good_matched_ids = json.load(f)
with open('good_combined_ids', 'r') as f:
    good_combined_ids = json.load(f)
with open('contact', 'r') as f:
    contact = json.load(f)    
with open('sequence', 'r') as f:
    sequence = json.load(f)

os.chdir('/home/leo/Documents/Database/Pipeline/All with peptide 5+ resolution 4A/structure')




#contact['1yy9']
#sequence['1mlc']['F'][65:67]
#print('executed')
##################################################################


#####################################################################

'''
Get the testing design matrix
'''

# Import the trained results
directory_paires = [['/home/leo/Documents/Database/Pipeline/Results/1_free',\
  '/home/leo/Documents/Database/Pipeline/Negative samples/1_free'],\
 ['/home/leo/Documents/Database/Pipeline/Results/0_free',\
                       '/home/leo/Documents/Database/Pipeline/Negative samples/0_free']]
        
wd_results = directory_paires[0][0]
wd_negative_samples = directory_paires[0][1]

#################################################################################### 

pdbid = '1'
#contact[pdbid]
matched_ids = good_matched_ids[pdbid]
matched_ids
combined_ids =good_combined_ids[pdbid]
combined_ids

chain = 'C'
opposite_chain='B'
#sequence[pdbid][chain][4]
#######################################
mutation_match_parameter = {}
#Need to be filled
mutation_match_parameter['mutations']=[96] 
mutation_match_parameter['opposite_chain']=opposite_chain 
mutation_match_parameter['mutations_aa']=['ALA' ] 
mutation_match_parameter['fold_change'] = 85
Ab_Ag = mutation_match_parameter['Ab_Ag'] = 'Ag'
# Load the mutation parameter

mutation_match_parameter['pdbid']=pdbid
mutation_match_parameter['mutation_chain']=chain 
mutation_match_parameter['matched_ids']=matched_ids 
mutation_match_parameter['combined_ids']=combined_ids




os.chdir('/home/leo/Documents/Database/Pipeline/All with peptide 5+ resolution 4A/structure')    
mutation_match_parameter = Paire_select(mutation_match_parameter)
mutation_match_parameter['paires']

mutation_match_parameter = Original_mutation_sets(mutation_match_parameter)
mutation_match_parameter['original_mutation_sets']

original_mutated_pred = Compare(wd_results, wd_negative_samples, mutation_match_parameter)
original_mutated_pred = np.array(original_mutated_pred)
original_mutated_pred





###################################################################################  
def Combination(sequence):
    combinations = []
    for i in range(len(sequence)-1):
        for j in range(i+1, len(sequence)):
            combinations.append([i,j])
    return combinations


total = 0
correct = 0
for i in range(21):
    print(i)
    testing = testings[i]
    order = orders[i]
    combination = Combination(order)
    predictions = Predition(testing, wd_results,  wd_negative_samples)  
    # Permute the order
    for comb in combination:
        product = (predictions[comb[0]]-predictions[comb[1]])*(order[comb[0]]-order[comb[1]])
        if product > 0:
            correct += 1
            total +=1
        elif product < 0:
            total += 1
            
correct
total
#predictions = Predition(testing2, wd_results,  wd_negative_samples)   
   
#############################################################
''''''
''''''
''''''    
''''''
''''''
''''''
''''''
''''''
''''''
'''
Check the code
'''
''''''
''''''
''''''
for i in range(len(testings)):
    if len(testings[i]) == len(orders[i]):
        print('yes')
        
testing_set = testings[0]    
header = str(len(testing_set[0][0]))+'_'+str(len(testing_set[0][1]))
#header

#header
os.chdir(wd_results)
with open(header+'_aa_train', 'r') as f:
    positive_training_set = json.load(f)
    
with open(header+'_test_results_contact', 'r') as f:
    results_contact = json.load(f)
#results_contact.keys()
#len(results_contact['centers'])

with open(header +'_test_results', 'r') as f:
    results_binary = json.load(f)
#results_binary.keys()
#len(results_binary['centers'])

    
os.chdir(wd_negative_samples)
with open(header+'_train_negative', 'r') as f:
    negative_training_set = json.load(f)

    
#results_contact.keys()
#    results = results_binary # Can be adjusted to see which one is better, binary or contact
#len(results_contact['centers'])
centers = results_contact['centers']
# Calculate the testing design matrix
training_set = copy.deepcopy(positive_training_set)
training_set.extend(negative_training_set)
selected_centers = []
for i in centers:
    selected_centers.append(training_set[i])
#selected_centers = np.array(training_set)[centers].flatten().tolist()
#selected_centers[-6:]

testing_design_matrix = np.zeros((len(testing_set), len(selected_centers)))
testing_distance_matrix = np.zeros_like(testing_design_matrix)
for i in range(len(testing_set)):
    for j in range(len(selected_centers)):
        Ab_seq1 = To_seq(testing_set[i][0])
        Ag_seq1 = To_seq(testing_set[i][1])
#        print(selected_centers[j])
        Ab_seq2 = To_seq(selected_centers[j][0])
        Ag_seq2 = To_seq(selected_centers[j][1])
        
        testing_distance_matrix[i,j] = Multiplication_distance(Ab_seq1, Ab_seq2, Ag_seq1, Ag_seq2)
testing_design_matrix = np.exp(- testing_distance_matrix)
testing_design_matrix = np.hstack((testing_design_matrix, np.ones((len(testing_set),1))))
testing_design_matrix.shape
len(centers)
#testing_design_matrix.shape

# Take out the coefficient from the testresults
#results_contact.keys()
#type(results_contact['all_coeff'])
#len(results_contact['all_coeff'][0])
#results_contact['all_coeff'][0].keys()
coeff = results_contact['all_coeff'][0]['coeff']
coeff = np.array(coeff).reshape((-1,1))
#Calculate the prediction results
predictions = testing_design_matrix.dot(coeff)
    
    