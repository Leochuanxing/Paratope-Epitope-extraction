import copy
import numpy as np
import json
import os
import math
os.chdir('/home/leo/Documents/Database/Pipeline/Codes and testing data')
from AAC_2 import Coordinates
from AAC_2 import Get_contact, Chain_seq
from FrameConstraint import Add, Sub_seq, Get_consecutive

from RBFN_coverage import To_seq
from RBFN_coverage import Multiplication_distance

from Logistic_Regression import Similarity_matrix, Truncated



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
def Select_contact_opposite(mutation_match_parameter, sequence):
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
#    with open(pdbid+'.pdb', 'r') as file:
#        pdbid_sequence = Chain_seq(file, combined_ids)
    pdbid_sequence = sequence[pdbid]
    # store the paires in paires
    cutoff = 4
    possible_opposite = []
    if Ab_Ag == 'Ab':
        chain_pos = 2
        opposite_chain_pos = 3
        aa_pos = 1
        opposite_aa_pos = 2
    elif Ab_Ag == 'Ag':
        chain_pos = 3
        opposite_chain_pos = 2
        aa_pos = 2
        opposite_aa_pos = 1
        
    while possible_opposite == [] and cutoff <=6:
        contact = Get_contact(cdn, matched_ids, cutoff)
        # Carry out the above process:
        # take out all the contact containing the chain_name
        selected_contact = []
        possible_opposite = []
        for i in contact:
            if chain == i[0][chain_pos] and i[aa_pos] in mutations:                
                if i[0][opposite_chain_pos] == opposite_chain:
                    selected_contact.append(i)
                    possible_opposite.append(i[opposite_aa_pos])                  
                
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
def Paire_select(mutation_match_parameter, sequence):

    # Extract the required information from the pdb file
    selected_contact, possible_opposite, pdbid_sequence = Select_contact_opposite(mutation_match_parameter, sequence)
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
        for choosen_opposite_pos in longest_possible_consecutives:
            # Correct the direction of the choosen_opposite_pos
            choosen_opposite_pos.sort()
            # define a small function to change the order of the paires
    #        if len(mutations) >= 2 and len(choosen_opposite_pos)>=2:
            if len(choosen_opposite_pos)>=2:
                if Ab_Ag == 'Ab':
                    choosen_opposite_pos = Order_Ag_sequence(mutations, choosen_opposite_pos, selected_contact)
                elif  Ab_Ag == 'Ag':
                    mutations = Order_Ag_sequence(choosen_opposite_pos,mutations, selected_contact)            
    
            # Load the amino acids to the paires according to the choosen_epitope_pos          
            original_aa = []
            for i in mutations:
                original_aa.append(pdbid_sequence[chain][i])
            opposite_aa = []
            for i in choosen_opposite_pos:
                opposite_aa.append(pdbid_sequence[opposite_chain][i]) 
            # Make a deep copy to be safe       
            kept_opposite_aa = copy.deepcopy(opposite_aa)
            kept_original_aa = copy.deepcopy(original_aa)
            # Here we have to make sure the Ab amino acids is the first element of 
            # the paires and the Ag amino acids is the second element of the paires.
            if Ab_Ag == 'Ab':
                paires.append([kept_original_aa, kept_opposite_aa])
            elif Ab_Ag == 'Ag':
                paires.append([kept_opposite_aa, kept_original_aa]) 
    
    # Load the results
    mutation_match_parameter['paires'] = paires
            
#    return mutation_match_parameter
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
    if paires != []:
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
    else:
        mutation_match_parameter['original_mutation_sets'] = []
    
#    return mutation_match_parameter
           
#######################################################################
'''
Prediction_RBFN_coverage:
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
def Predition_RBFN_coverage(testing_set, wd_results, wd_negative_samples):
#testing_set = testing5
    header = str(len(testing_set[0][0]))+'_'+str(len(testing_set[0][1]))
    #header
    os.chdir(wd_results)
    with open(header+'_aa_train', 'r') as f:
        positive_training_set = json.load(f)        
#    with open(header+'_test_results_contact', 'r') as f:
#        results_contact = json.load(f)   
    with open(header +'_test_results', 'r') as f:
        results = json.load(f)
    #results_binary.keys()
    #len(results_binary['centers'])        
    os.chdir(wd_negative_samples)
    with open(header+'_train_negative', 'r') as f:
        negative_training_set = json.load(f)    
    

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
Prediction_LogisticRegression:
    this function is to make prediction by using the Logistic regression model
Input: 
    testing_set, wd_results, wd_negative_samples
    They are the same as in the Prediction_RBFN_coverage
Output:
    predictions:
        the same as in the Prediction_RBFN_coverage.
'''
def Prediction_LogisticRegression(testing_set, wd_results, wd_negative_samples):
    # Create a header
    header = str(len(testing_set[0][0]))+'_'+str(len(testing_set[0][1]))
    #header
    os.chdir(wd_results)
    with open(header+'_aa_train', 'r') as f:
        positive_training_set = json.load(f)         
    with open(header +'_results_LogisticRegression', 'r') as f:
        results_Logistic = json.load(f)
    #results_binary.keys()
    #len(results_binary['centers'])        
    os.chdir(wd_negative_samples)
    with open(header+'_train_negative', 'r') as f:
        negative_training_set = json.load(f) 
        
    # Take out the values from the results_Logistic
    best_percentage = results_Logistic['best_percentage']
    coeff = results_Logistic['coeff']
    
    # Calculate the similarity matrix
    training_set = copy.deepcopy(positive_training_set)
    training_set.extend(negative_training_set)
    similarity_matrix = Similarity_matrix(testing_set, training_set, square = False)
    # Truncate the Logistic model
    p = len(positive_training_set)
    n = len(negative_training_set)
    number = math.floor((n+p)*best_percentage)
    truncated_sum  = Truncated(similarity_matrix, p, n, number)
    # Make prediction
    truncated_sum_augement = np.hstack((truncated_sum, np.ones((len(testing_set),1))))
    predictions = truncated_sum_augement.dot(coeff)
    
    # We can also calculate the logits and return the logits as well
    logits = np.exp(predictions)
    
    return logits
    
#########################################################################
'''
Compare:
    a function to compare the predicted values of the original data set and
    the mutated data set
Input:
     
     wd_results: The same as above
     wd_negative_samples:The same as above
     mutation_match_parameter['fold_change']: a float number, give the fold
     of the Kd changes
     model: a string, it takes the value of either 'RBFN' or 'Logistic'
     if it takes the value of 'RBFN', we use the RBFN model to make prediction
     if it takes the value of 'Logistic', we use the Logistic regression model to make prediction.
Output:
    compare:
        a list of float numbers, [original_pred_sum, mutation_pred_sum]
'''

def Compare(wd_results, wd_negative_samples, mutation_match_parameter, model = 'RBFN'):
    # Take out the values
#    fold_change = mutation_match_parameter['fold_change']
    original_mutation_sets = mutation_match_parameter['original_mutation_sets']
    # Get the original sets and the mutated sets
    if original_mutation_sets != []:
        original_sets = original_mutation_sets[0]
        mutated_sets = original_mutation_sets[1]
        # Make predictions for each set
        if model == 'RBFN':
            original_pred = Predition_RBFN_coverage(original_sets, wd_results, wd_negative_samples)
            mutated_pred = Predition_RBFN_coverage(mutated_sets, wd_results, wd_negative_samples)
        if model == 'Logistic':
            original_pred = Prediction_LogisticRegression(original_sets, wd_results, wd_negative_samples)
            mutated_pred = Prediction_LogisticRegression(mutated_sets, wd_results, wd_negative_samples)
        original_mutation_score = [np.sum(original_pred), np.sum(mutated_pred)]
    else:
        return 'Empty Match'        
        
    return original_mutation_score

######################################################################

####################################################################

'''
Preprocessing

Input:
    one_element_from_selected_affinity
    good_matched_ids
    good_combined_ids
    Ab_Ag
Output:
    list_1: 
        a list of mutation_match_parameter dictionaries 
    list_2:
        similary as list_1
    # A list is a basic unit to compare mutations with the originals. A list 
    # contains all the mutations combined together to form a complex, whose affinity
    # is measured.
    
'''
def Preprocessing(one_element_from_selected_affinity, good_matched_ids, good_combined_ids, Ab_Ag):
    list_1 = []
    list_2 = []
    # Take out the pdbids
    pdbid_1 = one_element_from_selected_affinity['affinities'][0][0]
    pdbid_2 = one_element_from_selected_affinity['affinities'][0][1]
    # Calculate the fold_changes
    affinity_1 = float(one_element_from_selected_affinity['affinities'][1][0])
    affinity_2 = float(one_element_from_selected_affinity['affinities'][1][1])
    
    fold_change_1 =[affinity_1, affinity_2, affinity_1/affinity_2]
    fold_change_2 =[affinity_2, affinity_1, affinity_2/affinity_1]
    
    # Take out the mutations
    mutations = one_element_from_selected_affinity['mutations']

    # Load the data to the mutaion_match_parameter dictionary
    for mutation in mutations:
        if mutation[0] == pdbid_1:
            for opposite in mutation[-1]:
                list_1.append(Load_mutation_match_parameter(mutation,\
                                              opposite, pdbid_1, good_matched_ids,\
                                              good_combined_ids, fold_change_1,pdbid_2, Ab_Ag))
        if mutation[0] == pdbid_2:
            for opposite in mutation[-1]:
                list_2.append(Load_mutation_match_parameter(mutation,\
                                              opposite, pdbid_2, good_matched_ids,\
                                              good_combined_ids, fold_change_2, pdbid_1,Ab_Ag))    
    
    return list_1, list_2
'''
This is a subfunction of Preprocessing
'''
def Load_mutation_match_parameter(mutation, opposite, pdbid, good_matched_ids,\
                                  good_combined_ids,fold_change, mutate_to, Ab_Ag):                
        temp_mutation_match_parameter = {}
        # Load the Ab_Ag
        temp_mutation_match_parameter['Ab_Ag'] = Ab_Ag
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
        # Load the fold_change
        temp_mutation_match_parameter['fold_change'] = fold_change
        # Load the mutate_to
        temp_mutation_match_parameter['mutate_to'] = mutate_to
        ''' mutate_to gives which complex it mutates to, if we don't know the 
             pdbid of the complex it mutates to, mutate_to = 'No'. '''
        return temp_mutation_match_parameter      

####################################################################################                    
'''
Mutation_list_prediction:
    to make prediction for a mutation list
input:
    one_list: 
        a list containing a set of mutations corresponding to measured affinity
Output:
    scores:
        either[original_score, mutation_score]
        or ['Empty Match']
'''

def Mutation_list_prediction(one_list, wd_results, wd_negative_samples, sequence,  model = 'RBFN'):  

    # Select the the paires and Generate the original and mutation sets
#    list_parameters = []
    for mutation_match_parameter in one_list:
#        mutation_match_parameter = Paire_select(mutation_match_parameter,sequence)
#        list_parameters.append(copy.deepcopy(Original_mutation_sets(mutation_match_parameter)))
        Paire_select(mutation_match_parameter,sequence)
        Original_mutation_sets(mutation_match_parameter)
    # Make prediction
    # When we make prediction, we simply assume each paratope-epitope paire contribute equally to 
    # the final affinity
    all_original_mutation_score= []
    for parameter in one_list:
        original_mutation_score = Compare(wd_results, wd_negative_samples, parameter, model)
        # Everytime after Compare we have to switch the working directory back. THis is not good.
        os.chdir('/home/leo/Documents/Database/Pipeline/Complexes with Affinity/structure')
        # Load the results to scores
        all_original_mutation_score.append(original_mutation_score)

    # Check if there are empty match, if so this prediction doesn't count 
    for j in all_original_mutation_score:
        if j == 'Empty Match':
            # If there is one empty match, the whole prediction ends and return the value 'Empty Match'
            return [j]

    # Add up the scores
    scores_original=0
    scores_mutation = 0
    for original_mutaion in all_original_mutation_score:
        scores_original += original_mutaion[0]
        scores_mutation += original_mutaion[1]
    
    return [scores_original, scores_mutation]

######################################################################
#os.chdir('/home/leo/Documents/Database/Pipeline/Complexes with Affinity')
#with open('good_matched_ids', 'r') as f:
#    good_matched_ids = json.load(f)
#with open('good_combined_ids', 'r') as f:
#    good_combined_ids = json.load(f)
#with open('contact', 'r') as f:
#    contact = json.load(f)    
#with open('sequence', 'r') as f:
#    sequence = json.load(f)

#os.chdir('/home/leo/Documents/Database/Pipeline/Complexes with Affinity/structure')


##################################################################
'''
List_list:
    a function to create a list of lists
Input: 
    workable_dicts_Ab_fixed:
        as explained in Select_affinity
    workable_dicts_Ag_fixed:
        as explained in Select_affinity
Output:
    list_list: 
        a list, each element in this list is a list of mutation_match_parameter dictionaries.
        The elements are the basic unit to make predictions, as explained in the Output of Preprocessing
'''
def List_list(workable_dicts_Ab_fixed, workable_dicts_Ag_fixed):
    # Create an empty container
    list_list = []
    # Deal with workable_dicts_Ab_fixed
    for i in range(len(workable_dicts_Ag_fixed)):
        one_element_from_selected_affinity = workable_dicts_Ag_fixed[i]
        list_1, list_2 = Preprocessing(one_element_from_selected_affinity,\
                                  good_matched_ids, good_combined_ids, Ab_Ag = 'Ab')
        
        list_list.append(copy.deepcopy(list_1))
        list_list.append(copy.deepcopy(list_2))
        
    #Deal with workable_dicts_Ab_fixed
    for i in range(len(workable_dicts_Ab_fixed)):
        one_element_from_selected_affinity = workable_dicts_Ab_fixed[i]
        list_1, list_2 = Preprocessing(one_element_from_selected_affinity,\
                                  good_matched_ids, good_combined_ids, Ab_Ag = 'Ag')
        
        list_list.append(copy.deepcopy(list_1))
        list_list.append(copy.deepcopy(list_2))
    return list_list
##################################################################
'''
Initiate_mutation_match_parameter:
    a function to initiate the mutation_match_parameter for later usage
Input:
    good_matched_ids:
        The same as we always use.
    good_combined_ids:
        As above
    one_mutation:
        a dictionary, with keys 'mutations', 'affinities'. This dictionay is one of
        the element of the recorded mutaions from the literatures
Output:
    unit_list:
        a list of mutation_match_parameters, and this list is a basic prediction unit
'''
def Initiate_mutation_match_parameter(good_matched_ids, good_combined_ids, one_mutation):
    # Take out the values from the one_mutations
    mutations_info = one_mutation['mutations']
    affinities = one_mutation['affinities']
    mutate_to = affinities[0][1]
    # Calculate the fold change
    fold = float(affinities[1][0])/float(affinities[1][1])
    #Load up the fold_change
    fold_change = [affinities[1][0], affinities[1][1], affinities[2], fold]
    # Finde the pdbid
    pdbid = mutations_info[0][0]
    # Find the combined ids
    combined_ids = good_combined_ids[pdbid]
    # Find the matched ids
    matched_ids = good_matched_ids[pdbid][0]
    # create an empty list
    unit_list = []
    # Load up the list
    for i in range(len(mutations_info)):
        sub_mutaion = mutations_info[i]
        # Find the mutaion chain
        mutation_chain = sub_mutaion[1]    
        # Find the value of Ab_Ag
        for i in range(3):
            if matched_ids[i] == mutation_chain and i == 2:
                Ab_Ag = 'Ag'
                opposite_chains = [matched_ids[0], matched_ids[1]]
            else:
                Ab_Ag = 'Ab'
                opposite_chains = [matched_ids[2]]
        # Find the opposite chains
#        opposite_chains = sub_mutaion[4]
        for opposite in opposite_chains:
            # Creat an empty dictionary
            temp_mutation_match_parameter = {}
            # Load the Ab_Ag
            temp_mutation_match_parameter['Ab_Ag'] = Ab_Ag
            # Load the pdbid
            temp_mutation_match_parameter['pdbid'] = pdbid
            # Load the mutation chain
            temp_mutation_match_parameter['mutation_chain'] = mutation_chain
            # Load the mutations
            temp_mutation_match_parameter['mutations'] = sub_mutaion[2]
            # Load the mutations_aa
            temp_mutation_match_parameter['mutations_aa'] = sub_mutaion[3]
            # Load the opposite chain 
            temp_mutation_match_parameter['opposite_chain'] = opposite
            # Load the mutate to
            temp_mutation_match_parameter['mutate_to'] = mutate_to
            # Load the fold change
            temp_mutation_match_parameter['fold_change'] = fold_change
            # Load the matched_ids
            temp_mutation_match_parameter['matched_ids'] = [matched_ids]
            # Load the combined_ids
            temp_mutation_match_parameter['combined_ids'] = combined_ids            
            
            # Add the dictionary to the unit_list
            unit_list.append(copy.deepcopy(temp_mutation_match_parameter))
    
    return unit_list
############################################################
'''
List_list_other:
    a function to generate similar output as List_list, the only difference is that 
    the input is a little bit different from List_list
Input:
    other_mutations: mutations collected from papers
    good_matched_ids:
    good_combined_ids:
    
Output:
    list_list_other:
        in the same form as the output of List_list
'''
def List_list_other(other_mutations, good_matched_ids, good_combined_ids):
    list_list_other = []
    # Initiate the mutation_match_parameters
    for one_mutation in other_mutations:
        unit_list = Initiate_mutation_match_parameter(good_matched_ids, good_combined_ids, one_mutation)
#        for mutation_match_parameter in unit_list:
#            Paire_select(mutation_match_parameter)
#            Original_mutation_sets(mutation_match_parameter)
#            # Load to the list_list_other
        list_list_other.append(copy.deepcopy(unit_list))
        
    return list_list_other
    
#################################################################


################################################################################
'''
Predict_list_list:
    a function to make prediction on the given list_list
Input:
    list_list:
        a list of  basic prediction unit, which is list of mutation_match_parameter dictionaries
Output:
    results_list: a list, each element is a 
        a dictionary in the following form:
        results['original'] = original pdbid
        results['mutate_to'] = the pidbid of the complex the original complex is mutated to
        results['fold_change'] = the fold change of the affinity  original_affinity/mutation_affinity
        results['predictions'] = [original_prediction_value, mutation_prediction_value]
'''
def Predict_list_list(list_list, wd_results, wd_negative_samples,sequence, model):
    results_list = []
    for i in range(len(list_list)):
        list_1 = list_list[i]
        # Change the working directory
        os.chdir('/home/leo/Documents/Database/Pipeline/All with peptide 5+ resolution 4A/structure')
               
        # Create a dictionary to contain the detailed prediction results
        results_1 = {}
        # Make prediction for the lists
        if list_1 != []:
            prediction_1 = Mutation_list_prediction(list_1, wd_results, wd_negative_samples,sequence, model)
            # Store the results
            results_1['original'] = list_1[0]['pdbid']
            results_1['mutate_to'] = list_1[0]['mutate_to']
            results_1['fold_change'] = list_1[0]['fold_change']
            results_1['predictions'] = prediction_1
       
        if results_1 != {}:
            results_list.append(results_1)
    return results_list

#####################################################################
'''
Right_or_wrong:
    a function to tell whether our  prediction is right or wrong
Input:
    fold: a float, shows how many times the affinity changes
    affinity_type: string, either 'Ka' or 'Kd'
    prediction1: a float, give the score of the first prediction
    predictin2: a float, give the score of the second prediction
Output:
    right_or_wrong:
        a boolean, if the prediction is correct, it takes the values of True
                   if the prediction is incorrect, it takes the values of False
'''
def Right_or_wrong(fold, affinity_type, prediction1, prediction2):
    # In the case when the affinity_type is Kd
    if affinity_type == 'Kd':
        if fold > 1 and prediction1 < prediction2:
            right_or_wrong = True
        elif fold < 1 and prediction1 > prediction2:
            right_or_wrong = True
        else:
            right_or_wrong = False
            
    # In the case when the affinity_type is Ka      
    elif affinity_type == 'Ka':
        if fold > 1 and prediction1 > prediction2:
            right_or_wrong = True
        elif fold < 1 and prediction1 < prediction2:
            right_or_wrong = True
        else:
            right_or_wrong = False
    
    return right_or_wrong
########################################################################
'''
Count the number of predictions and the number of correct predictions
'''
def Count_correct(results_list, fold_change_cut = 1):
    total_prediction = 0
    correct_prediction = 0
    
    for results in results_list:
        # Take out the values
        fold_change = results['fold_change']
        fold = fold_change[3]
        affinity_type = fold_change[2]
        predictions = results['predictions']
        
        if predictions[0] != 'Empty Match' and fold >= fold_change_cut: 
            total_prediction += 1
            if Right_or_wrong(fold, affinity_type, predictions[0], predictions[1]):
                correct_prediction += 1
                
        elif predictions[0] != 'Empty Match' and fold <= 1/fold_change_cut:
            total_prediction += 1
            if Right_or_wrong(fold, affinity_type, predictions[0], predictions[1]):
                correct_prediction += 1

    return total_prediction, correct_prediction
###############################################################################
'''
Do the prediction, let us use one paires of complexes as an example to build up the 
steps
'''
#  Assign the working directories
directory_paires = [['/home/leo/Documents/Database/Pipeline/Results/1_free',\
  '/home/leo/Documents/Database/Pipeline/Negative samples/1_free'],\
 ['/home/leo/Documents/Database/Pipeline/Results/0_free',\
                       '/home/leo/Documents/Database/Pipeline/Negative samples/0_free']]
        
wd_results = directory_paires[0][0]
wd_negative_samples = directory_paires[0][1]

# Open the workable list
os.chdir('/home/leo/Documents/Database/Pipeline/Affinity/All_structures')
with open('mutations', 'r') as f:
    mutations = json.load(f)
with open('combined_ids', 'r') as f:
    good_combined_ids = json.load(f)
with open('matched_ids', 'r') as f:
    good_matched_ids = json.load(f)
with open('mutations_mix', 'r') as f:
    mutations_mix = json.load(f)
with open('sequence', 'r') as f:
    sequence = json.load(f)
#len(mutations)
############################################################
if __name__ == '__main__':
    list_list = List_list_other(mutations, good_matched_ids, good_combined_ids)
    results_list = Predict_list_list(list_list, wd_results, wd_negative_samples,sequence, model = 'RBFN')
    total_prediction, correct_prediction = Count_correct(results_list, fold_change_cut = 1)   
    total_prediction
    correct_prediction
    results_list

correct_prediction
total_prediction
results_list
len(results_list)
n_pred = 0
for results in results_list:
    if results['predictions'][0] != 'Empty Match':
        n_pred += 1

n_pred
os.chdir('/home/leo/Documents/Database/Pipeline/Affinity/All_structures')
with open('results_list', 'w') as f:
    json.dump(results_list, f)
