###############################################################
# Import the modules
#import random
import numpy as np
import os
import json
import math
import copy
from matplotlib import pyplot as plt


###########################################################
os.chdir('/home/leo/Documents/Database/Pipeline/Codes and testing data')
from RBFN_coverage import Design_matrix, Distance_matrix, Generate_testing_design_matrix,\
                            Train_RBFN_BFGS, Multiplication_distance, Loss
from RBFN_coverage import Mrakov, To_seq
#############################################################
# Define  distances
from Bio import Align
from Bio.SubsMat.MatrixInfo import blosum62
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = blosum62
aligner.open_gap_score = -5
aligner.extend_gap_score = -1
aligner.mode = 'global'
  

################################################################################
'''
parameter contains
    parameter['centers']
    parameter['design_matrix']
    parameter['observed_values']
    parameter['coeff']
    parameter['reg']
    parameter['list_n_centers']
    parameter['list_centers']
    parameter['list_coeff']
    parameter['list_reg']
    parameter['training_set']
    parameter['non_redundent_training_set']
'''
#######################################################################

'''
One_step_reduce_centers:
    a function to reduce the number of centers by 1
Input:
    parameter:
        a dictionary, contains all the parameters related to the model
        
        parameter['centers']:
            a list of centers, this list will be reduced by 1 element each time
        parameter['design_matrix']:
            a design matrix, with the columns corresponding to the centers
                    
      ##### Pay attention that the column number is one larger than the row number.#####
            
        parameter['observed_values']:
            an array in the shape of (np.shape(design_matrix), ), with length equal to the number of the rows of the design_matrix
        parameter['coeff']:
            An array in the shape of (len(centers)+1, 1), gives the coefficient of 
            corresponding to the centers and the constant term
        parameter['reg']:
            an array in the shape of (len(centers)+1, 1), gives the coefficient of the
            regularization term.

Output:
    parameter, a dictionary with all the values updated
'''
'''
We can add a control to speed up the process
'''
def One_step_reduce_centers(parameter):
    # Take out the values
    centers = parameter['centers']
    design_matrix = parameter['design_matrix']
    observed_values = parameter['observed_values']

#    ratio = 3000/len(centers)
    if len(centers) > 1500 :
        termination = 1.5*len(centers)
    elif len(centers) <= 1500 and len(centers) >1000:
        termination =  len(centers)
    elif len(centers)<= 1000 and len(centers)>500: 
        termination = 10
    elif len(centers) <= 500:
        termination = len(centers)/1000
#    termination =10* len(centers)

        
    parameter, loss = Train_RBFN_BFGS(design_matrix, observed_values, rho=0.85, c = 1e-3, termination=termination,\
                                      parameter_inheritance = True, parameter=parameter)
    coeff = parameter['coeff']
    # match the coeff up with the index of the centers
    # Take the absolute value
    coeff_list = np.abs(coeff)
    # match up the abs of coeff with other indxes
    match_up = []
    for i in range(len(centers)):
        match_up.append([centers[i], coeff_list[i], i])
    # sort the match up according to the value of coeff
    match_up.sort(key = lambda x:x[1])
    # To speed up the process, we remove more than one centers if it is the number 
    # of centers is larger
    removed_sets = []
    if len(match_up) > 2000:
        n_remove = math.floor(len(centers)*0.004)
    elif len(match_up) > 1500 and len(match_up) <= 2000:
        n_remove = math.floor(len(centers)*0.002)
    else:
        n_remove = 1
        
    for i in range(n_remove):
        removed_sets.append(match_up[i])

    
    # Update the centers
    removed_col = []
    for to_be_removed in removed_sets:    
        centers.remove(to_be_removed[0])
        # load the removed_col
        removed_col.append(to_be_removed[2])
        # Update the coefficient
    coeff = np.delete(coeff, removed_col)
    coeff = np.reshape(coeff, (-1, 1))
    # Update the design matrix
    design_matrix = np.delete(design_matrix, removed_col, 1)
    # Load the updated values to the parameter
    parameter['centers']= centers
    parameter['coeff'] = coeff
    parameter['reg'] = np.ones((len(coeff), 1))
    parameter['design_matrix'] = design_matrix
    return parameter

###############################################################################

'''
Remove the duplicates before removing the centers
'''
def Remove_duplicates(training_set):
    centers = []
    pure_parepi = []
    non_redundent_training_set = []
    for i in range(len(training_set)):
        parepi = training_set[i]
        if [parepi[0], parepi[1]] not in pure_parepi:            
            pure_parepi.append([parepi[0], parepi[1]])
            centers.append(i)
            non_redundent_training_set.append(training_set[i])
    return centers, non_redundent_training_set

#centers = Remove_duplicates(training_set)

#############################################################

'''
Generate_coeff_centers:
    a function to generate the coeff centers with the length of the centers equal to a given length
Input:
    parameter:
        a dictionary with one more value
        parameter['list_n_centers']: a list of integers, gives the number of centers to be kept
#    control_coeff:
#        a float, to speed up the selection of centers, when removing the centers, we remove the ones with 
#        the absolute values fo the coefficients less than control_coeff, otherwise, we remove the one with 
#        the least absolute value of the coefficient.
Output:
    parameter:
        with one more value
        parameter['list_centers']: a list with each element a list of centers
        parameter['list_coeff']: a list of coeff, each coeff is corresponding to a list of centers
        parameter['list_reg']: similar as above
'''

def Coeff_select_centers(parameter):
    # Take out the values from the parameter
    list_n_centers = parameter['list_n_centers']
    centers = parameter['centers']
    n_centers = len(centers)
    # Create empty lists
    list_centers = []
    list_coeff = []
    list_reg = []
    # Reduce the centers according the list_n_centers
    for n in list_n_centers:
        while n < n_centers:
            parameter = One_step_reduce_centers(parameter)
            # update the n_centers
            centers = parameter['centers']
            n_centers = len(centers)
            print('Length of centers: ', n_centers)
        list_centers.append(copy.deepcopy(centers))
        # Load list_coeff and list_reg
        coeff = parameter['coeff']
        reg = parameter['reg']
        list_coeff.append(copy.deepcopy(coeff))
        list_reg.append(copy.deepcopy(reg))
    
    # Load list_centers, list_coeff and list_reg to the parameter
    parameter['list_centers'] = list_centers
    parameter['list_coeff'] = list_coeff
    parameter['list_reg'] = list_reg
    
    return parameter
##################################################################################################
def Testing_design_matrix(testing_set, training_set):
    test_design_matrix = np.zeros((len(testing_set), len(training_set)))
    for i in range(len(testing_set)):
        test = testing_set[i]
        for j in range(len(training_set)):
            train = training_set[j]
            Ab_seq1 = To_seq(test[0])
            Ag_seq1 = To_seq(test[1])
            Ab_seq2 = To_seq(train[0])
            Ag_seq2 = To_seq(train[1])
            distance = Multiplication_distance(Ab_seq1, Ab_seq2, Ag_seq1, Ag_seq2)
            test_design_matrix[i, j] = Mrakov(distance, radius = 1)
            
    test_design_matrix = np.hstack((test_design_matrix, np.ones((len(testing_set), 1))))
    
    return test_design_matrix
            
#####################################################################################
'''
Area_calculation:
    a function to caltulate the area under the precision_recall curve
Input:
    parameter: a dictionary as given abobe
    testing_set: a list of testing samples
    observed_values: a list of observed values corresponding to testing_set
Output:
    parameter with one more value
    parameter['list_recall_precision'], a list of areas   
'''
def Recall_precision_calculation(parameter, testing_set, observed_values):
    list_coeff = parameter['list_coeff']
    list_centers = parameter['list_centers']
    training_set = parameter['training_set']
    list_recall_precision = []
#    list_testing_design_matrix = []
    # Calculate the relative position of the centers
    big_centers = list_centers[0]
    list_relative_pos = []
    for i in range(len(list_centers)):
        centers_pos = []           
        for j in range(len(list_centers[i])):
            for k in range(len(big_centers)):
                if big_centers[k] == list_centers[i][j]:
                    centers_pos.append(k)
        list_relative_pos.append(copy.deepcopy(centers_pos))
    # Calculate the big_testing_design_matrix
    big_training_set = []
    for i in big_centers:
        big_training_set.append(training_set [i])
    big_test_design_matrix = Testing_design_matrix(testing_set, big_training_set)
    # Go into one set of centers
    for i in range(len(list_centers)):
        centers_pos = list_relative_pos[i]
        coeff = list_coeff[i]
        centers_pos.append(-1)
        test_design_matrix = big_test_design_matrix[:, centers_pos]
        
        pred = test_design_matrix.dot(coeff)
        # Match up the observed values with the pred
        match_up = []
        for i in range(len(observed_values)):
            match_up.append([pred[i], observed_values[i]])
        match_up.sort(key = lambda x:x[0], reverse = True)
            
        
        # Calculate the denominator
        denominator = 0
        for i in observed_values:
            if i == 1:
                denominator += 1
        # Calculate the recall and precision        
        n_positive = 0
        n_negative = 0
        precision = []
        recall = []
        for match in match_up:
            if match[1] == 1:
                n_positive += 1
            else:
                n_negative += 1
            precision.append(n_positive/(n_positive + n_negative))
            recall.append(n_positive/denominator)
        list_recall_precision.append([copy.deepcopy(recall), copy.deepcopy(precision)])
        
    parameter['list_recall_precision'] = list_recall_precision
########################################################################################3    

  
    
#####################################################################################
# Record the time

def main(wd_results, wd_negative_samples, header):
    os.chdir(wd_results)
    with open(header+'_aa_train', 'r') as f:
        positive_training_set = json.load(f)
    with open(header+'_aa_test', 'r') as f:
        positive_testing_set = json.load(f)
        
    os.chdir(wd_negative_samples)
    with open(header+'_train_negative', 'r') as f:
        negative_training_set = json.load(f)
    with open(header+'_test_negative', 'r') as f:
        negative_testing_set = json.load(f)  
    
    # set the observed value as 1 for the positive trainging set.    
    for parepi in positive_training_set:
        parepi[2] = 1
    for parepi in positive_testing_set:
        parepi[2] = 1
    

    
    training_set = copy.deepcopy(positive_training_set)
    training_set.extend(negative_training_set)
    # Remove the duplicates
    centers, non_redundent_training_set = Remove_duplicates(training_set)
    # Calculate the distance matrix and the total design matrix
    print('Calculating the distance matrix')
    distance_matrix = Distance_matrix(non_redundent_training_set, mode='Multiplication') # mode is an ajustable parameter
    design_matrix = Design_matrix(distance_matrix, basis_function = 'Markov', radius = 1) # basis function is adjustable

    # Remove the duplicates
    # Load the parameter
    parameter = {}
    parameter['design_matrix'] = design_matrix
    parameter['centers'] = centers
    parameter['non_redundent_training_set'] = non_redundent_training_set
    parameter['training_set'] = training_set
    len(centers)
    # Initiate the coeff and reg
    parameter['coeff'] = np.zeros((len(centers)+1,1))
    parameter['reg'] = np.ones((len(centers)+1,1))
        
    # Generate the observed values
    observed_values = np.zeros((len(non_redundent_training_set ), 1))
    for i in range(len(non_redundent_training_set)):
        observed_values[i,0]=non_redundent_training_set [i][2]
    # Load the parameter
    parameter['observed_values'] = observed_values
    len(observed_values)
    # Generate the testing design matrix
    testing_set = copy.deepcopy(positive_testing_set)
    testing_set.extend(negative_testing_set)

    percentage = [0.5, 0.2,0.1,0.05,0.02,0.01,0.008,0.005,0.002,0.001]
    list_n_centers = []
    for i in percentage:
        n= math.floor(len(non_redundent_training_set)*i)
        if n>= 1:
            list_n_centers.append(n)
        
    parameter['list_n_centers'] = list_n_centers
    
    parameter = Coeff_select_centers(parameter)



    testing_set = copy.deepcopy(positive_testing_set)
    testing_set.extend(negative_testing_set)
    observed_values = []
    for i in testing_set:
        observed_values.append(i[2])
    Recall_precision_calculation(parameter, testing_set, observed_values)
    
    return parameter
#len(parameter['list_recall_precision'])
#list_recall_precision = parameter['list_recall_precision']
#os.chdir("/home/leo/Documents/Database/Pipeline/Results/1_free")
#with open('2_2_test_results')as f:
#    coverage_results = json.load(f)
#########################################################################
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
######################################################
ordered_headers = ['2_2','3_3', '2_3', '3_2', '3_4', '4_3','2_4', '4_2','4_4',\
                   '1_1','1_2', '2_1', '1_3', '3_1', '1_4', '4_1']
############################################################
if __name__ == '__main__':
    directory_paires = [['/home/leo/Documents/Database/Pipeline/Results/1_free',\
      '/home/leo/Documents/Database/Pipeline/Negative samples/1_free'],\
     ['/home/leo/Documents/Database/Pipeline/Results/0_free',\
                           '/home/leo/Documents/Database/Pipeline/Negative samples/0_free']]
    
    # Let check if all the corresponding length of the positive training and the negative training
    # are the same
    for i in range(2):
        wd_results = directory_paires[i][0]
        wd_negative_samples = directory_paires[i][1]
        n = 0
        for i in range(1,5):
            for j in range(1, 5):
                header = str(i)+'_'+str(j)
                os.chdir(wd_results)
                with open(header+'_aa_train', 'r') as f:
                    positive_training_set = json.load(f)
                    
                os.chdir(wd_negative_samples)
                with open(header+'_train_negative', 'r') as f:
                    negative_training_set = json.load(f)
                
                if len(positive_training_set) != len(negative_training_set):
                    print(header+' positive and negative are not of the same length')
                    print(wd_results)
                    n = 1
                    break
    if n == 0:
        print('Now, we are going to the large amount of data processing. It may '+\
              'spend a lot of time.')
            
    # Get into the work
    # Fist let us set the working oder of different files. working on the most interesting file 
    # first then gradually reduce to the most uninteresting file, so the even if there are some
    # problems pop up in the middle, we have finished the most interesting files.
    for i in range(2):
        wd_results = directory_paires[i][0]
        wd_negative_samples = directory_paires[i][1]
        for header in ordered_headers:
            print('Working on '+ header)
            results =  main(wd_results, wd_negative_samples, header)
            
             # Save the results
            os.chdir(wd_results)
            with open(header+'_test_results_Coeff_RBFN', 'w') as f:
                json.dump(results, f,  cls=NumpyEncoder)
############################################################
#os.chdir("/home/leo/Documents/Database/Pipeline/Results/1_free")
#with open('2_2_test_results_coeff_RBFN','r') as f:
#    test_results_coeff_RBFN = json.load(f)
#test_results_coeff_RBFN.keys()
#coverage_recall = test_results_coeff_RBFN['recall']
#list_recall_precision = test_results_coeff_RBFN['list_recall_precision']
#len(test_results_coeff_RBFN['centers'])
#for i in test_results_coeff_RBFN['list_n_centers']:
#    print(i)
#for i in range(len(list_recall_precision)):
#    plt.plot(list_recall_precision[i][0], list_recall_precision[i][1])
#    plt.show()
#####################################################################################





