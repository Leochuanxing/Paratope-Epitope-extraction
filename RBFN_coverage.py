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

#############################################################
# Define  distances
from Bio import Align
from Bio.SubsMat.MatrixInfo import blosum62
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = blosum62
aligner.open_gap_score = -5
aligner.extend_gap_score = -1
aligner.mode = 'global'

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

def Multiplication_distance(Ab_seq1, Ab_seq2, Ag_seq1, Ag_seq2):
    l_Ab = len(Ab_seq1)
    l_Ag = len(Ag_seq1)
    distance = (1 - (4*l_Ab + aligner.score(Ab_seq1, Ab_seq2))/(15*l_Ab)) *  (1 - (4*l_Ag + aligner.score(Ag_seq1, Ag_seq2))/(15*l_Ag))
    return distance
    
def Addition_distance(Ab_seq1, Ab_seq2, Ag_seq1, Ag_seq2):
    l_Ab = len(Ab_seq1)
    l_Ag = len(Ag_seq1)
    distance = (1 - (4*l_Ab + aligner.score(Ab_seq1, Ab_seq2))/(15*l_Ab)) + (1 - (4*l_Ag + aligner.score(Ag_seq1, Ag_seq2))/(15*l_Ag))
#    distance = -(aligner.score(Ab_seq1, Ab_seq2) + aligner.score(Ag_seq1, Ag_seq2))/(l_Ab*11 + l_Ag * 11)
    return distance       
#############################################################
'''
Generate the appropriate distance matrix
'''
def Distance_matrix(samples, mode = 'Multiplication'):
    n = len(samples)
    distance_matrix = np.zeros((n,n))
    for i in range(n):
        for j in range(i, n):
            if mode == 'Addition':
                distance_matrix[i,j] = Addition_distance(To_seq(samples[i][0]),To_seq(samples[j][0]),
                               To_seq(samples[i][1]),To_seq(samples[j][1]))
                distance_matrix[j,i] =  distance_matrix[i,j]
            if mode == 'Multiplication':
                distance_matrix[i,j] = Multiplication_distance(To_seq(samples[i][0]),To_seq(samples[j][0]),
                               To_seq(samples[i][1]),To_seq(samples[j][1]))
                distance_matrix[j,i] =  distance_matrix[i,j]
    # Lets scale the distance matrix. No, do not scale here, scale the design matrix
#    std = np.std(distance_matrix, axis = 0, keepdims=True)
#    distance_matrix = distance_matrix / np.reshape(std, (1, n))
    return distance_matrix
#a = np.arange(0, 12, 1)
#np.reshape(a, (3, 4))
###################################################################################

'''
Pruning:
    to return the prunned centers
Input: 
      distance_matrix: 
          the distance matrix of all the training samples. It is better that the 
          distance_matrix is normalized or scaled by dividing the standard deviation.
      cover_radius: 
          gives the radius a center can cover
Output: 
    centers:
        in the same form as the training set
    samples_to_centers_matrix:
        A matrix gives the distance between the samples and the centers, it 
        can be used directly to calculate the design matrix            
'''

def Generate_coverage_centers(distance_matrix, cover_radius):
    n = np.shape(distance_matrix)[0]
    centers = []
    left_over = list(range(n))
    while len(left_over) != 0:
        centers.append(left_over[0])
        to_be_removed = []
        left_over.remove(left_over[0])
        for i in left_over:
            if distance_matrix[i, centers[-1]] <= cover_radius:
                to_be_removed.append(i)
        for j in to_be_removed:
            left_over.remove(j)
    # Do we need to sort the centers? I don't think so.
#    print(len(centers))
#    samples_to_centers_matrix = distance_matrix[:, centers]
    return centers
#####################################################################################

# Define Radial Basis Functions
'''
Gaussian: 
    The Gaussian radial basis function
Inputs:
    distance: 
        a positive real number
    radius:
        float, a real positive number, gives the scaling factor of the distance between x and the center
        
    distance_mode:
        a string, takes the value of either 'Multiplication_distance' or 'Addition_distance',
        Gives different way of combining the Ab_Ab distances and the Ag_Ag distances.        
'''
def Gaussian(distance, radius):
    return math.exp(-distance**2/radius**2)
    
def Mrakov(distance, radius):
    return math.exp(-distance/radius)
'''
beta is a positive number
'''
def Inverse_Multi_Quadric(distance,c, beta):
    return (distance**2 + c**2)**(-beta)

def Design_matrix(distance_matrix, basis_function = 'Markov', radius = 1):
    nrow = np.shape(distance_matrix)[0]
    ncol = np.shape(distance_matrix)[1]
    
    design_matrix = np.zeros_like(distance_matrix)
    for i in range(nrow):
        for j in range(ncol):
            if basis_function == 'Gaussian':
                design_matrix[i, j] = Gaussian(distance_matrix[i, j], radius)
            elif basis_function == 'Markov':
                design_matrix[i, j] = Mrakov(distance_matrix[i, j], radius)
            elif basis_function == 'Inverse_Multi_Quadric':
                design_matrix[i, j] = Inverse_Multi_Quadric(distance_matrix[i, j], c=1, beta=2)  
    # Since we will use the linear conbination of this design matrix, it is reasonable to normalize it
    average = np.average(design_matrix, axis= 0)
    design_matrix = design_matrix - average
    std = np.std(design_matrix, axis = 0)
    design_matrix = design_matrix / std
    # Should we add a constant column?  Lets add.
    design_matrix = np.hstack((design_matrix, np.ones((nrow,1))))
           
    return design_matrix
#####################################################################################

##############################################################################
'''
Loss:
    To return the loss, the gradient, and the parameters
Input:
        The return value of Design_matrix, in the shape of (m,n)
    observed_values:
        A vector gives the observed valuse, in the shape of (m,1)
    parameter:
        A dictionary, 
        parameter['coeff'] contains the vector of the coefficients, in the shape of (n,1)
        parameter['reg'] contains the vector of the regularization coefficients, in the shape of (n,1)
        Here we had better design this loss function so that it can be used in the 
        ridge selection of the centers as well.        
Output:
    loss:
        A real number
    gradient:
        a dictionary
        gradient['coeff'] is the gradient of the parameter['coeff']
        gradient['reg']  is the gradient of the parameter['reg']
'''
def Loss(design_matrix, observed_values, parameter):
    # Unpack the dictionary 
    coeff = parameter['coeff']
    reg = parameter['reg']
    
    coeff_square = coeff*coeff
    diff = design_matrix.dot(coeff) - observed_values
    loss = (diff.T).dot(diff)
    loss += (coeff_square.T).dot(reg)
    grad_coeff = 2*(design_matrix.T).dot(diff) + 2 * coeff * reg
    # Pack up the results
    gradient = {}
    gradient['coeff'] = grad_coeff
    
    return loss, gradient
####################################################################################

##################################################################################
'''
Train_RBFN_BFGS:
    A function to train the RBFN by using the BFGS method
Input:
    design_matrix:
        The return value of Design_matrix, in the shape of (m,n)
    observed_values:
        A vector gives the observed valuse, in the shape of (m,1)
    rho:
        The contraction factor in Ramijo backtracking
    c:
        the Goldstein coefficient c
    termination:
        The termination condition with the norm of the gradient < termination
Output:
    parameter:
        A dictionary contains
        parameter['coeff'], the coefficient after training
        parameter['reg'], the regularization coefficient    
'''


def Train_RBFN_BFGS(design_matrix, observed_values, rho=0.9, c = 1e-3, termination = 1e-2,\
                    parameter_inheritance = False, parameter=None):
    
    nrow = np.shape(design_matrix)[0]
    ncol = np.shape(design_matrix)[1]
        
    # Give the initial Hessian H. The variables are the coeff and reg
    H = np.eye(ncol)/(10*nrow)
    # Check if it inherit the parameter from somewhere else.
    if not parameter_inheritance :
        # Set the starting point
        coeff = np.zeros((ncol,1))
        parameter = {}
        parameter['coeff'] = coeff
        #The reg should not be negative. It is better that reg > delta, a small positive number
        reg = np.ones((ncol,1)) * 1
        parameter['reg'] = reg

    # BFGS algorithm
    loss, gradient = Loss(design_matrix, observed_values, parameter)
    grad_coeff = gradient['coeff']
    grad_square = (grad_coeff.T).dot(grad_coeff)
    while grad_square >= termination**2:        
        p = - H.dot(grad_coeff)        
        # Find the next coeff
        parameter_new = {}
        parameter_new['coeff'] = p + parameter['coeff']
        parameter_new['reg'] = parameter['reg']
        
        new_loss, new_gradient = Loss(design_matrix, observed_values, parameter_new)
        # Ramijo Back-tracking
        while new_loss > loss + c * (grad_coeff.T).dot(p):
            p *= rho
            parameter_new['coeff'] = p + parameter['coeff']            
            new_loss, new_gradient = Loss(design_matrix, observed_values, parameter_new)
        
        # update H
        s = p
        new_grad = new_gradient['coeff']
        y = new_grad - grad_coeff
        r = (y.T).dot(s)
        if r != 0:
            r = 1/r
            I = np.eye(ncol)
            H = (I - r*s.dot(y.T)).dot(H).dot(I - r*y.dot(s.T)) + r*s.dot(s.T)
        else:
            H = I
        # Update loss, grad, grad_square and paramter
        loss = new_loss
        grad_coeff = new_grad
        parameter['coeff'] = parameter_new['coeff']
        grad_square = (grad_coeff.T).dot(grad_coeff)
        print('loss  ', loss, '    ','grad_square   ', grad_square)
    return parameter, loss
###########################################################################################################

###############################################################################
'''
Generating the testing design matrix
'''
#len(training_set)
def Generate_testing_design_matrix(testing_set, centers, training_set, mode = 'Multiplication', basis_function = 'Markov'):
    centers_parepi = [training_set[i] for i in centers]
    nrow = len(testing_set)
    ncol = len(centers_parepi)
    distance_matrix = np.zeros((nrow,ncol))
    for i in range(nrow):
        for j in range(ncol):
            if mode == 'Addition':
                distance_matrix[i,j] = Addition_distance(To_seq(testing_set[i][0]),To_seq(centers_parepi[j][0]),
                               To_seq(testing_set[i][1]),To_seq(centers_parepi[j][1]))
            if mode == 'Multiplication':
                distance_matrix[i,j] = Multiplication_distance(To_seq(testing_set[i][0]),To_seq(centers_parepi[j][0]),
                               To_seq(testing_set[i][1]),To_seq(centers_parepi[j][1]))
                
    testing_design_matrix = Design_matrix(distance_matrix, basis_function, radius = 1)
    
    return testing_design_matrix

#a = [1, 3, 5]
#b = list(range(15))
#c = [b[i] for i in a]
#c
################################################################################
def Prediction(testing, parameter, testing_design_matrix):
    coeff = parameter['coeff']
    prediction = testing_design_matrix.dot(coeff)
    n = 0
    for i in range(len(testing)):
        if prediction[i,0] > 0 and testing[i][2] == 1:
            n += 1
        elif prediction[i,0] < 0 and testing[i][2] == -1:
            n += 1
    rate = n / len(testing)
    return round(rate, 3)

#rate = Prediction(testing_set, parameter, testing_design_matrix)

#################################################################################
'''
Define a small function
Input:
    one_testing_sample:
        is a list in the form of [Ab, Ag, m, n]
    prediction_pool:
        A list in the form of [[Ab, Ag, m, n], ....], where our prediction is chosen from.
        Lets assume there may be redundency in this prediction_pool.
    prediction_direction:
        a string, either 'Ab_to_Ag' or 'Ag_to_Ab'
Output:
    combination:
        if the prediction_direction is 'Ab_to_Ag', the conbination is all possible combinations
        between the Ab from the one_testing_sample and the Ag from the prediction_pool
        
        if the prediction_direction is 'Ag_to_Ab', the conbination is all possible combinations
        between the Ag from the one_testing_sample and the Ab from the prediction_pool.
        
    closest_index:
        an integer, gives the position index of the element in the combination with 
        the closest Ag to the one_testing_sample Ag if the prediction_direction is 'Ab_to_Ag'.
        
        Otherwise it gives the position index of the element in the combination with the 
        closest Ab to the one_testing_sample Ab.        
'''
def Pre_concrete_prediction(one_testing_sample, prediction_pool, prediction_direction= 'Ag_to_Ab'):
    # Lets creat the nonredundent pool
    Ag_pool = []
    Ab_pool = []
    for parepi in prediction_pool:
        if parepi[0] not in Ab_pool:
            Ab_pool.append(parepi[0])
        if parepi[1] not in Ag_pool:
            Ag_pool.append(parepi[1])
    # Find the combination and the closest index
    testing_Ab = one_testing_sample[0]
    testing_Ag = one_testing_sample[1]
    combination = []
    largest_similarity = -100000
    if prediction_direction == 'Ag_to_Ab':
        for i in range(len(Ab_pool)):
            combination.append([Ab_pool[i], testing_Ag])
            score = aligner.score(To_seq(testing_Ab), To_seq(Ab_pool[i]))
            if score > largest_similarity:
                largest_similarity = score
                closest_index = i
    elif prediction_direction == 'Ab_to_Ag':
        for i in range(len(Ag_pool)):
            combination.append([testing_Ab, Ag_pool[i]])
            score = aligner.score(To_seq(testing_Ag), To_seq(Ag_pool[i]))
            if score > largest_similarity:
                largest_similarity = score
                closest_index = i 
    return combination, closest_index


#######################################################################################
'''
THis function is to generate the prediction_pool\

Input:
    select_from:
        a string, takes the value of either 'positive_centers', 'positive_training_set'
        or 'all_possible'
    Ab_length:
        an integer, gives the length of the Ab
    Ag_length:
        an integer, gives the length of the Ag
    centers:
        a list of integer, gives the index of the center in the list of training_set
    positive_training_set:
        the same as above
Output:
    prediction_pool:
        a list as an input of Pre_concrete_prediction
'''
def Generate_prediction_pool(select_from, Ab_length, Ag_length,\
                              positive_training_set):

    if select_from == 'positive_training_set':
        prediction_pool = positive_training_set
    elif select_from == 'all_possible':
        Ab_all_possible = All_possible(Ab_length)
        Ag_all_possible = All_possible(Ag_length)
        n = max(len(Ab_all_possible), len(Ag_all_possible))
        m = min(len(Ab_all_possible), len(Ag_all_possible))
        if len(Ab_all_possible) == n:
            for i in range(n):
                if i >= m:
                    j = i - m
                else:
                    j = i
                prediction_pool.append([Ab_all_possible[i], Ag_all_possible[j]])
                
        if len(Ag_all_possible) == n:
            for i in range(n):
                if i >= m:
                    j = i - m
                else:
                    j = i
                prediction_pool.append([Ab_all_possible[j], Ag_all_possible[i]])
    return prediction_pool

#############################################################################
'''
Define a function to generate the all possible combinations with the given lenth
'''

def All_possible(length, all_possible = []):
    # 20 aa
    AA=['TYR', 'LYS', 'ASP', 'ASN', 'TRP', 'PHE', 'GLN', 'GLU', 'PRO','GLY',\
    'THR', 'SER','ARG','HIS','LEU','ILE', 'CYS','ALA','MET','VAL']
    
    if all_possible == []:
        for i in AA:
            all_possible.append([i])
    new_all_possible = []
    if length > 1:
        for i in AA:
            for j in all_possible:
                new = copy.deepcopy(j)
                new.append(i)                
                if new not in new_all_possible:
                    new_all_possible.append(new)
        length -= 1
        return All_possible(length, new_all_possible)
    else:
        return all_possible

#############################################################################3


'''
Pre_processing:
    a function to processing the data, so that they can be easily used latter
Input:
    positive_testing_set
    positive_training_set
    negative_testing_set
    negative_training_set
Output:
    data_dict:
        a dictionary, contains
        data_dict['training_set']
        data_dict['testing_set']
        data_dict['observed_values']
        data_dict['design_matrix']
        data_dict['distance_matrix']
'''

def Pre_processing(positive_training_set, positive_testing_set, negative_training_set, negative_testing_set):    
# set the observed value as 1 for the positive trainging set.    
    for parepi in positive_training_set:
        parepi[2] = 1
    for parepi in positive_testing_set:
        parepi[2] = 1
        
    training_set = copy.deepcopy(positive_training_set)
    training_set.extend(negative_training_set)
    
    #len(training_set)
    # Calculate the distance matrix and the total design matrix
    distance_matrix = Distance_matrix(training_set, mode='Multiplication') # mode is an ajustable parameter
    design_matrix = Design_matrix(distance_matrix, basis_function = 'Markov', radius = 1) # basis function is adjustable
    #distance_matrix.shape
    #design_matrix.shape
    # Generate the observed values
    observed_values = np.zeros((len(training_set), 1))
    for i in range(len(training_set)):
        observed_values[i,0]=training_set[i][2]    
      
    # Generate the testing design matrix
    testing_set = copy.deepcopy(positive_testing_set)
    testing_set.extend(negative_testing_set)
    
    # Load the results to the data_dict
    data_dict = {}
    data_dict['training_set'] = training_set
    data_dict['testing_set'] = testing_set
    data_dict['observed_values'] = observed_values
    data_dict['design_matrix'] = design_matrix
    data_dict['distance_matrix'] = distance_matrix
    
    return data_dict
    
#####################################################################################
 

######################################################################################
    
def Coverage_centers_binary_prediction(distance_matrix, design_matrix, training_set, testing_set, observed_values, n_centers): 

    # Generate the centers
    centers_to_centers_matrix = copy.deepcopy(distance_matrix)  
    
    all_centers = []
    m = len(training_set) + 1# Make sure the it goes to the while loop for the first n in n_centers 
    cover_radius = 0
    for n in n_centers:
        while m > n:
            centers = Generate_coverage_centers(centers_to_centers_matrix, cover_radius)
            centers_to_centers_matrix = distance_matrix[centers, :][:, centers]
            cover_radius += 0.001# the step increase of the centers, this stepsize can be adjusted
            m = len(centers)
    #        print(m)
            if m <= n:
                selected_centers = copy.deepcopy(centers)
                all_centers.append(selected_centers)
    # Calculate the prediction accuracy rate, parameters, testing design matrix, and the centers
    # and store all of them in a dictionary called results
    all_rate = []
    all_parameters=[]
    for i in range(len(all_centers)):
        centers = all_centers[i]
        # Calculate the design_matrix
        centers.append(-1)
        design_matrix1 = design_matrix[:, centers]
    #    design_matrix1.shape
        centers.remove(centers[-1])   
        # Train
        parameter, loss = Train_RBFN_BFGS(design_matrix1, observed_values, rho=0.9, c = 1e-3, termination = len(centers)/1000,\
                            parameter_inheritance = False, parameter=None) 
        #Calculate testing_design_matrix
        testing_design_matrix = Generate_testing_design_matrix(testing_set, centers,\
                                                           training_set, mode = 'Multiplication', basis_function = 'Markov')
        # Binary prediction
        rate = Prediction(testing_set, parameter, testing_design_matrix)
        all_rate.append(rate)
        # Load the parameters
        stored_parameter = copy.deepcopy(parameter)
        all_parameters.append(stored_parameter)

           
    all_rate
    results = {}
    results['mode'] = 'Multiplication_Markov'
    results['n_centers'] = n_centers
    results['rate'] = all_rate
    results['all_parameters'] = all_parameters
    results['all_centers']= all_centers
    
    return results
################################################################################



####################################################################
'''
Generate_cross_testing_training:
    A function, generates the testing and training positions for the given number of 
    cross validation. The return values can be used to extract the design matrix from
    the big design matrix
Input:
    training_set:
        a list contains all the training set
    cross_number:
        an integer, gives the number of crosses
Output:
    cross_training_set:
        a list contains 'cross_number' of lists of positions of the training_set
    cross_testing_set:
        a list contains 'cross_number' of lists of positions of the training_set
        
    The above two list should be corresponding with each other, that is they are 
    ordered correspondingly.
'''
def Generate_cross_testing_training(training_set, cross_number):
    cross_training_set = []
    cross_testing_set = []
    n_total = len(training_set)
    size = math.floor(n_total/cross_number)
    # Cut into pieces
    for i in range(cross_number-1):
        lower=i*size; upper = (i+1) * size
        cross_testing_set.append(list(range(lower, upper)))
    # Deal with the last block
    lower = size * (cross_number-1)
    cross_testing_set.append(list(range(lower, n_total)))
    # Load the cross_training_set
    for one_test_set in cross_testing_set:
        one_train_set = []
        for j in range(n_total):
            if j not in one_test_set:
                one_train_set.append(j)
        cross_training_set.append(one_train_set)   
            
    return cross_testing_set, cross_training_set

#a = list(range(10))
#cross_test, cross_train = Generate_cross_testing_training(a, 3)
#cross_test
#cross_train
######################################################################################
'''
Stack:
    this function is to adjust the returned values from the Generate_cross_testing_training, 
    so that the indices is consistant with the design_matrix, observed_values, and the training_set
Input:
    positive_cross_test: the returned value of Generate_cross_testing_training(positive_training_set, cross_number)
    positive_cross_train: As above
    negative_cross_test: the returned value of Generate_cross_testing_training(negative_training_set, cross_number)
    negative_cross_train: As above
    ith: an integer, gives with validation set is used, in range(cross_number)
Output:
    cross_train_indices: A list of indices correspondint to the design_matrix, observed_values, and the training_set
    cross_test_indices: the same as above
'''
def Stack(positive_cross_test, positive_cross_train, negative_cross_test, negative_cross_train, ith):
    # Raise the negative_cross_training and negative_cross_testing
    n=len(positive_cross_test[ith])+len(positive_cross_train[ith])
    
    raised_negative_test = (np.array(negative_cross_test[ith]) + n).flatten().tolist()
    
    raised_negative_train = (np.array(negative_cross_train[ith]) + n).flatten().tolist()
    
    # Stack the indices
    cross_train_indices = copy.deepcopy(positive_cross_train[ith])
    cross_train_indices.extend(raised_negative_train)
    
    cross_test_indices = copy.deepcopy(positive_cross_test[ith])
    cross_test_indices.extend(raised_negative_test)
    return cross_train_indices, cross_test_indices
###########################################################################################

'''
Cross_pre_processing:
    Pre processing the cross validation data
Input:
    cross_train_indices: The return value of Stack
    cross_test_indices: The return value of Stack
    design_matrix: the design matrix of training set
    data_dict:
        a dictionary, contains
        data_dict['training_set']
        data_dict['design_matrix']
        data_dict['distance_matrix']
        
 
Output:
    cross_data_dict: a dictionary, contains:
        
        cross_data_dict]['design_matrix']:
        design matrix corresponding to the cross_train, which is the combination 
                  one_positive_cross_train and one_negative_cross_train
                  
        cross_data_dict['observed_values']:
        The observed values corresponding to the cross_design_matrix
        
        cross_data_dict['testing_set']
        Tha parepis correponding to the cross_test_indices
        
        cross_data_dict['training_set']
        The parepies corresponding to the cross_train_indices
        
        cross_data_dict['distance_matrix']
        The distance matrix corresponding to the cross_data_dict['training_set']
            
'''
def Cross_pre_processing(cross_train_indices, cross_test_indices, data_dict):
    # Take out the data from the data_dict
    training_set = data_dict['training_set']
    design_matrix = data_dict['design_matrix']
    distance_matrix = data_dict['distance_matrix']
    # Get the cross_design_matrix
    cols = copy.deepcopy(cross_train_indices)
    cols.append(-1)
    cross_design_matrix = design_matrix[cross_train_indices,:][:,cols]
    # Generate the cross distance matrix
    cross_distance_matrix = distance_matrix[cross_train_indices,:][:, cross_train_indices]
    # Get the cross_test
    cross_test = []
    for i in cross_test_indices:
        cross_test.append(training_set[i])
    # Get the cross_train
    cross_train = []
    for i in cross_train_indices:
        cross_train.append(training_set[i])
    # Get the observed values
    cross_observed_values = []
    for i in cross_train_indices:
        cross_observed_values.append(training_set[i][2])
    cross_observed_values = np.array(cross_observed_values).reshape((-1, 1))
        
    # Load the values to the cross_data_dict:
    cross_data_dict ={}
    cross_data_dict['design_matrix'] = cross_design_matrix
    cross_data_dict['observed_values'] = cross_observed_values
    cross_data_dict['testing_set'] = cross_test
    cross_data_dict['training_set'] = cross_train
    cross_data_dict['distance_matrix'] = cross_distance_matrix
    
    return cross_data_dict


#################################################################################
## Use a samll samples set to test the above function
#positive_training = positive_training_set[:10]
#negative_training = negative_training_set[:10]
#training_set = copy.deepcopy(positive_training)
#training_set.extend(negative_training)
#
#small_distance_matrix = Distance_matrix(training_set, mode='Multiplication') # mode is an ajustable parameter
#design_matrix = Design_matrix(small_distance_matrix, basis_function = 'Markov', radius = 1)
#
#cross_number = 3
#positive_cross_test, positive_cross_train = Generate_cross_testing_training(positive_training, cross_number)
#negative_cross_test, negative_cross_train = Generate_cross_testing_training(negative_training, cross_number)
#
#ith = 1
#cross_train_indices, cross_test_indices = Stack(positive_cross_test, positive_cross_train,\
#                                                negative_cross_test, negative_cross_train, ith)
#cross_design_matrix, cross_observed_values, cross_test, cross_train = \
#Cross_pre_processing(cross_train_indices, cross_test_indices, design_matrix, training_set)
#
##cross_design_matrix
#cross_observed_values
#cross_test
#cross_train
#
#for i in range(len(cross_train_indices)):
#    indi=cross_train_indices[i]
#    for j in range(len(cross_train_indices)):
#        indj = cross_train_indices[j]
#        if design_matrix[indi, indj] != cross_design_matrix[i,j]:
#            print('different')
#            break
#        
#cross_train_indices
#cross_test_indices
#training_set
########################################################################
'''
Generate_coverage_centers:
    This function is to generate the coverage centers
Input:
    cross_data_dict
    cross_parameter:
        right now, it contains
        cross_parameter['center_percentages']:
            gives a list of the percentages the number of centers to the training_set
Output:
    cross_parameter['centers']:
        a list, contains lists of centers
        
'''
def Cross_coverage_centers(cross_data_dict, cross_parameter):
    # Take out the values
    distance_matrix = cross_data_dict['distance_matrix']
    training_set = cross_data_dict['training_set']
    
    center_percentages = cross_parameter['center_percentages']
    
    # Calculate the number of centers
    n_centers = []
    for percentage in center_percentages:
        n_centers.append(math.floor(len(training_set)*percentage))
        
    # Creat the initial centers_to_centers_matrix for coverage center selection
    centers_to_centers_matrix = copy.deepcopy(distance_matrix)  
    # Start the process of selecting centers
    all_centers = []
    m = len(training_set) + 1# Make sure the it goes to the while loop for the first n in n_centers 
    cover_radius = 0
    for n in n_centers:
        while m > n:
            centers = Generate_coverage_centers(centers_to_centers_matrix, cover_radius)
            centers_to_centers_matrix = distance_matrix[centers, :][:, centers]
            cover_radius += 0.001# the step increase of the centers, this stepsize can be adjusted
            m = len(centers)
            if m <= n:
                selected_centers = copy.deepcopy(centers)
                all_centers.append(selected_centers)
    # Load the results
    cross_parameter['centers'] = all_centers
    
    return cross_parameter
########################################################################################
'''
Get_cross_coeff:
    train the model, to get the coefficient for the validation data.
Input:
    cross_parameter
    cross_data
Ouput:
    cross_parameter['all_coeff']:
        a list of coeff arrays, corresponding to cross_parameter['centers']
'''
def Get_cross_coeff(cross_parameter, cross_data_dict):
    # Take out the values
    all_centers = cross_parameter['centers']
    observed_values = cross_data_dict['observed_values'] 
    design_matrix = cross_data_dict['design_matrix']
    # Begin the training
    all_parameters=[]
    for i in range(len(all_centers)):
        centers = all_centers[i]
        # Calculate the design_matrix
        centers.append(-1)
        design_matrix_selected = design_matrix[:, centers]
        centers.remove(-1) 
#        print(len(centers))
        # Train
        parameter, loss = Train_RBFN_BFGS(design_matrix_selected, observed_values,\
                                          rho=0.9, c = 1e-3, termination = len(centers)/1000,\
                                          parameter_inheritance = False, parameter=None) 
 
        # Load the parameters
        stored_parameter = copy.deepcopy(parameter)
        all_parameters.append(stored_parameter)
        
        cross_parameter['all_coeff'] = all_parameters
    
    return cross_parameter
##########################################################################################
'''
Hyperparameter_evaluation:
    A function to calculate the scaled area between the recall-precision curve and 
    the line with precision = 0.5. 
Input:
    cross_data_dict:
    cross_parameter
Output:
    cross_parameter['areas']
      A list, gives the areas for different set of coefficient in cross_parameter['all_coeff']
'''

def Hyperparameter_evaluation(cross_data_dict, cross_parameter):
    # Take out the values
    training_set = cross_data_dict['training_set']
    testing_set = cross_data_dict['testing_set']
    centers = cross_parameter['centers']
    all_coeff = cross_parameter['all_coeff']
    
    # Iterate over different set of centers
    areas = []
    for i in range(len(centers)):
        one_centers = centers[i]
        
        # Find the predicted values
        testing_design_matrix = Generate_testing_design_matrix(testing_set, one_centers, training_set,\
                                       mode = 'Multiplication', basis_function = 'Markov')
        coeff = np.array(all_coeff[i]['coeff']).reshape((-1,1))
        prediction = testing_design_matrix.dot(coeff)
        
        # Related the observed testing values with the predicted values
        pred_observed = []
        for i in range(len(testing_set)):
            pred_observed.append([prediction[i,0], testing_set[i][2]])
        pred_observed.sort(key=lambda x:x[0], reverse=True)
        
        #Calculate the area
#        precision=[]
#        recall = []
#        denominator = 0.5 * len(testing_set)
        area = 0
        n_positive = 0
        n_negative = 0
        for i in pred_observed:
            if i[1] == 1:
                n_positive += 1
            else:
                n_negative += 1
            area += n_positive/(n_positive+n_negative)
#            precision.append(n_positive/(n_positive+n_negative))
#            recall.append(n_positive/denominator)
        areas.append(area)
    
    #Load the results
    cross_parameter['areas'] = areas
    
    return cross_parameter
#########################################################################
'''
Select_the_best_hyperparameter:
    This function is to select the best parameter accroding to the cross validation
Input:
    hyperparameters:
        a list, gives all the possible percentages
    cross_number: 
        an integer, gives the number of cross validation
    data_dict
Output:
    best_hyperparameter: a float, chosen from hperparameters
    all_cross_parameter: A dictionary of dictionaries, which are cross_parameters
'''
def Select_the_best_hyperparameter(hyperparameters, cross_number, data_dict):
    # Take out the values
    positive_training_set = data_dict['positive_training_set']
    negative_training_set = data_dict['negative_training_set']
    # Make sure the positive_training set and the negative training set are of 
    # the same length
    if len(positive_training_set) != len(negative_training_set):
        print('Positive and negative training set should be of the same length')
        return

    # Generate the cross validation data
    positive_cross_test, positive_cross_train =\
    Generate_cross_testing_training(positive_training_set, cross_number)
    
    negative_cross_test, negative_cross_train =\
    Generate_cross_testing_training(negative_training_set, cross_number)
    
    # Go to the different validations
    # Store all the parameters  in a dictionary called all_cross_parameters
    all_cross_parameter = {}
    # Record the areas
    total_area = np.zeros((cross_number, 1))
    for i in range(cross_number):
        # Load the hyperparameter, which is the center percentages to cross_parameter
        cross_parameter = {}
        cross_parameter['center_percentages'] = hyperparameters
        # From now on we are in one validation
        # Fist lets generate the cross_data_dict
        cross_train_indices, cross_test_indices =\
        Stack(positive_cross_test, positive_cross_train, negative_cross_test, negative_cross_train, i)
        
        cross_data_dict = Cross_pre_processing(cross_train_indices, cross_test_indices, data_dict)
        
        # Lets generate the centers in this cross validation
        print('Generating coverage centers')
        cross_parameter = Cross_coverage_centers(cross_data_dict, cross_parameter)
        
        # Lets train the model in this validation to get the coefficients for different hyperparameters
        cross_parameter = Get_cross_coeff(cross_parameter, cross_data_dict)
        
        # Evaluate the coefficent in this validation data set
        print('Evaluate the hyperparameters')
        cross_parameter = Hyperparameter_evaluation(cross_data_dict, cross_parameter)
        
        # We should take out the area, and add them up.
        all_cross_parameter[str(i)] = copy.deepcopy(cross_parameter)
        
        # Add up the area
        total_area += np.array(cross_parameter['areas'])
    
    # Select the best parameter according to the parameter
    best_hyperparameter = hyperparameters[np.argmax(total_area)]
    # Let us store the best_hyperparameter in all_cross_parameter
    all_cross_parameter['best_hyperparameter'] = best_hyperparameter
    return best_hyperparameter, all_cross_parameter
    


###########################################################################
   


########################################################################

'''
Define a main function, with input the working directory
'''

def main(wd_results, wd_negative_samples):
    os.chdir(wd_results)
    with open('2_2_aa_train', 'r') as f:
        positive_training_set = json.load(f)
    with open('2_2_aa_test', 'r') as f:
        positive_testing_set = json.load(f)
        
    os.chdir(wd_negative_samples)
    with open('2_2_train_negative', 'r') as f:
        negative_training_set = json.load(f)
    with open('2_2_test_negative', 'r') as f:
        negative_testing_set = json.load(f)  
    
    print('Pre processing the data')
    data_dict= Pre_processing(positive_training_set, positive_testing_set,\
                              negative_training_set, negative_testing_set)
    
    data_dict['positive_training_set'] = positive_training_set
    data_dict['negative_training_set'] = negative_training_set
    # set the hyper parameters for different length of training data
    
    if len(data_dict['training_set'])>=1000:
        hyperparameters = [0.005, 0.01, 0.02,0.04, 0.08, 0.16, 0.32, 0.64]
    if len(data_dict['training_set'])<1000:
        hyperparameters = [ 0.01, 0.02,0.04, 0.08, 0.16, 0.32, 0.64]
    # The percentages should be in a decreasint order, so that the code to find the coverage centers can work.
    hyperparameters.sort(reverse=True)
    # Set the cross number
    cross_number = 6
    
    best_hyperparameter, all_cross_parameter =\
    Select_the_best_hyperparameter(hyperparameters, cross_number, data_dict)
    
    return best_hyperparameter,all_cross_parameter
#best_hyperparameter
#all_cross_parameter.keys()
#type(all_cross_parameter['2']['all_coeff'])
#type(all_cross_parameter['2']['all_coeff'][0])
#len(all_cross_parameter['2']['all_coeff'])
#all_cross_parameter['2']['all_coeff'][0].keys()
#all_cross_parameter['2']['all_coeff'][0]['coeff']
###################################################################################
'''
Run the main and Save the results
'''
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
###################################################################################
        
if __name__ == '__main__':
    wd_results = '/home/leo/Documents/Database/Pipeline/Results/0_free'
    wd_negative_samples = '/home/leo/Documents/Database/Pipeline/Negative samples/0_free'
    
    best_hyperparameter,all_cross_parameter =  main(wd_results, wd_negative_samples)
   
# Save the results
    os.chdir(wd_results)
    with open('2_2_cross_results', 'w') as f:
        json.dump(all_cross_parameter, f,  cls=NumpyEncoder)
###################################################################################



