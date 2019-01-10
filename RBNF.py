###############################################################
# Import the modules
import random
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
#def Generate_random_negative(sample_size, data, Ab_lenth, Ag_length):
#    
#    TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
#                    ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
#                    ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]
#    AA = []
#    for aa in TripleSingle:
#        AA.append(aa[0])
#        
#    Ab_Ag = []
#    for parepi in data:
#        Ab_Ag.append([parepi[0], parepi[1]])
#        
#    negative_samples = []
#    while len(negative_samples) < sample_size:
#        r_Ab_r_Ag = []
#        while r_Ab_r_Ag == []:        
#            r_Ab = random.sample(AA, Ab_lenth)
#            r_Ag = random.sample(AA, Ag_length)
#            r_Ab_r_Ag  = [r_Ab, r_Ag]
#            if r_Ab_r_Ag in Ab_Ag:
#                r_Ab_r_Ag  = []
#        negative_samples.append([r_Ab, r_Ag, -1, -1])
#    return negative_samples
#####################################################################################

############################################################
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


################################################################
##################################################################
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
#
#design_matrix = Design_matrix(samples_to_centers_matrix, basis_function = 'Markov', radius = 1)    
#design_matrix.shape
#design_matrix[:5, :5]
#
#design_matrix = Design_matrix(distance_matrix, basis_function = 'Markov', radius = 1)
#design_matrix.shape
#import numpy as np
#a = np.arange(0, 16)
#design_matrix = np.reshape(a, (4, 4))
#average = np.average(design_matrix, axis= 0)
#design_matrix = design_matrix - average
#design_matrix
#std = np.std(design_matrix)
#design_matrix /= std
#design_matrix
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
Do a concrete predict. Don't just predict whether they contact or not, but also the 
exact paratope

Input:
    combination:
        the output of the pre concrete prediction
    complete_testing_design_matrix:
        a testing design_matrix generated by take all the combination as the input and 
        all the training_set as the centers
    centers:
        a list in the form of [1, 5, ...], where the elements indicate the position of those 
        centers in the training_set
    parameter: 
        the trained parameters corresponding to the centers       

Output:
    a boolean, True or False, indicates whether the prediction is right or wrong.
'''



def Top_x_percent_prediction(testing_design_matrix, parameter, closest_index, top_x_percent = 0.1):
    
    # Prediction
    coeff = parameter['coeff']
    coeff = np.array(coeff)
    prediction = testing_design_matrix.dot(coeff)
    prediction = prediction.flatten().tolist()
    # find the predictin score for the closest
    score_for_closest = prediction[closest_index]
    # find the lowest prediction score for the top 10 percent
    n_pool = testing_design_matrix.shape[0]
    n_top = math.floor(n_pool * top_x_percent)
    prediction.sort(reverse = True)
    n_top_score = prediction[n_top]
    # Tell whether our prediction is correct
    if score_for_closest >= n_top_score:
        return 1
    else:
        return 0
################################################################################################3

################################################################################
'''
The following block is to generate the centers through the ridge method.
The basic idea is:
    Fist lets fix the regulization coefficient the same for all the centers, then 
    we do the optimization with the regularization term part of the loss function.
    Then we check which coefficent is among the smallest ones, we delete the corresponding
    centers. 
    
    We repeat the process, until we get the required number of centers.
    
    Don't know whether the above method works or not. Let's try it out.
'''
'''
Define a step deletion function
Input:
    centers:
        a list of centers
    design_matrix:
        a design matrix in with all the training samples as centers.
        this design matrix is fixed and will be used again and again.
        
        Pay attention that the column number is one larger than the row number.
        
    observed_values:
        a list, in the same order of the row of the design_matrix
    reduce_number:
        An integer, gives how many centers should be deleted from the centers.
Output:
    centers:
        The length of this list is shorter than the input centers.
'''
def One_step_reduce_centers(centers, design_matrix, observed_values, parameter):
    # Don't forget to add the last column of the design matrix, which is for the biss.
    select = copy.deepcopy(centers)
    select.append(-1)
    new_design_matrix = design_matrix[:, select]
    # To speed up the search for centers
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

        
    parameter, loss = Train_RBFN_BFGS(new_design_matrix, observed_values, rho=0.85, c = 1e-3, termination=termination,\
                                      parameter_inheritance = True, parameter=parameter)
    coeff = parameter['coeff']
    # match the coeff up with the index of the centers
    
    coeff_list = np.abs(coeff)

    match_up = []
    for i in range(len(centers)):
        match_up.append([centers[i], coeff_list[i], i])
    # sort the match up according to the value of coeff
    match_up.sort(key = lambda x:x[1])
    # Remove the first one
    to_be_removed = match_up[0]
    # Remove the corresponding value from the center and the coeff
    centers.remove(to_be_removed[0])
    removed_coeff = coeff[to_be_removed[2], 0]
    coeff = np.delete(coeff, [to_be_removed[2]])
    coeff = coeff.reshape((-1, 1))
    # Return the parameter for further traiing and removed coeff value for control
    # and the centers for further deletion
    parameter['coeff'] = coeff
    parameter['reg'] = np.ones((len(coeff), 1))
    
    return centers, parameter, removed_coeff
    

###############################################################################

'''
Remove the duplicates before removing the centers
'''
def Remove_duplicates(training_set):
    centers = []
    pure_parepi = []
    for i in range(len(training_set)):
        parepi = training_set[i]
        if [parepi[0], parepi[1]] not in pure_parepi:            
            pure_parepi.append([parepi[0], parepi[1]])
            centers.append(i)
    return centers

#centers = Remove_duplicates(training_set)

#############################################################

'''
Input:
    control_coeff:
        float, remove all the centers with coeff larger than control_coeff
'''

def Coeff_select_centers(parameter,centers, training_set, design_matrix, observed_values, control_coeff, centers_inherit = False):
    # Initiate teh centers and teh removed_coeff
    if not centers_inherit:
        centers = Remove_duplicates(training_set)
        ncol = len(centers) + 1
        # Set the starting point
        coeff = np.zeros((ncol,1))
        parameter = {}
        parameter['coeff'] = coeff
        #The reg should not be negative. It is better that reg > delta, a small positive number
        reg = np.ones((ncol,1))  
        parameter['reg'] = reg

    #Make sure to add the last column
    removed_coeff = 0
    while abs(removed_coeff) < control_coeff:
        centers, parameter, removed_coeff = One_step_reduce_centers(centers, design_matrix, observed_values, parameter)
        print(removed_coeff, '   ', len(centers))
    cutoff_coeff = removed_coeff
#    print('cutoff coeff = d%', cutoff_coeff)# Here is to monitor the process.
    
    return centers, parameter, cutoff_coeff
##################################################################################################################

#n
#############################################################################################################################


os.chdir("/home/leo/Documents/Database/Pipeline/Ready_2_2_1_1")


with open('training_2_2_1_1_all_jump', 'r') as f:
    positive_training_set = json.load(f)
with open('testing_2_2_1_1_all_jump', 'r') as f:
    positive_testing_set = json.load(f)
with open('training_negative', 'r') as f:
    negative_training_set = json.load(f)
with open('testing_negative', 'r') as f:
    negative_testing_set = json.load(f) 
    
# set the observed value as 1 for the positive trainging set.    
for parepi in positive_training_set:
    parepi[2] = 1
for parepi in positive_testing_set:
    parepi[2] = 1
    

    
training_set = copy.deepcopy(positive_training_set)
training_set.extend(negative_training_set)

# Calculate the distance matrix and the total design matrix
distance_matrix = Distance_matrix(training_set, mode='Multiplication') # mode is an ajustable parameter
design_matrix = Design_matrix(distance_matrix, basis_function = 'Markov', radius = 1) # basis function is adjustable

# Generate the observed values
observed_values = np.zeros((len(training_set), 1))
for i in range(len(training_set)):
    observed_values[i,0]=training_set[i][2]

  
# Generate the testing design matrix
testing_set = copy.deepcopy(positive_testing_set)
testing_set.extend(negative_testing_set)

percentage = [1,0.8,0.4,0.2,0.1,0.05,0.02,0.01,0.008,0.005,0.002,0.001]
n_centers = []
for i in percentage:
    n_centers.append(math.floor(len(training_set)*i))
n_centers
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

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

new_RBFN_coverage_results = Coverage_centers_binary_prediction(\
             distance_matrix, design_matrix, training_set, testing_set, observed_values, n_centers)
new_RBFN_coverage_results['all_centers']
new_RBFN_coverage_results['n_centers'] 
new_RBFN_coverage_results['rate']
new_RBFN_coverage_results['all_parameters']
len(new_RBFN_coverage_results['all_parameters'])
# Save the results
with open('RBFN_coverage_results', 'w') as f:
    json.dump(new_RBFN_coverage_results, f,  cls=NumpyEncoder)
#
#with open('RBFN_coverage_results', 'r') as f:
#    res = json.load(f)

###########################################################################################

##################################################################################################
'''
Define a number controled remove
'''
def Coef_centers_binary_prediction(training_set, testing_set, design_matrix, observed_values):
    # Load the results from the RBFN_coverage_results'
    with open('RBFN_coverage_results', 'r') as f:
        RBFN_coverage_results = json.load(f)
    # Take out the number of centers, make sure the centers generated form the Coef 
    # method is the same as the coverage method.
    n_centers = RBFN_coverage_results['n_centers']    
    # Remove the duplicates in the training set, and use the unredundent dataset as
    # the centers.
    centers = Remove_duplicates(training_set)    
    #Initiate the parameters
    ncol = len(centers) + 1
    # Set the starting point
    coeff = np.zeros((ncol,1))
    parameter = {}
    parameter['coeff'] = coeff
    reg = np.ones((ncol,1))  
    parameter['reg'] = reg
    # Pay attention to the length of the centers after removing the duplicates
    # Creat a bunch of empty variables to store the results.    
    coef_all_parameter=[]
    coef_all_centers = []
    coef_all_rate = []
    coef_n_centers = []
    RBFN_coeff_results = {}
    for n in n_centers[1:]:
        if len(centers)>=n:
            while len(centers)>n:
                centers, parameter, removed_coeff = One_step_reduce_centers(centers, design_matrix, observed_values, parameter)
                print(removed_coeff, '   ', len(centers))        
        
            #Calculate testing_design_matrix
            testing_design_matrix = Generate_testing_design_matrix(testing_set, centers,\
                                                               training_set, mode = 'Multiplication', basis_function = 'Markov')
            # Binary prediction
            rate = Prediction(testing_set, parameter, testing_design_matrix)
            # Store the results. Pay attention to the mutable variables
            coef_all_rate.append(rate)
            # Deep copy the stuff
            stored_centers = copy.deepcopy(centers)
            stored_parameter = copy.deepcopy(parameter)
            # the centers
            coef_all_centers.append(stored_centers)
            # the parameters
            coef_all_parameter.append(stored_parameter)
            # the length of the centers
            coef_n_centers.append(len(centers))
    # Load the results into a dictionary, and return it.    
    RBFN_coeff_results['all_centers'] = coef_all_centers
    RBFN_coeff_results['all_parameters'] = coef_all_parameter
    RBFN_coeff_results['rate'] = coef_all_rate
    RBFN_coeff_results['n_centers'] = coef_n_centers
    RBFN_coeff_results['mode'] = 'Multiplication_Markov'
    
    return RBFN_coeff_results

########################################################################
New_RBFN_coeff_results =  Coef_centers_binary_prediction(training_set, testing_set, design_matrix, observed_values) 
New_RBFN_coeff_results['all_centers']
New_RBFN_coeff_results['all_parameters']
New_RBFN_coeff_results['rate'] 
New_RBFN_coeff_results['n_centers']
New_RBFN_coeff_results['mode']
######################################################

    
with open('RBFN_coeff_results', 'w') as f:
    json.dump(New_RBFN_coeff_results, f, cls=NumpyEncoder)    
    
with open('RBFN_coeff_results', 'r') as f:
    RBFN_coeff_results = json.load(f) 
RBFN_coeff_results.keys()
RBFN_coeff_results['rate']
RBFN_coeff_results['all_parameters']
RBFN_coeff_results['all_centers']
RBFN_coeff_results['top_x_percent_rate']
type(RBFN_coeff_results['all_parameters'][1]['coeff'])
##################################################################
def Plot_binary_prediction():
    # Load the results
    with open('RBFN_coeff_results', 'r') as f:
        RBFN_coeff_results = json.load(f) 
    with open('RBFN_coverage_results', 'r') as f:
        RBFN_coverage_results = json.load(f)
    
    
    x_coverage = RBFN_coverage_results['n_centers']
    y_coverage = RBFN_coverage_results['rate']
    x_coefficent = [x_coverage[0]]
    x_coefficent.extend(RBFN_coeff_results ['n_centers'])
    y_coefficient = [y_coverage[0]]
    y_coefficient.extend(RBFN_coeff_results ['rate'])
    
    
    
       
    log_percentage_coverage = []
    n = len(training_set)
    for i in x_coverage:
        log_percentage_coverage.append(-math.log(i/n))
        
    log_percentage_coefficient = []
    for i in x_coefficent:
        log_percentage_coefficient.append(-math.log(i/n))
    
    
    plt.figure(figsize = (8, 6))
    plt.plot(log_percentage_coverage, y_coverage, 'r--')
    plt.plot(log_percentage_coefficient, y_coefficient)
    plt.plot(log_percentage_coverage, y_coverage, 'go')
    plt.plot(log_percentage_coefficient, y_coefficient, 'bo')
    plt.legend(['Coverage Prunning', 'Coefficient Prunning'])
    
    plt.xlabel('-log(percentage)', fontsize = 20)
    plt.ylabel('Accuracy', fontsize = 20)
    plt.ylabel()
    plt.show()
    
    
    # Compare the RBFN results with the Logistic regression results
    os.chdir("/home/leo/Documents/Database/Pipeline/Ready_2_2_1_1")

    with open('Logistic_results', 'r') as f:
        result = json.load( f) 
    plt.figure(figsize = (8, 6))
    plt.plot(log_percentage_coverage, y_coverage, 'r--')
    plt.plot(log_percentage_coefficient, y_coefficient)
    
    percentage_logistic = result['truncate_percentage']
    log_percentage_logistic = []
    for i in percentage_logistic:
        log_percentage_logistic.append(-math.log(i))
#    plt.plot(log_percentage_logistic, result['Addition'], )
    plt.plot(log_percentage_logistic, result['Multiplication'])
    #plt.plot(log_percentage, result['Addition'], 'yo')
    #plt.plot(log_percentage, result['Multiplication'], 'co')
    plt.legend(['Coverage Prunning', 'Coefficient Prunning','Multiplication_distance'])
    plt.xlabel('-log(percentage)', fontsize = 20)
    plt.ylabel('Accuracy', fontsize = 20)
    plt.ylabel()
    plt.show()
    
    return

##################################################################

Plot_binary_prediction()

###########################################################################
'''
Faster method to calculate the testing design matrix. Too complicated, Let give up, for if the compution is not too overwhelming
,it is ok.
'''
def Pre_faster_testing_design_matrix(testing_set, training_set, positive_training_set, select_from = 'positive_training_set'):
    # Generate the prediction pool
    prediction_pool = Generate_prediction_pool(select_from=select_from, Ab_length = 2, Ag_length =2,\
                          positive_training_set=positive_training_set)
    # Generate the matrices between the testing and the training data
    testing_Ag_training_Ag_distance_matrix = np.zeros((len(testing_set), len(training_set)))
    testing_Ab_training_Ab_distance_matrix = np.zeros_like(testing_Ag_training_Ag_distance_matrix)
    l_Ag = len(training_set[0][1])
    l_Ab = len(training_set[0][0])
    for i in range(len(testing_set)):
        for j in range(len(training_set)):
            testing_Ag_training_Ag_distance_matrix[i,j] =  1 - (4*l_Ag + aligner.score(To_seq(testing_set[i][1]), To_seq(training_set[j][1])))/(15*l_Ag)
            testing_Ab_training_Ab_distance_matrix[i,j] =  1 - (4*l_Ab + aligner.score(To_seq(testing_set[i][0]), To_seq(training_set[j][0])))/(15*l_Ab)      

    # Generate the prediction_Ab_pool and the prediciton_Ag_pool 
    prediction_Ab_pool = []
    prediction_Ag_pool = []
    for i in prediction_pool:
        if i[0] not in prediction_Ab_pool:
            prediction_Ab_pool.append(i[0])
        if i[1] not in prediction_Ag_pool:
            prediction_Ag_pool.append(i[1])
            
    #  Generate the matrices between the pool and the training data
    pool_Ag_training_Ag_distance_matrix = np.zeros((len(prediction_Ag_pool), len(training_set)))
    pool_Ab_training_Ab_distance_matrix = np.zeros((len(prediction_Ab_pool), len(training_set)))
    for i in range(len(prediction_Ag_pool)):
        for j in range(len(training_set)):
            pool_Ag_training_Ag_distance_matrix[i,j]=1 - (4*l_Ag + aligner.score(To_seq(prediction_Ag_pool[i]), To_seq(training_set[j][1])))/(15*l_Ag)
            
    for i in range(len(prediction_Ab_pool)):
        for j in range(len(training_set)):
            pool_Ab_training_Ab_distance_matrix[i,j]=1 - (4*l_Ab + aligner.score(To_seq(prediction_Ab_pool[i]), To_seq(training_set[j][0])))/(15*l_Ab)  

    # Is it better to return the value in terms of dictionary
    many_matrics = {}
    many_matrics['testing_Ag_training_Ag_distance_matrix'] = testing_Ag_training_Ag_distance_matrix
    many_matrics['testing_Ab_training_Ab_distance_matrix'] = testing_Ab_training_Ab_distance_matrix
    many_matrics['pool_Ag_training_Ag_distance_matrix'] = pool_Ag_training_Ag_distance_matrix
    many_matrics['pool_Ab_training_Ab_distance_matrix'] = pool_Ab_training_Ab_distance_matrix
    # We have to generate the correct result here, so that we are able to tell whether our prediction is right or wrong.
    closest_Ab = []
    for i in range(len(testing_set)):
        testing_Ab = testing_set[i][0]
        Ab_score = -1000000
        closest = 0
        for j in range(len(prediction_Ab_pool)):
            pool_Ab = prediction_Ab_pool[j]
            score = aligner.score(To_seq(testing_Ab),To_seq(pool_Ab))
            if score > Ab_score:
                Ab_score = score
                closest = j
        closest_Ab.append(closest)
        
    closest_Ag = []
    for i in range(len(testing_set)):
        testing_Ag = testing_set[i][1]
        Ag_score = -1000000
        closest = 0
        for j in range(len(prediction_Ag_pool)):
            pool_Ag = prediction_Ag_pool[j]
            score = aligner.score(To_seq(testing_Ag),To_seq(pool_Ag))
            if score > Ag_score:
                Ag_score = score
                closest = j
        closest_Ag.append(closest)
        
    
    return many_matrics, closest_Ab, closest_Ag
            
            
'''
Define a function to find the complete design matrix when one testing sample is given. This design matrix with 
no row [1, 1, 1, ...].T attached

Input: 
    one_testing_sample:
        an integer, gives the position of the testing sample in the testing set, and this should be consistant with 
        the Ag_Ag or Ab_Ab matrix
'''
def One_complete_testing_design_matrix(one_testing_sample, many_matrices,\
                                       prediction_direction = 'Ag_to_Ab',basis_function = 'Markov', distance_mode = 'Multiplication'):
    # Unpack the matrices
    testing_Ag_training_Ag_distance_matrix = many_matrices['testing_Ag_training_Ag_distance_matrix']
    testing_Ab_training_Ab_distance_matrix = many_matrices['testing_Ab_training_Ab_distance_matrix']

    pool_Ag_training_Ag_distance_matrix = many_matrices['pool_Ag_training_Ag_distance_matrix']
    pool_Ab_training_Ab_distance_matrix = many_matrices['pool_Ab_training_Ab_distance_matrix']

    
    if prediction_direction == 'Ag_to_Ab':
        sliced_test_train= testing_Ag_training_Ag_distance_matrix [one_testing_sample, :]
        directioned_pool_train = pool_Ab_training_Ab_distance_matrix
    
    if prediction_direction == 'Ab_to_Ag':
        sliced_test_train= testing_Ab_training_Ab_distance_matrix [one_testing_sample, :]
        directioned_pool_train = pool_Ag_training_Ag_distance_matrix

    
    if basis_function == 'Markov':
        if distance_mode == 'Multiplication' :
            complete_testing_design_matrix = np.exp(- (directioned_pool_train * (sliced_test_train.reshape((1, -1)))))
        if distance_mode == 'Addition':
            complete_testing_design_matrix = np.exp(- (directioned_pool_train + (sliced_test_train.reshape((1, -1)))))
    if basis_function == 'Gaussian':
        if distance_mode == 'Multiplication':
            b_square = directioned_pool_train*directioned_pool_train
            a_square = sliced_test_train*sliced_test_train
            complete_testing_design_matrix = np.exp(- (b_square * (a_square.reshape((1, -1)))))
        if distance_mode == 'Addition':
            a_plus_b = directioned_pool_train + (sliced_test_train.reshape((1, -1)))
            complete_testing_design_matrix = np.exp(- (a_plus_b*a_plus_b))
    
    # Since we will use the linear conbination of this design matrix, it is reasonable to normalize it
    average = np.average(complete_testing_design_matrix, axis= 0)
    complete_testing_design_matrix = complete_testing_design_matrix - average
    std = np.std(complete_testing_design_matrix, axis = 0)
    complete_testing_design_matrix = complete_testing_design_matrix / std
            
    n_row = complete_testing_design_matrix.shape[0]
    complete_testing_design_matrix = np.hstack((complete_testing_design_matrix, np.ones((n_row, 1))))
        
    return complete_testing_design_matrix

###############################################################################################
# Lets test whether the above method works by comparing with the ordinary method           
'''
Define an ordinary method
'''            
def Oridnary_method_testing_design_matrix(combination, training_set, basis_function = 'Markov',
                                          distance_mode = 'Multiplication'):
    n_row = len(combination)
    n_col = len(training_set)
    # calculate the distance matrix
    distance_matrix = np.zeros((n_row, n_col))
    for i in range(n_row):
        for j in range(n_col):
            if distance_mode  == 'Addition':
                distance_matrix[i,j] = Addition_distance(To_seq(combination[i][0]),To_seq(training_set[j][0]),
                               To_seq(combination[i][1]),To_seq(training_set[j][1]))

            if distance_mode  == 'Multiplication':
                distance_matrix[i,j] = Multiplication_distance(To_seq(combination[i][0]),To_seq(training_set[j][0]),
                               To_seq(combination[i][1]),To_seq(training_set[j][1]))
                
    complete_testing_design_matrix = Design_matrix(distance_matrix, basis_function = 'Markov', radius = 1) # basis function is adjustable    
    
    return complete_testing_design_matrix



################################################################


#################################################################################################
# Do the concrete prediction, take the rank/percentile into consideration


'''
Input:
    positive_testing_position:
        An integer, gives the position of the testing sample in the positive testing set.
    RBFN_results:
        there are two types of RBFN_results, coeff and the coverage, for each type,
        the function will return corresponding correct rate
    complete_testing_design_matrix:
        a testing design matric corresponding to the testing sample with index given by 
        positive_testing_position
    closest:
        a list with length the same as the length of the all_centers i the RBFN_results.
        gives the correct position of the Ab of Ag to be predicted
Output:
    pred:
        a list with the same length as the closest. The elements in it take values of
        either 0 or 1. With 1 means a correct prediction and 0 a wrong prediction.
          
    
'''    


def Concrete_pred(RBFN_results, complete_testing_design_matrix, closest_index):

    all_centers = RBFN_results['all_centers']
    concrete_pred = []
    for i in range(len(all_centers)):
        
        centers = all_centers[i]    
        parameter = RBFN_results['all_parameters'][i]
        
        centers.append(-1)
        testing_design_matrix = complete_testing_design_matrix[:, centers]
        centers.remove(-1)
        
        res = Top_x_percent_prediction(testing_design_matrix, parameter, closest_index, top_x_percent = 0.1)
        
        concrete_pred.append(res)
    
    return concrete_pred
###########################################################################
'''
define a function to do all the concrete_pred for different types of RBFN results

Inputs:
    testing_set:
        As above
    positive_training_set:
        As above
    prediction_parameter:
        a dictionary, contains
        prediction_parameter['select_from']
        prediction_parameter['prediction_direction']
        prediction_parameter['basis_function']
        prediction_parameter['distance_mode']
Outputs:
    RBFN_coeff_results:
        a dictionary, with one more key called 'correct_rate', which is a list 
        with the same length as all_centers, and the elements give the corresponding
        correct rate
    RBFN_coverage_results :
        The same as above.    
'''

def Packed_top_x_percent_prediction(testing_set, training_set, positive_training_set, prediction_parameter):
    # Load the results from the above calculations
    with open('RBFN_coeff_results', 'r') as f:
        RBFN_coeff_results = json.load(f) 
    with open('RBFN_coverage_results', 'r') as f:
        RBFN_coverage_results = json.load(f)
    # Calculate the many matrices used for calculating the testing design matrices   
    many_matrices, closest_Ab, closest_Ag = Pre_faster_testing_design_matrix(testing_set,\
                                                    training_set, positive_training_set,\
                                                     select_from= prediction_parameter['select_from'])
    # Creat empty nd arrays to store the results
    coeff_all_centers = RBFN_coeff_results['all_centers']
    coverage_all_centers = RBFN_coverage_results['all_centers']
    pred_coeff =  np.zeros((1, len(coeff_all_centers)))
    pred_coverage = np.zeros((1, len(coverage_all_centers)))
    # Do the prediction
    for i in range(len(testing_set)):    
        # Calculate the complete testing design matrix
        complete_testing_design_matrix = One_complete_testing_design_matrix(i,\
                                    many_matrices, prediction_direction = prediction_parameter['prediction_direction'],\
                                       basis_function = prediction_parameter['basis_function'],\
                                       distance_mode = prediction_parameter['distance_mode'])
        
        # Do the prediction for the one testing sample under different centers
        if prediction_parameter['prediction_direction'] == 'Ag_to_Ab':
            one_test_closest = closest_Ab[i]
        elif prediction_parameter['prediction_direction'] == 'Ab_to_Ag':
            one_test_closest = closest_Ag[i]
        # 
        coeff_one_pred = Concrete_pred(RBFN_coeff_results,complete_testing_design_matrix,one_test_closest)
        coverage_one_pred = Concrete_pred(RBFN_coverage_results,complete_testing_design_matrix,one_test_closest)
        # Change the above into nd array and add to pred_coeff and prec_coverage
#        print(coeff_one_pred)
        pred_coeff += np.array(coeff_one_pred)
        pred_coverage += np.array(coverage_one_pred)
        # Keep the row number of the complete_testing_design_matrix to calculate the correct rate
        n_pool = complete_testing_design_matrix.shape[0]
    
    # Calculate the rate
    pred_coeff = np.round(pred_coeff/n_pool, 3)
    pred_coverage = np.round(pred_coverage/n_pool, 3)
    
    # Load the results to the dictionaries
    RBFN_coeff_results['top_x_percent_rate'] = pred_coeff
    RBFN_coverage_results['top_x_percent_rate'] = pred_coverage
    
    return RBFN_coeff_results, RBFN_coverage_results

# Creat a prediction parameter to feed the function
prediction_parameter = {}
prediction_parameter['select_from'] = 'positive_training_set'
prediction_parameter['prediction_direction'] = 'Ag_to_Ab'
prediction_parameter['basis_function'] = 'Markov'
prediction_parameter['distance_mode'] = 'Multiplication'
RBFN_coeff_results, RBFN_coverage_results = Packed_top_x_percent_prediction(testing_set,\
                                            training_set, positive_training_set, prediction_parameter)     
RBFN_coeff_results ['top_x_percent_rate']
RBFN_coverage_results['top_x_percent_rate']       
RBFN_coeff_results ['n_centers']

# Plot the result
def Concrete_plot():
    x_coverage = RBFN_coverage_results['n_centers']
    y_coverage = RBFN_coverage_results['top_x_percent_rate'][0]
    x_coefficent = [x_coverage[0]]
    x_coefficent.extend(RBFN_coeff_results ['n_centers'])
    y_coefficient = [y_coverage[0]]
    y_coefficient.extend(RBFN_coeff_results ['top_x_percent_rate'][0])
    
    len(x_coefficent)
    len(y_coefficient)
    
       
    log_percentage_coverage = []
    n = len(training_set)
    for i in x_coverage:
        log_percentage_coverage.append(-math.log(i/n))
        
    log_percentage_coefficient = []
    for i in x_coefficent:
        log_percentage_coefficient.append(-math.log(i/n))
    
    
    plt.figure(figsize = (8, 6))
    plt.plot(log_percentage_coverage, y_coverage, 'r--')
    plt.plot(log_percentage_coefficient, y_coefficient)
    plt.plot(log_percentage_coverage, y_coverage, 'go')
    plt.plot(log_percentage_coefficient, y_coefficient, 'bo')
    plt.legend(['Coverage Prunning', 'Coefficient Prunning'])
    
    plt.xlabel('-log(percentage)', fontsize = 20)
    plt.ylabel('Accuracy', fontsize = 20)
    plt.ylabel()
    plt.show()    
    return
#######################################################################
'''
What if we train according to the centers and predict
'''
def Train_coeff_centers_prediction(positive_testing_set, training_set, positive_training_set, prediction_parameter):
    # Load the results from the above calculations
    with open('RBFN_coeff_results', 'r') as f:
        RBFN_coeff_results = json.load(f) 
    with open('RBFN_coverage_results', 'r') as f:
        RBFN_coverage_results = json.load(f)
    # Calculate the many matrices used for calculating the testing design matrices   
    many_matrices, closest_Ab, closest_Ag = Pre_faster_testing_design_matrix(positive_testing_set,\
                                                    training_set, positive_training_set,\
                                                     select_from= prediction_parameter['select_from'])
    # Creat empty nd arrays to store the results
    coeff_all_centers = RBFN_coeff_results['all_centers']
    coverage_all_centers = RBFN_coverage_results['all_centers']
    pred_coeff =  np.zeros((1, len(coeff_all_centers)))
    pred_coverage = np.zeros((1, len(coverage_all_centers)))
    RBFN_coeff_results['all_parameters'] = []
    for i in range(len(coeff_all_centers)):
        centers = coeff_all_centers[i]
        # Calculate the design_matrix
        centers.append(-1)
        design_matrix1 = design_matrix[:, centers]
    #    design_matrix1.shape
        centers.remove(centers[-1])   
        # Train
        parameter, loss = Train_RBFN_BFGS(design_matrix1, observed_values, rho=0.9, c = 1e-3, termination = len(centers)/1000,\
                            parameter_inheritance = False, parameter=None) 
        # load the centers
        RBFN_coeff_results['all_parameters'].append(copy.deepcopy(parameter))
        
    # Do the prediction
    for i in range(len(positive_testing_set)):    
        # Calculate the complete testing design matrix
        complete_testing_design_matrix = One_complete_testing_design_matrix(i,\
                                    many_matrices, prediction_direction = prediction_parameter['prediction_direction'],\
                                       basis_function = prediction_parameter['basis_function'],\
                                       distance_mode = prediction_parameter['distance_mode'])
        
        # Do the prediction for the one testing sample under different centers
        if prediction_parameter['prediction_direction'] == 'Ag_to_Ab':
            one_test_closest = closest_Ab[i]
        elif prediction_parameter['prediction_direction'] == 'Ab_to_Ag':
            one_test_closest = closest_Ag[i]
        # 
        coeff_one_pred = Concrete_pred(RBFN_coeff_results,complete_testing_design_matrix,one_test_closest)
        coverage_one_pred = Concrete_pred(RBFN_coverage_results,complete_testing_design_matrix,one_test_closest)
        # Change the above into nd array and add to pred_coeff and prec_coverage
#        print(coeff_one_pred)
        pred_coeff += np.array(coeff_one_pred)
        pred_coverage += np.array(coverage_one_pred)
        # Keep the row number of the complete_testing_design_matrix to calculate the correct rate
        n_pool = complete_testing_design_matrix.shape[0]
    
    # Calculate the rate
    pred_coeff = np.round(pred_coeff/n_pool, 3)
    pred_coverage = np.round(pred_coverage/n_pool, 3)
    
    # Load the results to the dictionaries
    RBFN_coeff_results['top_x_percent_rate'] = pred_coeff
    RBFN_coverage_results['top_x_percent_rate'] = pred_coverage
    
    return RBFN_coeff_results, RBFN_coverage_results

RBFN_coeff_results, RBFN_coverage_results = Train_coeff_centers_prediction(positive_testing_set,\
                                            training_set, positive_training_set, prediction_parameter)
RBFN_coeff_results['top_x_percent_rate']
RBFN_coverage_results['top_x_percent_rate']
########################################################################
'''
The slow prediction method
Input:
    design_matrix:
        it is the big design matrix given by training_set by training_set
    positive_training_set:
        The same as above.
    observed_values:
        The same as avove
Output:
    rate:
        a list corresponding to the set of centers with the first few center set excluded.
'''
def Slow_concrete_pred(design_matrix, positive_training_set, observed_values):
    with open('RBFN_coeff_results', 'w') as f:
        json.dump(RBFN_coeff_results,f,  cls=NumpyEncoder)
    with open('RBFN_coverage_results', 'w') as f:
        json.dump(RBFN_coverage_results,f,  cls=NumpyEncoder)
    # Concrete prediction, the slow method
    # Generate the prediction pool
    prediction_pool = Generate_prediction_pool(select_from='positive_training_set', Ab_length=2, Ag_length=2,\
                                  positive_training_set=positive_training_set)
    rate = []
    for centers in RBFN_coeff_results['all_centers'][3:]:
            # Calculate the design_matrix
        centers.append(-1)
        design_matrix1 = design_matrix[:, centers]
    #    design_matrix1.shape
        centers.remove(centers[-1])   
        # Train
        parameter, loss = Train_RBFN_BFGS(design_matrix1, observed_values, rho=0.9, c = 1e-3, termination = len(centers)/1000,\
                            parameter_inheritance = False, parameter=None)
        
        n = 0
        for i in range(len(positive_testing_set)):
            print(i)
            one_test = positive_testing_set[i]
            combination, closest_index = Pre_concrete_prediction(one_test,\
                                                                 prediction_pool, prediction_direction= 'Ag_to_Ab')
     
            #Calculate testing_design_matrix
            testing_design_matrix = Generate_testing_design_matrix(combination, centers,\
                                                               training_set, mode = 'Multiplication', basis_function = 'Markov')
            pred = testing_design_matrix.dot(parameter['coeff'])
            pred = pred.flatten().tolist()
            closest_score = pred[closest_index]
            pred.sort(reverse=True)
            n_top = math.floor(len(combination)*0.1)
            n_top_score = pred[n_top]
            if n_top_score <= closest_score:
                n += 1
        rate.append(round(n/len(combination), 3))
    return rate
rate = Slow_concrete_pred(design_matrix, positive_training_set, observed_values)  
    
rate











