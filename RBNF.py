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
def Generate_random_negative(sample_size, data, Ab_lenth, Ag_length):
    
    TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
                    ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
                    ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]
    AA = []
    for aa in TripleSingle:
        AA.append(aa[0])
        
    Ab_Ag = []
    for parepi in data:
        Ab_Ag.append([parepi[0], parepi[1]])
        
    negative_samples = []
    while len(negative_samples) < sample_size:
        r_Ab_r_Ag = []
        while r_Ab_r_Ag == []:        
            r_Ab = random.sample(AA, Ab_lenth)
            r_Ag = random.sample(AA, Ag_length)
            r_Ab_r_Ag  = [r_Ab, r_Ag]
            if r_Ab_r_Ag in Ab_Ag:
                r_Ab_r_Ag  = []
        negative_samples.append([r_Ab, r_Ag, -1, -1])
    return negative_samples
#####################################################################################
    
os.chdir("/home/leo/Documents/Database/Pipeline/Ready_2_2_1_1")
with open('ready_2_2_1_1__cn_gate_1_all', 'r') as f:
    data = json.load(f)

negative_samples = Generate_random_negative(len(data), data, 2, 2)
len(negative_samples)
negative_samples[:6]
'''
change the value of the positive samples to 1
'''
def Processing_positive(data):
    positive_data = []
    for parepi in data:
        positive_data.append([parepi[0], parepi[1], 1, 1])
    return positive_data
        
positive_data = Processing_positive(data)        
len(positive_data)        
positive_data[:6]  
      
def Generate_testing_set(data):
    n = math.floor(len(data) * 0.1)
    testing_set = random.sample(data, n)
    for i in testing_set:
        data.remove(i)
    return data, testing_set
positive_training_set, positive_testing_set = Generate_testing_set(positive_data)
negative_training_set, negative_testing_set = Generate_testing_set(negative_samples)
#len(positive_training_set)
#positive_training_set[:6]
#len(positive_testing_set)
#positive_testing_set[:6]
#len(negative_testing_set)
#negative_training_set[:6]
#len(negative_testing_set)
#negative_testing_set[:6]

training_set = copy.deepcopy(positive_training_set)
training_set.extend(negative_training_set)
len(training_set)

distance_matrix = Distance_matrix(training_set, mode='Multiplication')
distance_matrix.shape
distance_matrix[:5, :5]
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

def Generate_centers(distance_matrix, cover_radius):
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
    samples_to_centers_matrix = distance_matrix[:, centers]
    return centers, samples_to_centers_matrix
#####################################################################################
    
n_samples = distance_matrix.shape[0]
cover = np.arange(0.1, 5, 0.1)
cover_radius = cover[0]
centers, samples_to_centers_matrix = Generate_centers(distance_matrix, cover_radius)
len(centers)
n_centers = math.floor(n_samples * 0.1)
i = 0
while len(centers) > n_centers and i < len(cover)-1:
    i += 1
    print(len(centers), '   ', i)
    cover_radius = cover[i]
    centers, samples_to_centers_matrix = Generate_centers(distance_matrix, cover_radius)
    
len(centers)
centers    
for i in centers:
    for j in centers:
        if i != j and distance_matrix[i,j] <= cover_radius:
            print('Something is wrong'\
                  )
            break  
samples_to_centers_matrix.shape
#import numpy as np
#matrix = np.zeros((6, 6))
#for i in range(6):
#    for j in range(6):
#        if i != j:
#            matrix[i, j] = 1/abs(i-j)
#        else:
#            matrix[i,j] = 0
#centers, samples_to_centers_matrix = Generate_centers(matrix, cover_radius=0.3)
#centers        
#samples_to_centers_matrix.round(3)
#matrix.round(3)
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
    std = np.std(design_matrix)
    design_matrix = design_matrix / std
    # Should we add a constant column?  Lets add.
    design_matrix = np.hstack((design_matrix, np.ones((nrow,1))))
           
    return design_matrix
#####################################################################################

design_matrix = Design_matrix(samples_to_centers_matrix, basis_function = 'Markov', radius = 1)    
design_matrix.shape
design_matrix[:5, :5]

design_matrix = Design_matrix(distance_matrix, basis_function = 'Markov', radius = 1)
design_matrix.shape
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
    design_matrix:
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
observed_values = np.zeros((len(training_set), 1))
for i in range(len(training_set)):
    observed_values[i,0]=training_set[i][2]
len(observed_values)
observed_values[:5]    
training_set[:5]
#ncol = np.shape(design_matrix)[1]
#coeff = np.zeros((ncol,1))
#parameter = {}
#parameter['coeff'] = coeff
##The reg should not be negative. It is better that reg > delta, a small positive number
#reg = np.ones((ncol,1)) * 0.1
#parameter['reg'] = reg
#loss, gradient = Loss(design_matrix, observed_values, parameter)

###################################################################################

#import numpy as np
#a = np.arange(0, 6)
#a = a.reshape((6,1))
#a.shape
#c = np.arange(0, 16)
#c = c.reshape((-1, 1))
#c[:8, 0].reshape((-1, 1))
#c[8:16,0]
#np.hstack((a, np.ones((2, 1))))
#b = np.ones((6,1))*2
#a*b
#a*a
#a
#np.where(a<3, 3, a)
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
        print(loss, '    ', grad_square)
    return parameter, loss
###########################################################################################################
parameter, loss = Train_RBFN_BFGS(design_matrix, observed_values, rho=0.85, c = 1e-3, termination = 1e-2)
parameter['coeff'].shape
###############################################################################
'''
Generating the testing design matrix
'''
len(training_set)
def Generate_testing_design_matrix(testing_data, centers, training_set, mode = 'Multiplication', basis_function = 'Markov'):
    centers_parepi = [training_set[i] for i in centers]
    nrow = len(testing_data)
    ncol = len(centers_parepi)
    distance_matrix = np.zeros((nrow,ncol))
    for i in range(nrow):
        for j in range(ncol):
            if mode == 'Addition':
                distance_matrix[i,j] = Addition_distance(To_seq(testing_data[i][0]),To_seq(centers_parepi[j][0]),
                               To_seq(testing_data[i][1]),To_seq(centers_parepi[j][1]))
            if mode == 'Multiplication':
                distance_matrix[i,j] = Multiplication_distance(To_seq(testing_data[i][0]),To_seq(centers_parepi[j][0]),
                               To_seq(testing_data[i][1]),To_seq(centers_parepi[j][1]))
                
    testing_design_matrix = Design_matrix(distance_matrix, basis_function, radius = 1)
    
    return testing_design_matrix

testing_data = copy.deepcopy(positive_testing_set)
testing_data.extend(negative_testing_set)
len(testing_data)
testing_data[:6]
testing_data[-6:-1]
testing_design_matrix = Generate_testing_design_matrix(testing_data, centers,\
                                                       training_set, mode = 'Multiplication', basis_function = 'Markov')
testing_design_matrix.shape
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

rate = Prediction(testing_data, parameter, testing_design_matrix)

#################################################################################
'''
Do a concrete predict. Don't just predict whether they contact or not, but also the 
exact paratope
'''
len(positive_testing_set)
len(negative_testing_set)
len(testing_data)
def Concrete_prediction(one_positive_sample, centers, positive_training_set, parameter, \
                        training_set, mode = 'Multiplication', basis_function = 'Markov'):
    # get the positive matches with the one_positive_sample, and the true closest.
    closest = 0
    score = -10000
    positive_parepi = []
    for ind in range(len(centers)):
        i = centers[ind]
        if i <= len(positive_training_set)-1:
            positive_parepi.append([positive_training_set[i][0], one_positive_sample[1]])
            align_score = aligner.score(To_seq(one_positive_sample[0]), To_seq(positive_training_set[i][0]))
            if align_score > score:
                score = align_score
                closest = ind
            
    #Calculate the testing_design matrix
    testing_design_matrix = Generate_testing_design_matrix(positive_parepi, centers, training_set, \
                                   mode, basis_function)
    # Prediction
    coeff = parameter['coeff']
    prediction = testing_design_matrix.dot(coeff)
    prediction = prediction.flatten().tolist()
    # find the predictin score for the closest
    score_for_closest = prediction[closest]
    # find the lowest prediction score for the top 10 percent
    n_top = math.floor(len(positive_parepi) * 0.1)
    prediction.sort(reverse = True)
    n_top_score = prediction[n_top]
    # Tell whether our prediction is correct
    if score_for_closest >= n_top_score:
        return True
    else:
        return False
################################################################################################3
'''
concrete_prediction_rate:
    A function gives the concrete prediction rate, which can be compared with the hierachical clustering
'''
def Concrete_prediction_rate(positive_testing_set, centers, positive_training_set, training_set,\
                             mode = 'Multiplication', basis_function = 'Markov'):
    n = 0
    j = 0
    for i in positive_testing_set:
        j += 1
        print(j)
        pred = Concrete_prediction(i, centers, positive_training_set, parameter, training_set,\
                                   mode, basis_function)
        if pred:
            n += 1
    return round(n/len(positive_testing_set), 2)        
#a = np.arange(0, 8)
#a = a.reshape((-1, 1))
#a
#a = a.flatten().tolist()
#a
#a.sort(reverse=True)
#a
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
    ratio = 3000/len(centers)
    termination = 10**(-3)

        
    parameter, loss = Train_RBFN_BFGS(new_design_matrix, observed_values, rho=0.9, c = 1e-3, termination=termination,\
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
    

import numpy as np
a = np.arange(-3, 5)
a = a.reshape((-1, 1))
a.flatten().tolist()
a
b = [1, 3]
select = copy.deepcopy(b)
select.append(-1)
c = a[:, select]
c
b
a
select
np.delete(a, [1])
d = [-1, 0, 2]
np.abs(a)
a
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

centers = Remove_duplicates(training_set)

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
for i in range(1):   
    centers, parameter, cutoff_coeff = Coeff_select_centers(parameter, centers, training_set, design_matrix,\
                                                            observed_values, control_coeff=1e-1, centers_inherit = True)
centers, parameter, cutoff_coeff = Coeff_select_centers(parameter, centers, training_set, design_matrix,\
                                                        observed_values, control_coeff=1e-1, centers_inherit = False)
len(centers)
centers
testing_data = copy.deepcopy(positive_testing_set)
testing_data.extend(negative_testing_set)
len(testing_data)
testing_data[:6]
testing_data[-6:-1]
testing_design_matrix = Generate_testing_design_matrix(testing_data, centers,\
                                                       training_set, mode = 'Multiplication', basis_function = 'Markov')
testing_design_matrix.shape
rate = Prediction(testing_data, parameter, testing_design_matrix)    
rate        
Concrete_prediction_rate(positive_testing_set, centers, positive_training_set, training_set,\
                             mode = 'Multiplication', basis_function = 'Markov')


len(positive_training_set)
n = 0
for i in centers:
    if i <= len(positive_training_set)-1:
        n +=1
n

