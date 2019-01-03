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


def Train_RBFN_BFGS(design_matrix, observed_values, rho=0.9, c = 1e-3, termination = 1e-2):
    nrow = np.shape(design_matrix)[0]
    ncol = np.shape(design_matrix)[1]
    
    # Give the initial Hessian H. The variables are the coeff and reg
    H = np.eye(ncol)/(10*nrow)
    # Set the starting point
    coeff = np.zeros((ncol,1))
    parameter = {}
    parameter['coeff'] = coeff
    #The reg should not be negative. It is better that reg > delta, a small positive number
    reg = np.ones((ncol,1)) * 1
    parameter['reg'] = reg
    # Stack up the variables
#    stacked_coeff = np.vstack((coeff, reg))
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
    
n = 0
j = 0
for i in positive_testing_set:
    j += 1
    print(j)
    pred = Concrete_prediction(i, centers, positive_training_set, parameter, training_set,\
                               mode = 'Multiplication', basis_function = 'Markov')
    if pred:
        n += 1
round(n/len(positive_testing_set), 2)        
#a = np.arange(0, 8)
#a = a.reshape((-1, 1))
#a
#a = a.flatten().tolist()
#a
#a.sort(reverse=True)
#a
################################################################################

class RBFN(object):
    def __init__(self, training, testing, basis_function = 'Markov', distance_mode = 'Multiplication_distance'):
        self.training = training
        self.testing = testing
        self.distance_mode = distance_mode
        self.basis_function = basis_function
        self.params = {}
        observed_contact = []
        for parepi in training:
            observed_contact.append(parepi[2])
        self.observed_contact = np.asanyarray(observed_contact)
        
    def Centers_from_NewDataAnalysis(self, center_mode='positive_only'):
        if center_mode == 'positive_only':
            with open('RBFN_centers_Ab_in_Ag_Mouse', 'r') as f:
                centers = json.load(f)
            self.centers = centers
            self.params['W'] =  0.1 * np.random.randn(len(self.centers), 1)
        elif center_mode == 'augmented':
            with open('Augmented_centers_Mouse', 'r') as f:
                centers = json.load(f)
            self.centers = centers
            self.params['W'] =  0.1 * np.random.randn(len(self.centers), 1)
            


    
    
    '''
    Pruning:
        to return the prunned centers
    Input: 
          distance_matrix: the distance matrix of all the training samples
          cover_radius: gives the radius a center can cover
    Output: 
        centers:
            in the same form as the training set
        Samples_to_center_matrix:
            A matrix gives the distance between the samples and the centers, it 
            can be used directly to calculate the design matrix            
    '''

        
        
        


        
        
        

    def Formula_ridge_find_reg(self):
        basis_function = self.basis_function
        distance_mode = self.distance_mode
        training = self.training
        centers = self.centers
        
        REG_vector = [0.1, 1, 5, 10, 50, 60, 70, 80, 90, 100, 150]
        Radi = [0.01, 0.04, 0.16, 0.64, 0.81, 1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 144, 225]
#        GCV_vector = []
        GCV_lowest = 1000000
        reg_best = 0
        
        for radius in Radi:
            
            design_matrix = Design_matrix(training, centers, basis_function, radius, distance_mode)
            
            for reg in REG_vector:
            
                A = (design_matrix.T).dot(design_matrix) + np.eye(len(centers))*reg
                
                inverse_A = np.linalg.inv(A)
                
                W = inverse_A.dot((design_matrix.T)).dot(self.observed_contact)
                
                P = np.eye(len(training)) - design_matrix.dot(inverse_A).dot(design_matrix.T)
                
                GCV = len(self.training) * ((self.observed_contact).T).dot(P.dot(P)).dot(self.observed_contact)
                
                GCV /= (np.trace(P))**2
                
#                GCV_vector.append(GCV)
                
                if GCV <= GCV_lowest:
                    GCV_lowest = GCV
                    reg_best = reg
                    self.radius = radius
                    self.params['W'] = W
        
        self.GCV_lowest = GCV_lowest
        self.reg_best = reg_best
            
#        plt.title('GCV Vs cut REG')
#        plt.xlabel('REG')
#        plt.ylabel('GCV')
#        plt.plot(REG_vector, GCV_vector)
#        plt.show()
#        plt.close() 
        
    '''
     predictions:
         are given in the form of [[[1, 2], [6, 8]], ...]
         1, 2 are the predicted centers that the Ab belong to, 6, 8 are the predicted number
         of contact.
    '''
    def Formula_ridge_predict_Ab(self):
        W = self.params['W']
        testing = self.testing
        centers = self.centers
        distance_mode = self.distance_mode
        basis_function = self.basis_function
        radius = self.radius
        positive_centers = []
        for i in centers:
            if i[2] > 0:
                positive_centers.append(i)
        self.positive_centers = positive_centers
        top_n = math.floor(len(positive_centers) * 0.1)
          
        correct_prediction = 0
        for test_Ag in testing:
            # do prediction for test_Ag
            # get the combination and the distance scores between the observed and the centers
            test_Ag_positive_centers_scores = []
            test_Ag_positive_centers = []
            for i in range(len(positive_centers)) :
                # conbine the test_Ag with different positive centers Ab, and the predicted contact 
                # of the combination will be caltulated later
                test_Ag_positive_centers.append([positive_centers[i][0], test_Ag[1]])
                # calculate the distance between the observed Ab and the Ab of the centers
                l_Ab = len(test_Ag[0])
                Ab_seq1 = To_seq(test_Ag[0])
                Ab_seq2 = To_seq(positive_centers[i][0])
                distance = (1 - (4*l_Ab + aligner.score(Ab_seq1, Ab_seq2))/(15*l_Ab))
                test_Ag_positive_centers_scores.append([i, distance])
            # sort the test_Ag_positive_centers_scores, and chose the nearest center index
            test_Ag_positive_centers_scores.sort(key = lambda x:x[1])
            predicted_center = test_Ag_positive_centers_scores[0][0]
                
            design_matrix = Design_matrix(test_Ag_positive_centers, centers,
                                          basis_function, radius, distance_mode)
            # calculate the scores
            scores = design_matrix.dot(W.T)
            # indexed scores
            indexed_scores = []
            for i in range(len(scores)):
                indexed_scores.append([i, scores[i]])
            # sort the indexed scores to find the top n
            indexed_scores.sort(key = lambda x:x[1])
            #find the top_n prediction
            top_n_predictions = []
            for i in range(top_n):
                top_n_predictions.append(indexed_scores[i][0])
            # Check if the prediction is correct
            if predicted_center in top_n_predictions:
                correct_prediction += 1
                
        self.correct_rate =round(correct_prediction/len(testing), 2)
        
    def Add_negative_samples(self):
        with open('negative_samples_RBFN_Mouse', 'r') as f:
            negative_samples = json.load(f)

        self.training.extend(negative_samples)
        
    def Feedforward_selection(self):
        # pick the one with the largest frequency as the starting vectors
        centers = []
        training = copy.deepcopy(self.training)
        sliced_training = []
        for parepi in training:
            sliced_training.append(parepi[:2])
        container = []
        for parepi in sliced_training:
            frequency = 0 
            for i in sliced_training:
                if i == parepi:
                    frequency += 1
            container.append([parepi, frequency])
        container.sort(key= lambda x:x[1], reverse = True)
        n = container[0][1]
        for parepi in container:
            if parepi[1] == n and [parepi[0][0], parepi[0][1], 0, 0] not in centers:
                centers.append([parepi[0][0], parepi[0][1], 0, 0])
#        centers.append(training[0])

                
        # calculate A, P, W, reg, C for the centers selected above
        reg = [0]; C = [0]
        
        basis_function = self.basis_function
        distance_mode = self.distance_mode
        training = self.training
        radius = self.radius
        observed_contact = self.observed_contact
        
        design_matrix = Design_matrix(training, centers, basis_function, radius, distance_mode)
        
        A = (design_matrix.T).dot(design_matrix) + np.eye(len(centers))*reg[-1]

        inverse_A = np.linalg.inv(A)
        
        W = inverse_A.dot((design_matrix.T)).dot(observed_contact)
        
        P = np.eye(len(training)) - design_matrix.dot(inverse_A).dot(design_matrix.T)
        
        C_update =  len(training) * ((observed_contact).T).dot(P.dot(P)).dot(observed_contact)
                
        C_update  /= (np.trace(P))**2
        
        C.append(C_update)
        
        reg_update = (observed_contact.T).dot(P.dot(P.T)).dot(observed_contact)*(np.trace(inverse_A))
        reg_update /= (W.T).dot(inverse_A).dot(W)*(np.trace(P))
        reg[0] = reg_update
        
        # do the iteration
        improve_all = [C_update]
        for i in range(len(training)):
            for parepi in training:
                largest_improve = -1000
                if parepi not in centers:
                    # new column in the design matrix
                    hm = Design_matrix(training, [parepi], basis_function, radius, distance_mode)
                    #update P first
                    delta = reg[-1] + (hm.T).dot(P).dot(hm)
                    temp_P = P - P.dot(hm).dot(hm.T).dot(P)/delta
    
                    # update Cost, find the improvement
                    C_update =  len(training) * ((observed_contact).T).dot(temp_P.dot(temp_P.T)).dot(observed_contact)
                    C_update  /= (np.trace(temp_P))**2
                    improve = C[-1] - C_update
                    if improve > largest_improve:
                        largest_improve = improve
                        appended_center = parepi
                        best_hm = hm
                        best_P = temp_P
                        best_delta = delta
                        
            # give a conditon to abort the loop 
#            if largest_improve <= 0:
#                break
            # update the chosen center and related parameters
            centers.append(appended_center)
            # append the improve
            improve_all.append(largest_improve)
            # delta
            delta = best_delta
            #Calculate a beta
            beta = inverse_A.dot(design_matrix.T).dot(hm)
            beta = np.vstack((beta, -1*np.ones((1 ,1))))
            # update design matrix
            hm = best_hm
            design_matrix =  np.hstack((design_matrix, hm))
            # updata P
            P = best_P
            # update the inverse_A
            # Enlarge the dimensions by 1
            (row_n, col_n) = np.shape(inverse_A)
            empty_col = np.zeros((row_n, 1))
            empty_row = np.zeros((1, col_n+1))
            inverse_A = np.hstack((inverse_A, empty_col))
            inverse_A = np.vstack((inverse_A, empty_row))
            # do the update 
            inverse_A += beta.dot(beta.T)/delta  
            # update reg
            W = inverse_A.dot(design_matrix.T).dot(observed_contact)
            denominator = (W.T).dot(inverse_A).dot(W)*np.trace(P)
            numerator = (observed_contact.T).dot(P).dot(P.T).dot(observed_contact) 
            multiplier = np.trace(inverse_A - reg[-1]* (inverse_A.dot(inverse_A.T)))
            numerator *= multiplier
            reg.append( numerator/denominator)

            
        #load the results
            self.forward_centers = centers
            self.forward_improve_all = improve_all
            self.forward_reg = reg
            self.forward_design_matrix = design_matrix
            self.forward_W = W
            
            if len(centers) == math.floor(len(training) * 0.05):
                break
                
                
                
      
            
                
        
        
               
            
            
  
        
#        for parepi in testing:
#            # find the top n predictions
               
        

#    def Loss(self, reg = 1):
#        W = self.params['W']
#        b = self.params['b']
#        distance_function = self.distance_function
#        training = self.training
#        centers = self.centers
#        design_matrix = np.zeros((len(training), len(centers)))
#        for (i, j), value in np.ndenumerate(design_matrix):
#            design_matrix[i, j] = distance_function(training[i][0], centers[i][0], training[i][1], centers[i][1])
#        self.design_matrix = design_matrix
#        
#        predicted_value = design_matrix.dot(W)
#        loss = 
#        
#        pass
#        
#    def Train(self):
#        W = self.params['W']
#
#        
#        pass
# Import the data
########################################################################################
######################################################################################
os.chdir("/home/leo/Documents/Database/Pipeline/Ready_2_2_1_1")
#os.listdir()


with open('training_2_2_1_1_RBFN_Mouse', 'r') as f:
    training = json.load(f)
with open('testing_2_2_1_1_Mouse', 'r') as f:
    testing = json.load(f)

training[:66]
testing[:66]

'''
Change the contact into 1 or zero, 1 means they form a core, 0 means they can not {'rate_positive': 0.528,
 'rate_negative': 0.525,
 'truncated_rate_positive': 0.663,
 'truncated_rate_negative': 0.637}
form a core.
'''
new_training = []
for parepi in training:
    new_training.append([parepi[0], parepi[1], 1, parepi[3]])
new_testing = []
for parepi in testing:
    new_testing.append([parepi[0], parepi[1], 1, parepi[3]])
new_training[:6]
len(new_training)
new_testing[:6]
#####################################################################################
######################################################################################

rbfn = RBFN(training,testing, basis_function= 'Markov', distance_mode = 'Multiplication_distance')
rbfn.Centers_from_NewDataAnalysis(center_mode='positive_only')
#rbfn.Pruning(delta= -0.2)
len(rbfn.centers)
#for i in [x*0.05 for x in range(21)]:
#    rbfn.Pruning(delta = i)
#    print(len(rbfn.centers))
rbfn.Formula_ridge_find_reg()
rbfn.GCV_lowest
rbfn.reg_best
rbfn.radius
rbfn.Formula_ridge_predict_Ab()
rbfn.correct_rate
###########################################################################
'''
This block is to do forward selection, it is a little bit slow, need to be polished 
further.
'''
##############################################################################
'''
Do prediction by using different forms of data, adding negative samples.
'''
#rbfn.Add_negative_samples()
len(rbfn.training)
augmented_training = rbfn.training
augmented_training[:6]
augmented_training[-6:-1]
len(augmented_training)
rbfn = RBFN(augmented_training, new_testing, basis_function= 'Markov', distance_mode = 'Multiplication_distance')
rbfn.Centers_from_NewDataAnalysis(center_mode = 'positive_only')
len(rbfn.centers)
rbfn.Formula_ridge_find_reg()
rbfn.GCV_lowest
rbfn.reg_best
rbfn.radius
rbfn.Formula_ridge_predict_Ab()
rbfn.correct_rate
###############################################################################





















