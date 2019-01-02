'''
Logistic Regression
'''
###############################################################
# Import the modules
import random
import numpy as np
import os
import json
import math
import copy
import matplotlib as plt


###########################################################
'''
Import the data
'''
os.chdir("/home/leo/Documents/Database/Pipeline/Ready_2_2_1_1")

with open('ready_2_2_1_1__cn_gate_1_all', 'r') as f:
    data = json.load(f)
len(data)
    
#with open('negative_samples', 'r') as f:
#    negative_samples = json.load(f)
#len(negative_samples)
#negative_samples = random.sample(negative_samples, 1000)
#len(negative_samples)
#############################################################
'''
Generate the testing the samples and the indicator for the testing samples
and the indicator for the testing samples.
'''  
def Generating_testing(positive_samples, negative_samples):
    p = len(positive_samples)
    tp = math.floor(p * 0.1)
    testing_positive = random.sample(positive_samples, tp)
    for i in testing_positive:
        positive_samples.remove(i)
        
    n = len(negative_samples)    
    tn = math.floor(n*0.1)
    testing_negative = random.sample(negative_samples, tn)
    for i in testing_negative:
        negative_samples.remove(i)            
    return positive_samples, negative_samples, testing_positive, testing_negative
####################################################################################################
    ###############################################################################################
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

'''
This function should return the vectors required for the input of the Loss function
'''
def Similarity_matrix(testing_positive, testing_negative, positive_samples, negative_samples, percentage = 0.1):
    p = len(positive_samples)
    n = len(negative_samples)
    pp_matrix = np.zeros((p, p))
    pn_matrix = np.zeros((p, n))
    nn_matrix = np.zeros((n, n))
    np_matrix = np.zeros((n, p))
    for i in range(p):
        for j in range(i, p):
            pp_matrix[i, j] = aligner.score(To_seq(positive_samples[i][0]), To_seq(positive_samples[j][0])) 
            + aligner.score(To_seq(positive_samples[i][1]), To_seq(positive_samples[j][1]))
            pp_matrix[j, i] = pp_matrix[i, j]
        for k in range(n):
            pn_matrix[i,k] = aligner.score(To_seq(positive_samples[i][0]), To_seq(negative_samples[k][0]))
            + aligner.score(To_seq(positive_samples[i][1]), To_seq(negative_samples[k][1]))
            
    for i in range(n):
        for j in range(i, n):
            nn_matrix[i, j] = aligner.score(To_seq(negative_samples[i][0]), To_seq(negative_samples[j][0])) 
            + aligner.score(To_seq(negative_samples[i][1]), To_seq(negative_samples[j][1]))
            nn_matrix[j, i] = nn_matrix[i, j]
        for k in range(p):
            np_matrix[i,k] = aligner.score(To_seq(negative_samples[i][0]), To_seq(positive_samples[k][0]))
            + aligner.score(To_seq(negative_samples[i][1]), To_seq(positive_samples[k][1]))
    # Stack the above matrices
    matrix_positive = np.vstack((pp_matrix, np_matrix))
    matrix_negative = np.vstack((pn_matrix, nn_matrix))
    # Calculator the vectors
    sum_positive = np.sum(matrix_positive, axis=1, keepdims = True)
    sum_negative = np.sum(matrix_negative, axis = 1, keepdims = True)
    sum_all = np.hstack((sum_positive, sum_negative))
    # Truncate 
    training_matrix = np.hstack((matrix_positive, matrix_negative))
    number = math.floor((n+p) * percentage)
    truncated_training_sum = Truncated(training_matrix, p, n, number)
    
    # Create the training_indicator vector
    positive = np.ones((p, 1))
    negative = np.zeros((n, 1))
    training_indicator = np.vstack((positive, negative))
    
    # Calculate the vectors for testing
    tp = len(testing_positive)
    testing_pp = np.zeros((tp, p))
    testing_pn = np.zeros((tp, n))
    for i in range(tp):
        for j in range(p):
            testing_pp[i, j] = aligner.score(To_seq(testing_positive[i][0]), To_seq(positive_samples[j][0])) 
            + aligner.score(To_seq(testing_positive[i][1]), To_seq(positive_samples[j][1]))
        for k in range(n):
            testing_pn[i,k] = aligner.score(To_seq(testing_positive[i][0]), To_seq(negative_samples[k][0]))
            + aligner.score(To_seq(testing_positive[i][1]), To_seq(negative_samples[k][1]))
            
    tn = len(testing_negative)
    testing_np = np.zeros((tn, p))
    testing_nn = np.zeros((tn, n))
    for i in range(tn):
        for j in range(p):
            testing_np[i, j] = aligner.score(To_seq(testing_negative[i][0]), To_seq(positive_samples[j][0])) 
            + aligner.score(To_seq(testing_negative[i][1]), To_seq(positive_samples[j][1]))
        for k in range(n):
            testing_nn[i,k] = aligner.score(To_seq(testing_negative[i][0]), To_seq(negative_samples[k][0]))
            + aligner.score(To_seq(testing_negative[i][1]), To_seq(negative_samples[k][1]))
            
    # Calculator the vectors
    testing_sum_pp = np.sum(testing_pp, axis=1, keepdims = True)
    testing_sum_pn = np.sum(testing_pn, axis = 1, keepdims = True)
    positive_testing_sum = np.hstack((testing_sum_pp, testing_sum_pn))

    
    testing_sum_np = np.sum(testing_np, axis=1, keepdims = True)
    testing_sum_nn = np.sum(testing_nn, axis = 1, keepdims = True)
    negative_testing_sum = np.hstack((testing_sum_np, testing_sum_nn))
    
    # Truncated
    number = math.floor((p+n) * percentage)
    positive_testing_matrix = np.hstack((testing_pp, testing_pn))
    negative_testing_matrix = np.hstack((testing_np, testing_nn))
    truncated_positive_testing_sum = Truncated(positive_testing_matrix, p, n, number)
    truncated_negative_testing_sum = Truncated(negative_testing_matrix, p, n, number)
    
    # Pack all the returned value in a dictionary
    package = {}
    package['training_sum'] = sum_all
    package['training_indicator'] = training_indicator
    package['positive_testing_sum'] = positive_testing_sum
    package['negative_testing_sum'] = negative_testing_sum
    package['truncated_positive_testing_sum'] = truncated_positive_testing_sum
    package['truncated_negative_testing_sum'] = truncated_negative_testing_sum
    package['truncated_training_sum'] = truncated_training_sum

           
    return package

def Distance_matrix(testing_positive, testing_negative, positive_samples, negative_samples, percentage):
    p = len(positive_samples)
    n = len(negative_samples)
    pp_matrix = np.zeros((p, p))
    pn_matrix = np.zeros((p, n))
    nn_matrix = np.zeros((n, n))
    np_matrix = np.zeros((n, p))
    for i in range(p):
        for j in range(i, p):
            pp_matrix[i, j] = Multiplication_distance(To_seq(positive_samples[i][0]), To_seq(positive_samples[j][0]),
                                     To_seq(positive_samples[i][1]), To_seq(positive_samples[j][1])) 
            pp_matrix[j, i] = pp_matrix[i, j]
        for k in range(n):
            pn_matrix[i,k] = Multiplication_distance(To_seq(positive_samples[i][0]), To_seq(negative_samples[k][0]), 
                     To_seq(positive_samples[i][1]), To_seq(negative_samples[k][1]))
            
    for i in range(n):
        for j in range(i, n):
            nn_matrix[i, j] = Multiplication_distance(To_seq(negative_samples[i][0]), To_seq(negative_samples[j][0]),
                     To_seq(negative_samples[i][1]), To_seq(negative_samples[j][1])) 
            nn_matrix[j, i] = nn_matrix[i, j]
        for k in range(p):
            np_matrix[i,k] = Multiplication_distance(To_seq(negative_samples[i][0]), To_seq(positive_samples[k][0]),
                     To_seq(negative_samples[i][1]), To_seq(positive_samples[k][1]))
    # Stack the above matrices
    matrix_positive = - np.vstack((pp_matrix, np_matrix))
    matrix_negative = - np.vstack((pn_matrix, nn_matrix))
    # Calculator the vectors
    sum_positive = np.sum(matrix_positive, axis=1, keepdims = True)
    sum_negative = np.sum(matrix_negative, axis = 1, keepdims = True)
    sum_all = np.hstack((sum_positive, sum_negative))
    # Truncate 
    training_matrix = np.hstack((matrix_positive, matrix_negative))
    number = math.floor((n+p) * percentage)
    truncated_training_sum = Truncated(training_matrix, p, n, number)
    

    # Create the training_indicator vector
    positive = np.ones((p, 1))
    negative = np.zeros((n, 1))
    training_indicator = np.vstack((positive, negative))
    
    # Calculate the vectors for testing

    tp = len(testing_positive)
    testing_pp = np.zeros((tp, p))
    testing_pn = np.zeros((tp, n))
    for i in range(tp):
        for j in range(p):
            testing_pp[i, j] = Multiplication_distance(To_seq(testing_positive[i][0]), To_seq(positive_samples[j][0]),
                      To_seq(testing_positive[i][1]), To_seq(positive_samples[j][1])) 

        for k in range(n):
            testing_pn[i,k] = Multiplication_distance(To_seq(testing_positive[i][0]), To_seq(negative_samples[k][0]),
                      To_seq(testing_positive[i][1]), To_seq(negative_samples[k][1]))

            
    tn = len(testing_negative)
    testing_np = np.zeros((tn, p))
    testing_nn = np.zeros((tn, n))
    for i in range(tn):
        for j in range(p):
            testing_np[i, j] = Multiplication_distance(To_seq(testing_negative[i][0]), To_seq(positive_samples[j][0]),
                      To_seq(testing_negative[i][1]), To_seq(positive_samples[j][1])) 
        for k in range(n):
            testing_nn[i,k] = Multiplication_distance(To_seq(testing_negative[i][0]), To_seq(negative_samples[k][0]),
                      To_seq(testing_negative[i][1]), To_seq(negative_samples[k][1]))
            
    # Calculator the vectors
    testing_pp = - testing_pp
    testing_pn = -testing_pn
    testing_sum_pp = np.sum(testing_pp, axis=1, keepdims = True)
    testing_sum_pn = np.sum(testing_pn, axis = 1, keepdims = True)
    positive_testing_sum = np.hstack((testing_sum_pp, testing_sum_pn))
    
    testing_np = - testing_np
    testing_nn = - testing_nn
    testing_sum_np = np.sum(testing_np, axis=1, keepdims = True)
    testing_sum_nn = np.sum(testing_nn, axis = 1, keepdims = True)
    negative_testing_sum = np.hstack((testing_sum_np, testing_sum_nn))
    
    # Truncated
    number = math.floor((p+n) * percentage)
    positive_testing_matrix = np.hstack((testing_pp, testing_pn))
    negative_testing_matrix = np.hstack((testing_np, testing_nn))
    truncated_positive_testing_sum = Truncated(positive_testing_matrix, p, n, number)
    truncated_negative_testing_sum = Truncated(negative_testing_matrix, p, n, number)
    
    # Pack all the returned value in a dictionary
    package = {}
    package['training_sum'] = sum_all
    package['training_indicator'] = training_indicator
    package['positive_testing_sum'] = positive_testing_sum
    package['negative_testing_sum'] = negative_testing_sum
    package['truncated_positive_testing_sum'] = truncated_positive_testing_sum
    package['truncated_negative_testing_sum'] = truncated_negative_testing_sum
    package['truncated_training_sum'] = truncated_training_sum
           
    return package
'''
Truncate the logistic regression, make sure the pair is only influenced by the nearby samples
Input:
    matrix, a matrix in the shape of (s, p+n), where s is the number of samples, p is the number of 
    positive training samples and n is the number of negative training samples.
    number, an integre, gives the number of neareast 'number' samples
Output:
    sum_all_truncated, in the same form of sum_all
'''
def Truncated(matrix, p, n, number):
    s = np.shape(matrix)[0]
    sum_all_p = np.zeros((s, 1))
    sum_all_n = np.zeros((s, 1))
    for i in range(s):
        similarity = matrix[i, :]
        # mark the samples with positive or negative signs
        similarity_signed = []
        for j in range(p):
            similarity_signed.append([similarity[j], 1])
        for k in range(p, p+n):
            similarity_signed.append([similarity[k], -1])
        similarity_signed.sort(key = lambda x:x[0], reverse = True)
        # select the top number
        similarity_selected = similarity_signed[:number]
        # Calculate the sums for different signs
        sum_p = 0
        sum_n = 0
        for smile in similarity_selected:
            if smile[1] == 1:
                sum_p += smile[0]
            elif smile[1] == -1:
                sum_n  += smile[0]
        # load to the sum_all_p and sum_all_n
        sum_all_p[i, 0] = sum_p  
        sum_all_n[i, 0] = sum_n 
    # Generate the sum_all_truncated
    sum_all_truncated = np.hstack((sum_all_p, sum_all_n))
    return sum_all_truncated
#######################################################################
#####################################################################
#import numpy as np
#a = list(range(16))
#b = np.reshape(a, (4, 4))
#np.shape(b)
#np.shape(b)[0]
#b[0,:]
#sum_all = Truncated(b, 2, 2, 3)
#b
#sum_all
######################################################################################
    #################################################################################

def Loss(sum_all, indicator, parameter):
    coeff = parameter['coeff']
    reg = parameter['reg']
    # Calculate the cross entrapy loss
    logit = (sum_all).dot(coeff)
    exp_logit = np.exp(logit)
    prob = exp_logit/(1+exp_logit)
    ones = np.ones_like(indicator)
    cross_entropy = - (indicator.T).dot(logit) - (ones.T).dot(np.log(1-prob))      
    # Calculate the regularization loss
    loss_reg = reg * (coeff.T).dot(coeff)
    # Calculate the total loss
    loss_total = cross_entropy + loss_reg
    
    # Calculate the up date gradient
    d_logit = - indicator + prob
    d_coeff = (sum_all.T).dot(d_logit)
    d_coeff += 2 * reg * coeff                          
    pass
    return loss_total, d_coeff


def Train_model(sum_all, indicator, parameter, l_rate = 1e-1, iteration = 10):
    # Normalize the input sum_all
    average = np.average(sum_all, axis = 0)
    std = np.std(sum_all, axis = 0, keepdims = True)
    sum_all -= average
    sum_all /= std
    ones = np.ones((np.shape(sum_all)[0], 1))
    sum_all = np.hstack((sum_all, ones))
    # take out the coeff
    coeff = parameter['coeff']     
    loss=[]
    for i in range(iteration):
        temp_loss, grad_coeff = Loss(sum_all, indicator, parameter)
        coeff -= l_rate * grad_coeff
        parameter['coeff'] = coeff
        if iteration%1000 == 0:
            print(temp_loss) 
            loss.append(temp_loss)
    return parameter, loss

def Train_BFGS(sum_all, indicator, reg = 1, rho=0.9, c = 1e-3, termination = 1e-2):
    # Normalize the input sum_all
    average = np.average(sum_all, axis = 0)
    std = np.std(sum_all, axis = 0, keepdims = True)
    sum_all -= average
    sum_all /= std
    ones = np.ones((np.shape(sum_all)[0], 1))
    sum_all = np.hstack((sum_all, ones))
    # Give the initial Hessian h
    H = np.eye(3)*1e-5
    # Set the starting point
    coeff = np.zeros((3, 1))
    parameter = {}
    parameter['coeff'] = coeff
    parameter['reg'] = reg
   # BFGS algorithm
    loss, grad = Loss(sum_all, indicator, parameter)
    grad_square = (grad.T).dot(grad)
    while grad_square >= termination**2:        
        p = - H.dot(grad)        
        # Find the next coeff
        parameter_new = {}
        parameter_new['reg'] = reg
        parameter_new['coeff'] = p + parameter['coeff']
        
        new_loss, new_grad = Loss(sum_all, indicator, parameter_new)

        while new_loss > loss + c * (grad.T).dot(p):
            p *= rho
            print(p, 'p')
            parameter_new['coeff'] = p + parameter['coeff']
            new_loss, new_grad = Loss(sum_all, indicator, parameter_new)
        
        # update H
        s = p
        y = new_grad - grad
        r = (y.T).dot(s)
        if r != 0:
            r = 1/r
            I = np.eye(3)
            H = (I - r*s.dot(y.T)).dot(H).dot(I - r*y.dot(s.T)) + r*s.dot(s.T)
        else:
            H = I
        # Update loss, grad, grad_square and paramter
        loss = new_loss
        grad = new_grad
        parameter['coeff'] = parameter_new['coeff']
        grad_square = (grad.T).dot(grad)
        print(loss, '    ', grad_square)
    return parameter, loss

def Plot(loss):
    iteration  = list(range(len(loss)))
    plt(iteration, loss)
    plt.show()
    pass
    return 

def Prediction(testing_sum_all, indicator,  parameter):
    # Normalize the testing_sum_all matrix
    sum_all = testing_sum_all
    average = np.average(sum_all, axis = 0)
    std = np.std(sum_all, axis = 0, keepdims = True)
    sum_all -= average
    sum_all /= std
    ones = np.ones((np.shape(sum_all)[0], 1))
    sum_all = np.hstack((sum_all, ones))
    # Do the prediction
    n_total = np.shape(testing_sum_all)[0]
    coeff = parameter['coeff']
    pred = sum_all.dot(coeff)
    n = 0
    for i in range(n_total):
        if pred[i] > 0 and indicator[i] == 1:
            n += 1
        elif pred[i] < 0 and indicator[i] == 0:
            n += 1
    rate = n / n_total
    
    return  round(rate, 3)
#######################################
'''
Generate random samples by random process
'''
def Generate_random_negative(sample_size, data, Ab_lenth, Ag_length):
    
    TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
                    ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
                    ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]
    AA = []
    for aa in TripleSingle:
        AA.append(aa[0])
        
    Ab = []
    Ag = []
    for parepi in data:
        Ab.append(parepi[0])
        Ag.append(parepi[1])
    Ab_Ag = copy.deepcopy(Ab)
    Ab_Ag.extend(Ag)
        
    negative_samples = []
    while len(negative_samples) < sample_size:
        r_Ab = ''
        while r_Ab == '':        
            r_Ab = random.sample(AA, Ab_lenth)
            if r_Ab in Ab_Ag:
                r_Ab = ''
        r_Ag = ''
        while r_Ag == '':
            r_Ag = random.sample(AA, Ag_length)
            if r_Ag in Ab_Ag:
                r_Ag = ''
        negative_samples.append([r_Ab, r_Ag, 0, 0])
    return negative_samples
##################################################################################
    #######################################################################
'''
This code block is used to generate training samples and testing samples, both negative
and positive
'''
with open('ready_2_2_1_1__cn_gate_1_all', 'r') as f:
    data = json.load(f)
positive_data = []
for parepi in data:
    if parepi[2] >= 16:
        positive_data.append(parepi)

len(positive_data)
positive_data[:6]

negative_samples = Generate_random_negative(800, positive_data, 2, 2)
len(negative_samples)
negative_samples[:6]

positive_samples, negative_samples, testing_positive, testing_negative = Generating_testing(positive_data, negative_samples)
len(positive_samples)
len(negative_samples)
len(testing_positive)
len(testing_negative) 
#####################################################################################3
#################################################################################### 


#parameter = {}
#parameter['reg'] = 1
#coeff = np.random.normal(0, 1, (3, 1))
#parameter['coeff'] = coeff
#sum_all = package['sum_all']
#training_indicator = package['training_indicator']
#parameter, loss = Train_model(sum_all, training_indicator, parameter, l_rate = 1e-6, iteration = 10000)
#parameter
#loss

##########################################
def Rate(package):
    parameter_BFGS, loss_BFGS = Train_BFGS(package['training_sum'], package['training_indicator'], 
                                           reg=0, rho=0.85, c=1e-2, termination=1e-5)
#    parameter_BFGS
#    loss_BFGS
    #Plot(loss)
    rate = {}
    positive_indicator = np.ones((len(testing_positive), 1))
    rate_positive = Prediction(package['positive_testing_sum'], positive_indicator, parameter_BFGS)
    rate['rate_positive'] = rate_positive
#    rate
#    prediction
    
    negative_indicator = np.zeros((len(testing_negative), 1))
    rate_negative = Prediction(package['negative_testing_sum'], negative_indicator, parameter_BFGS)
    rate['rate_negative'] = rate_negative
#    rate
    
    ##################################################################################
    '''
    Truncated logistic
    '''
    parameter_BFGS, loss_BFGS = Train_BFGS(package['truncated_training_sum'], package['training_indicator'], 
                                           reg=0, rho=0.85, c=1e-2, termination=1e-5)
#    parameter_BFGS
#    loss_BFGS
    #Plot(loss)
    positive_indicator = np.ones((len(testing_positive), 1))
    truncated_rate_positive = Prediction(package['truncated_positive_testing_sum'], positive_indicator, parameter_BFGS)   
    rate['truncated_rate_positive'] = truncated_rate_positive
    
    negative_indicator = np.zeros((len(testing_negative), 1))
    truncated_rate_negative = Prediction(package['truncated_negative_testing_sum'], negative_indicator, parameter_BFGS)
    rate['truncated_rate_negative'] = truncated_rate_negative
    
    return rate

#####################################################################################
package = Similarity_matrix(testing_positive, testing_negative, positive_samples, negative_samples, percentage=0.05)
package = Distance_matrix(testing_positive, testing_negative, positive_samples, negative_samples, percentage=0.05)
keys = list(package.keys())
keys
len(package['training_sum'])
rate = Rate(package)
rate
#######################################
######################################
'''
A little bit of experimental data
'''


sum_all_testing = np.random.normal(0, 1, (5, 2))*2 + 1
sum_all_tsting = np.ones((5, 2))
indicator_testing = np.ones((5, 1))
#Train_model(sum_all_testing, indicator_testing, reg = 0, l_rate = 1e-1, iteration = 150)
average = np.average(sum_all_testing, axis = 0)
std = np.std(sum_all_testing, axis = 0, keepdims = True)
sum_all_testing -= average
sum_all_testing /= std
ones = np.ones((np.shape(sum_all_testing)[0], 1))
sum_all_testing = np.hstack((sum_all_testing, ones))
sum_all_testing
# Initiate the parameters as zeros
reg = 0
parameter = {}
parameter['reg'] = reg
coeff = np.random.normal(0, 1, (3, 1))
parameter['coeff'] = coeff
coeff

temp_loss, grad_coeff = Loss(sum_all_testing, indicator_testing, parameter)
temp_loss

Train_model(sum_all_testing, indicator_testing, reg = 0, l_rate = 1e-1, iteration = 150)

iteration = 150
l_rate = 1e-1
loss=[]
for i in range(iteration):
    temp_loss, grad_coeff = Loss(sum_all_testing, indicator_testing, parameter)
    loss.append(temp_loss)
    coeff -= l_rate * grad_coeff
    parameter['coeff'] = coeff
    print(temp_loss) 


#for i in range(150):
#    temp_loss, grad_coeff = Loss(sum_all_testing, indicator_testing, parameter_testing)
#    loss_testing.append(temp_loss)
#    coeff -= 1e-1 * grad_coeff
#    parameter_testing['coeff'] = coeff
#    print(temp_loss) 
  
#
#np.hstack((a, b, c.T))
#np.sum(a, axis = 1, keepdims = True)


