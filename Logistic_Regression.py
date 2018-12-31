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
from matplotlib import pyplot as plt


###########################################################
'''
Import the data
'''
os.chdir("/home/leo/Documents/Database/Pipeline/Ready_2_2_1_1")

with open('ready_2_2_1_1__cn_gate_1_all', 'r') as f:
    data = json.load(f)
len(data)
    
with open('negative_samples', 'r') as f:
    negative_samples = json.load(f)
len(negative_samples)
#############################################################
'''
Generate the testing the samples and the indicator for the testing samples
and the indicator for the testing samples.
'''  
def Generating_testing(positive_samples, negative_samples):
    tp = math.floor(len(positive_samples) * 0.1)
    testing_positive = random.sample(positive_samples, tp)
    for i in testing_positive:
        positive_samples.remove(i)
        
    tn = math.floor(len(negative_samples)*0.1)
    testing_negative = random.sample(negative_samples, tn)
    for i in testing_negative:
        negative_samples.remove(i)
        
    testing = copy.deepcopy(testing_positive)
    testing.extend(testing_negative)
    
    positive = np.ones((tp, 1))
    negative = np.zeros((tn, 1))
    indicator = np.vstack((positive, negative))
    
    return positive_samples, negative_samples, testing, indicator
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

'''
This function should return the vectors required for the input of the Loss function
'''
def Similarity_matrix(testing, positive_samples, negative_samples):
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
            pp_matrix[j, i] = pp_matrix[i, j]
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
    
    # Calculate the vectors for testing

    t = len(testing)
    testing_matrix_positive = np.zeros((t, p))
    testing_matrix_negative = np.zeros((t, n))
    for i in range(t):
        for j in range(p):
            testing_matrix_positive[i, j] = aligner.score(To_seq(testing[i][0]), To_seq(positive_samples[j][0])) 
            + aligner.score(To_seq(testing[i][1]), To_seq(positive_samples[j][1]))
        for k in range(n):
            testing_matrix_negative[i,k] = aligner.score(To_seq(testing[i][0]), To_seq(negative_samples[k][0]))
            + aligner.score(To_seq(testing[i][1]), To_seq(negative_samples[k][1]))
            
    # Calculator the vectors
    testing_sum_positive = np.sum(testing_matrix_positive, axis=1, keepdims = True)
    testing_sum_negative = np.sum(testing_matrix_negative, axis = 1, keepdims = True)
    testing_sum_all = np.hstack((testing_sum_positive, testing_sum_negative))
    # Create the indicator vector
    indicator_positive = np.ones((p, 1))
    indicator_negative = np.zeros((n, 1))
    indicator = np.vstack((indicator_positive, indicator_negative))
            
    return sum_all, indicator, testing_sum_all

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

def Loss_BFGS():
    pass

def Train_model(sum_all, indicator, reg = 0, l_rate = 1e-1, iteration = 10):
    # Normalize the input sum_all
    average = np.average(sum_all, axis = 0)
    std = np.std(sum_all, axis = 0, keepdims = True)
    sum_all -= average
    sum_all /= std
    ones = np.ones((np.shape(sum_all)[0], 1))
    sum_all = np.hstack((sum_all, ones))
    # Initiate the parameters as zeros
    parameter = {}
    parameter['reg'] = reg
    coeff = np.random.normal(0, 1, (3, 1))
    parameter['coeff'] = coeff
    
    loss=[]
    for i in range(iteration):
        temp_loss, grad_coeff = Loss(sum_all, indicator, parameter)
        coeff -= l_rate * grad_coeff
        parameter['coeff'] = coeff
        if iteration%100 == 0:
            print(temp_loss) 
            loss.append(temp_loss)
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
    
    return pred, round(rate, 3)


positive_samples, negative_samples, testing, testing_indicator = Generating_testing(data, negative_samples)
len(positive_samples)
len(negative_samples)
len(testing)
len(testing_indicator)  
sum_all, indicator, testing_sum_all = Similarity_matrix(testing, positive_samples, negative_samples)
parameter, loss = Train_model(sum_all, indicator, reg = 1, l_rate = 1e-6, iteration = 100000)
parameter
#Plot(loss)
prediction, rate = Prediction(testing_sum_all, testing_indicator, parameter)
prediction
rate
##################################################################################
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


