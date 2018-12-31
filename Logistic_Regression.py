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
    
with open('negative_samples', 'r') as f:
    negative_samples = json.load(f)
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
'''
This function should return the vectors required for the input of the Loss function
'''
def Similarity_matrix(testing, training_positive_samples, training_negative_samples):
    all_data = copy.deepcopy(training_positive_samples)
    all_data.extend(training_negative_samples)
    n = len(all_data)
    np = len(training_positive_samples)
    nn = len(training_negative_samples)
    training_similarity_matrix_positive = np.zeros((n, np))
    training_similarity_matrix_negative =  np.zeros((n, nn))
    for i in all_data:
        for j in training_positive_samples:
            training_similarity_matrix_positive[i, j] = aligner.score(To_seq(i[:2]), To_seq(j[:2])) 
        for k in training_negative_samples:
            training_similarity_matrix_negative[i,k] = aligner.score(To_seq(i[:2], To_seq(j[:2])))
    # Calculator the vectors
    ones = np.ones((n, 1))
    sum_positive = np.sum(training_similarity_matrix_positive, axis=1)
    sum_negative = np.sum(training_similarity_matrix_negative, axis = 1)
    sum_all = np.hstack(sum_positive, sum_negative, ones)
    
    # Calculate the vectors for testing

    nt = len(testing)
    testing_similarity_matrix_positive = np.zeros((nt, np))
    testing_similarity_matrix_negative = np.zeros((nt, nn))
    for i in testing:
        for j in training_positive_samples:
            testing_similarity_matrix_positive[i, j] = aligner.score(To_seq(i[:2]), To_seq(j[:2])) 
        for k in training_negative_samples:
            testing_similarity_matrix_negative[i,k] = aligner.score(To_seq(i[:2], To_seq(j[:2])))
            
    # Calculator the vectors
    ones = np.ones((nt, 1))
    testing_sum_positive = np.sum(testing_similarity_matrix_positive, axis=1)
    testing_sum_negative = np.sum(testing_similarity_matrix_negative, axis = 1)
    testing_sum_all = np.hstack(testing_sum_positive, testing_sum_negative, ones)
    
    indicator_positive = np.ones((np, 1))
    indicator_negative = np.zeros((nn, 1))
    indicator = np.vstack(indicator_positive, indicator_negative)
            
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
    loss_reg = reg * coeff.dot(coeff)
    # Calculate the total loss
    loss_total = cross_entropy + loss_reg
    
    # Calculate the up date gradient
    d_logit = - indicator + prob
    d_coeff = (sum_all.T).dot(d_logit)
    d_coeff += 2 * reg * coeff                    
        
    pass
    return loss_total, d_coeff

def Train_model(sum_all, indicator, reg = 0, l_rate = 1e-3, iteration = 10):
    # Initiate the parameters as zeros
    parameter = {}
    parameter['reg'] = reg
    parameter['coeff'] = np.zeros((3, 0))
    
    loss=[]
    for i in range(iteration):
        temp_loss, grad_coeff = Loss(sum_all, indicator, parameter)
        loss.append(temp_loss)
        parameter['coeff'] -= l_rate * grad_coeff        
    return parameter, loss

def Plot(loss):
    iteration  = list(range(len(loss)))
    plt(iteration, loss)
    plt.show()
    pass
    return 

def Prediction():
    pass




