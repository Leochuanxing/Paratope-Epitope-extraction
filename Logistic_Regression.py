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
import matplotlib.pyplot as plt


###########################################################


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
def Similarity_matrix(testing_positive, testing_negative, positive_samples, negative_samples):
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
#    # Calculator the vectors
#    sum_positive = np.sum(matrix_positive, axis=1, keepdims = True)
#    sum_negative = np.sum(matrix_negative, axis = 1, keepdims = True)
#    sum_all = np.hstack((sum_positive, sum_negative))
    # Truncate 
    training_matrix = np.hstack((matrix_positive, matrix_negative))
#    number = math.floor((n+p) * percentage)
#    truncated_training_sum = Truncated(training_matrix, p, n, number)
    truncate_package = {}
    truncate_package['training_matrix'] = training_matrix
    truncate_package['p'] = p
    truncate_package['n'] = n

    
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
#    testing_sum_pp = np.sum(testing_pp, axis=1, keepdims = True)
#    testing_sum_pn = np.sum(testing_pn, axis = 1, keepdims = True)
#    positive_testing_sum = np.hstack((testing_sum_pp, testing_sum_pn))
#
#    
#    testing_sum_np = np.sum(testing_np, axis=1, keepdims = True)
#    testing_sum_nn = np.sum(testing_nn, axis = 1, keepdims = True)
#    negative_testing_sum = np.hstack((testing_sum_np, testing_sum_nn))
#    
    # Truncated
#    number = math.floor((p+n) * percentage)
    positive_testing_matrix = np.hstack((testing_pp, testing_pn))
    negative_testing_matrix = np.hstack((testing_np, testing_nn))
#    truncated_positive_testing_sum = Truncated(positive_testing_matrix, p, n, number)
#    truncated_negative_testing_sum = Truncated(negative_testing_matrix, p, n, number)
    truncate_package['positive_testing_matrix'] = positive_testing_matrix
    truncate_package['negative_testing_matrix'] = negative_testing_matrix
    
    # Pack all the returned value in a dictionary
#    package = {}
#    package['training_sum'] = sum_all
#    package['training_indicator'] = training_indicator
#    package['positive_testing_sum'] = positive_testing_sum
#    package['negative_testing_sum'] = negative_testing_sum
#    package['truncated_positive_testing_sum'] = truncated_positive_testing_sum
#    package['truncated_negative_testing_sum'] = truncated_negative_testing_sum
#    package['truncated_training_sum'] = truncated_training_sum

           
    return training_indicator, truncate_package

def Distance_matrix(testing_positive, testing_negative, positive_samples, negative_samples):
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
#    sum_positive = np.sum(matrix_positive, axis=1, keepdims = True)
#    sum_negative = np.sum(matrix_negative, axis = 1, keepdims = True)
#    sum_all = np.hstack((sum_positive, sum_negative))
    # Truncate 
    training_matrix = np.hstack((matrix_positive, matrix_negative))
#    number = math.floor((n+p) * percentage)
#    truncated_training_sum = Truncated(training_matrix, p, n, number)
    truncate_package = {}
    truncate_package['training_matrix'] = training_matrix
    truncate_package['p'] = p
    truncate_package['n'] = n

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
#    testing_sum_pp = np.sum(testing_pp, axis=1, keepdims = True)
#    testing_sum_pn = np.sum(testing_pn, axis = 1, keepdims = True)
#    positive_testing_sum = np.hstack((testing_sum_pp, testing_sum_pn))
    
    testing_np = - testing_np
    testing_nn = - testing_nn
#    testing_sum_np = np.sum(testing_np, axis=1, keepdims = True)
#    testing_sum_nn = np.sum(testing_nn, axis = 1, keepdims = True)
#    negative_testing_sum = np.hstack((testing_sum_np, testing_sum_nn))
    
    # Truncated
#    number = math.floor((p+n) * percentage)
    positive_testing_matrix = np.hstack((testing_pp, testing_pn))
    negative_testing_matrix = np.hstack((testing_np, testing_nn))
#    truncated_positive_testing_sum = Truncated(positive_testing_matrix, p, n, number)
#    truncated_negative_testing_sum = Truncated(negative_testing_matrix, p, n, number)
    truncate_package['positive_testing_matrix'] = positive_testing_matrix
    truncate_package['negative_testing_matrix'] = negative_testing_matrix
    
    # Pack all the returned value in a dictionary
#    package = {}
#    package['training_sum'] = sum_all
#    package['training_indicator'] = training_indicator
#    package['positive_testing_sum'] = positive_testing_sum
#    package['negative_testing_sum'] = negative_testing_sum
#    package['truncated_positive_testing_sum'] = truncated_positive_testing_sum
#    package['truncated_negative_testing_sum'] = truncated_negative_testing_sum
#    package['truncated_training_sum'] = truncated_training_sum
#           
    return training_indicator, truncate_package
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

##################################################################################



##################################################################################
    #######################################################
def main():
    
    os.chdir("/home/leo/Documents/Database/Pipeline/Ready_2_2_1_1")
    
    
    with open('training_2_2_1_1_all_jump', 'r') as f:
        positive_samples = json.load(f)
    with open('testing_2_2_1_1_all_jump', 'r') as f:
        testing_positive = json.load(f)
    with open('training_negative', 'r') as f:
        negative_samples = json.load(f)
    with open('testing_negative', 'r') as f:
        testing_negative = json.load(f) 
        
    results = {}#store the results
    truncate_percentage = [0.008, 0.005, 0.002, 0.001]
    for i in range(2):
        if i == 0:           
            training_indicator, truncate_package = Similarity_matrix(testing_positive,\
                                    testing_negative, positive_samples, negative_samples)
        else:
            training_indicator, truncate_package = Distance_matrix(testing_positive,\
                                    testing_negative, positive_samples, negative_samples)
    
        p = truncate_package['p']
        n = truncate_package['n']
        training_matrix = truncate_package['training_matrix']
        positive_testing_matrix = truncate_package['positive_testing_matrix']
        negative_testing_matrix = truncate_package['negative_testing_matrix']
        
        # Store the prediction rate in rate
        rate = []
        for percentage in truncate_percentage:
            number = math.floor((n+p) * percentage)
            truncated_training_sum = Truncated(training_matrix, p, n, number)
            truncated_positive_testing_sum = Truncated(positive_testing_matrix, p, n, number)
            truncated_negative_testing_sum = Truncated(negative_testing_matrix, p, n, number)
            
            
            parameter_BFGS, loss_BFGS = Train_BFGS(truncated_training_sum, training_indicator, 
                                                   reg=0, rho=0.85, c=1e-2, termination=1e-3)
            
            # stack the truncated positive and truncated negative testing matrix
            truncated_testing_sum = np.vstack((truncated_positive_testing_sum, truncated_negative_testing_sum))
            # Generate the indicator
            positive_indicator = np.ones((len(testing_positive), 1))
            negative_indicator = np.zeros((len(testing_negative), 1))
            testing_indicator = np.vstack((positive_indicator, negative_indicator))
            # Do the prediction.
            testing_rate = Prediction(truncated_testing_sum, testing_indicator, parameter_BFGS)   
            rate.append(testing_rate)
            
        if i == 0:       
            results['Addition'] = rate
        else:
            results['Multiplication'] = rate
        
    return truncate_percentage, results
        

if __name__ == '__main__':
    truncate_percentage, results = main()

#results
os.chdir("/home/leo/Documents/Database/Pipeline/Ready_2_2_1_1")
#results['truncate_percentage'] = truncate_percentage
#with open('Logistic_results', 'w') as f:
#    json.dump(result, f)
##
with open('Logistic_results', 'r') as f:
    result = json.load( f)  

result
#result['Addition'].extend(results['Addition'])
#result['Addition']
#result['Multiplication'].extend(results['Multiplication'])
#result['Multiplication']
#result['truncate_percentage'].extend(results['truncate_percentage'])
#result['truncate_percentage']
############################################################
'''
Graph the result
'''
#log_percentage = []
#for percentage in result['truncate_percentage']:
#    log_percentage.append(-math.log(percentage))
#log_percentage
#
#result['Addition']
#plt.figure(figsize = (8, 6))
#plt.plot(log_percentage, result['Addition'], 'r--')
#plt.plot(log_percentage, result['Multiplication'])
#plt.plot(log_percentage, result['Addition'], 'go')
#plt.plot(log_percentage, result['Multiplication'], 'bo')
#plt.legend(['Addition_distance', 'Multiplication_distance'])
#plt.xlabel('-log(percentage)', fontsize = 20)
#plt.ylabel('Accuracy', fontsize = 20)
#plt.ylabel()
#plt.show()


