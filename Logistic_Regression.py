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

def Multiplication_similarity(Ab_seq1, Ab_seq2, Ag_seq1, Ag_seq2):
    l_Ab = len(Ab_seq1)
    l_Ag = len(Ag_seq1)
    similarity = ((4*l_Ab + aligner.score(Ab_seq1, Ab_seq2))/(15*l_Ab)) *  ((4*l_Ag + aligner.score(Ag_seq1, Ag_seq2))/(15*l_Ag))
    return similarity
    
#def Addition_similarity(Ab_seq1, Ab_seq2, Ag_seq1, Ag_seq2):
#    l_Ab = len(Ab_seq1)
#    l_Ag = len(Ag_seq1)
#    similarity = ((4*l_Ab + aligner.score(Ab_seq1, Ab_seq2))/(15*l_Ab)) + ((4*l_Ag + aligner.score(Ag_seq1, Ag_seq2))/(15*l_Ag))
##    distance = -(aligner.score(Ab_seq1, Ab_seq2) + aligner.score(Ag_seq1, Ag_seq2))/(l_Ab*11 + l_Ag * 11)
#    return similarity
##############################################################################
'''
Similarity_matrix:
    a function to calculate the similarity_matrix by the function Multiplication_similarity
INput:
    row_set:
        a list gives the rows of the matrix
    column_set:
        a list gives the columns of the matrix
    square:
        boolean, if it is True, it means the row_set is the same as the column set. Otherwise
                 the row_set and the column_set are taken as the same
output:
    similarity_matrix:
        a matrix in the shape of (len(row_set), len(column_set))
'''
def Similarity_matrix(row_set, column_set, square = False):
    # Create the empty similarity_matrix
    n_row =  len(row_set)
    n_col = len(column_set)
    similarity_matrix = np.zeros((n_row, n_col))
    # Load the values to the matrix
    if square == False:
        for i in range(n_row):
            for j in range(n_col):
                Ab_seq1 = row_set[i][0]
                Ab_seq2 = column_set[j][0]
                Ag_seq1 = row_set[i][1]
                Ag_seq2 = column_set[j][1]
                similarity_matrix[i,j] = Multiplication_similarity(Ab_seq1, Ab_seq2, Ag_seq1, Ag_seq2)
    if square == True:
        for i in range(n_row):
            for j in range(i-1, n_row):
                Ab_seq1 = row_set[i][0]
                Ab_seq2 = column_set[j][0]
                Ag_seq1 = row_set[i][1]
                Ag_seq2 = column_set[j][1]
                similarity_matrix[i,j] = Multiplication_similarity(Ab_seq1, Ab_seq2, Ag_seq1, Ag_seq2)
                similarity_matrix[j,i] = similarity_matrix[i,j]
                
    return similarity_matrix
                
        
###########################################################################

def Load_truncate_package(positive_testing_set, negative_testing_set, positive_training_set, negative_training_set):
    p = len(positive_training_set)
    n = len(negative_training_set)
    # Combine the positive_training_set and negative_training_set
    training_set = copy.deepcopy(positive_training_set)
    training_set.extend(negative_training_set)
    # Load the training_matrix
    training_similarity_matrix = Similarity_matrix(training_set, training_set, square=True)


    # Create the training_indicator vector
    positive = np.ones((p, 1))
    negative = np.zeros((n, 1))
    training_indicator = np.vstack((positive, negative))
    
    # Calculating the testing similarity matrix
    negative_testing_matrix = Similarity_matrix(negative_testing_set, training_set, square = False)
    positive_testing_matrix = Similarity_matrix(positive_testing_set, training_set, square = False)

    # Load the truncate package    
    truncate_package = {}
    truncate_package['training_matrix'] = training_similarity_matrix
    truncate_package['p'] = p
    truncate_package['n'] = n
    truncate_package['positive_testing_matrix'] = positive_testing_matrix 
    truncate_package['negative_testing_matrix'] = negative_testing_matrix
    truncate_package['training_indicator'] = training_indicator
#    truncate_package['training_indicator'] = training_indicator           
    return truncate_package
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
        # To make the calculation more stable, we take the average of similarity values
        # Track the number of positive and negative neighbor by p and n respectively.
        p = 0
        n = 0
        for smile in similarity_selected:
            if smile[1] == 1:
                sum_p += smile[0]
                p += 1
            elif smile[1] == -1:
                sum_n  += smile[0]
                n += 1 
        # load to the sum_all_p and sum_all_n
        if p != 0:            
            sum_all_p[i, 0] = sum_p / p 
        if n != 0:
            sum_all_n[i, 0] = sum_n / n 
    # Generate the sum_all_truncated
    sum_all_truncated = np.hstack((sum_all_p, sum_all_n))
    return sum_all_truncated
#######################################################################

######################################################################################

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

################################################################################
    
def Train_BFGS(sum_all, indicator, reg = 1, rho=0.9, c = 1e-3, termination = 1e-2):
    # Normalize the input sum_all
#    average = np.average(sum_all, axis = 0)
#    std = np.std(sum_all, axis = 0, keepdims = True)
#    sum_all /= sum_all.shape[1] 
#    sum_all /= std
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

##################################################################

#######################################
def Prediction_area(testing_sum_all, indicator, parameter):
    # Normalize the testing_sum_all matrix
    sum_all = testing_sum_all
#    average = np.average(sum_all, axis = 0)
#    std = np.std(sum_all, axis = 0, keepdims = True)
#    sum_all /= sum_all.shape[0]
#    sum_all /= std
    ones = np.ones((np.shape(sum_all)[0], 1))
    sum_all = np.hstack((sum_all, ones))
    # Do the prediction
    n_total = np.shape(testing_sum_all)[0]
    coeff = parameter['coeff']
    pred = sum_all.dot(coeff)
    # Relate the pred with the indicator
    pred_indicator = []
    for i in range(n_total):
        pred_indicator.append([pred[i,0], indicator[i]])
    pred_indicator.sort(key=lambda x:x[0], reverse = True)         
        
    # Calculate the area        
    area = 0
    n_positive = 0
    n_negative = 0
    for i in pred_indicator:
        if i[1] == 1:
            n_positive += 1
        else:
            n_negative += 1
        area += n_positive/(n_positive+n_negative)
    
    return  area
#################################################################
def Precision_recall(testing_sum, testing_indicator, parameter):
    # Attach the constant term
    ones = np.ones((np.shape(testing_sum)[0], 1))
    sum_all = np.hstack((testing_sum, ones))
    # Take out the coeff
    coeff = parameter['coeff']
    # Make prediction
    pred = sum_all.dot(coeff)
    # Relate the pred with the indicator
    pred_indicator = []
    for i in range(len(testing_indicator)):
        pred_indicator.append([pred[i,0], testing_indicator[i]])
    pred_indicator.sort(key=lambda x:x[0], reverse = True)

    # Calculate the recall and precision
    precision=[]
    recall = []
    denominator = np.sum(np.array(testing_indicator))
    n_positive = 0
    n_negative = 0
    for i in pred_indicator:
        if i[1] == 1:# if we do the binary, we can do if i[1] ==1:
            n_positive += 1
        else:
            n_negative += 1
        precision.append(n_positive/(n_positive+n_negative))
        recall.append(n_positive/denominator)
        
    return precision, recall

##################################################################################


def main(wd_results, wd_negative_samples, header):
    # Indicate the process
    print('working on ' + header)    
    os.chdir(wd_results)
    with open(header+'_cross_results_LogisticRegression', 'r') as f:
        results = json.load(f)
    with open(header+'_aa_train', 'r') as f:
        positive_training_set = json.load(f)
    with open(header+'_aa_test', 'r') as f:
        positive_testing_set = json.load(f)
        
    os.chdir(wd_negative_samples)
    with open(header+'_train_negative', 'r') as f:
        negative_training_set = json.load(f)
    with open(header+'_test_negative', 'r') as f:
        negative_testing_set = json.load(f)  
    
    # Find the best parameter
    precentages = results['precentages']
    areas = results['areas']
    best_percentage =precentages[np.argmax(np.array(areas))]
    
    # Load the package
    truncate_package = Load_truncate_package(positive_testing_set,\
                                             negative_testing_set, positive_training_set, negative_training_set)
    # Truncate
    p = truncate_package['p']
    n = truncate_package['n']
    number = math.floor((n+p) * best_percentage)
    training_similarity_matrix = truncate_package['training_matrix']
    positive_testing_matrix = truncate_package['positive_testing_matrix']  
    negative_testing_matrix = truncate_package['negative_testing_matrix']
    training_indicator = truncate_package['training_indicator']
    
    truncated_training_sum = Truncated(training_similarity_matrix , p, n, number)
    truncated_positive_testing_sum = Truncated(positive_testing_matrix, p, n, number)
    truncated_negative_testing_sum = Truncated(negative_testing_matrix, p, n, number)
    
    # Train the model
    parameter_BFGS, loss_BFGS = Train_BFGS(truncated_training_sum, training_indicator, 
                                       reg=1, rho=0.85, c=1e-2, termination=1e-3)
    
    # Prepare the testing data
    truncated_testing_sum = np.vstack((truncated_positive_testing_sum, truncated_negative_testing_sum))
    testing_indicator = np.vstack((np.ones((len(positive_testing_set), 1)), np.zeros((len(negative_testing_set), 1))))
    
    # Make prediction about the testing data
    precision, recall = Precision_recall(truncated_testing_sum, testing_indicator, parameter_BFGS)
    
    # Load the results to a dictionary and return this dictionary
    results = {}
    results['best_percentage'] = best_percentage
    results['precision'] = precision
    results['recall'] = recall
    results['coeff'] = parameter_BFGS['coeff']
    
    return results
##################################################################### 
ordered_headers = ['2_2','3_3', '2_3', '3_2', '3_4', '4_3','2_4', '4_2','4_4',\
                   '1_1','1_2', '2_1', '1_3', '3_1', '1_4', '4_1']
############################################################################
'''
Run the main and Save the results
'''
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
###################################################################
if __name__ == '__main__':
    directory_paires = [['/home/leo/Documents/Database/Pipeline/Results/1_free',\
      '/home/leo/Documents/Database/Pipeline/Negative samples/1_free'],\
     ['/home/leo/Documents/Database/Pipeline/Results/0_free',\
                           '/home/leo/Documents/Database/Pipeline/Negative samples/0_free']]
    
    # Let check if all the corresponding length of the positive training and the negative training
    # are the same. And we have to make sure that all the files needed exist.
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
                with open(header+'_aa_test', 'r') as f:
                    positive_testing_set = json.load(f)
                with open(header+'_cross_results_LogisticRegression', 'r') as f:
                    results = json.load(f)
                    
                os.chdir(wd_negative_samples)
                with open(header+'_train_negative', 'r') as f:
                    negative_training_set = json.load(f)
                with open(header+'_test_negative', 'r') as f:
                    negative_testing_set = json.load(f)
                
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
            results =  main(wd_results, wd_negative_samples, header)
            
             # Save the results
            os.chdir(wd_results)
            with open(header+'_results_LogisticRegression', 'w') as f:
                json.dump(results, f,  cls=NumpyEncoder)    
    
    
