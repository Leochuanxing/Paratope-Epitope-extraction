import math
import numpy as np
import json
import os
import copy
'''
Import file to generate cross testing data and cross training data
'''
os.chdir('/home/leo/Documents/Database/Pipeline/Codes and testing data')
from RBFN_coverage import Generate_cross_testing_training
from Logistic_Regression import Truncated, Train_BFGS, Prediction_area,  Multiplication_similarity, To_seq

def main(wd_results, wd_negative_samples, header):
    # Indicate the process
    print('working on ' + header)
    
    os.chdir(wd_results)
    with open(header+'_aa_train', 'r') as f:
        positive_training_set = json.load(f)

        
    os.chdir(wd_negative_samples)
    with open(header+'_train_negative', 'r') as f:
        negative_training_set = json.load(f)

    # Calculate the big design matrix so that it can be used to extract smaller matrix
    # Combine the training set
    training_set = copy.deepcopy(positive_training_set)
    training_set.extend(negative_training_set)
    # Indicate the process
    print('Calculating the big multiplication matrix of ' + header + '.')
        print('Number of total training samples' , len(training_set))
    # Creata an empty matrix as a container
    big_multiplication_matrix = np.zeros((len(training_set), len(training_set)))
    n_row = np.shape(big_multiplication_matrix)[0]
    n_col = n_row
    # Load up the matrix
    for i in range(n_row):
        for j in range(i-1, n_col):
            Ab_seq1=To_seq(training_set[i][0])
            Ab_seq2 = To_seq(training_set[j][0])
            Ag_seq1 = To_seq(training_set[i][1])
            Ag_seq2 = To_seq(training_set[j][1])
            big_multiplication_matrix[i, j] = Multiplication_similarity(Ab_seq1, Ab_seq2, Ag_seq1, Ag_seq2)
            big_multiplication_matrix[j, i] = big_multiplication_matrix[i,j]
    
    # We use the Generate_cross_testing_training to generate the cross testing and cross training set
    cross_positive_testing_pos, cross_positive_training_pos =\
    Generate_cross_testing_training(positive_training_set, 6)
    
    cross_negative_testing_pos, cross_negative_training_pos =\
    Generate_cross_testing_training(negative_training_set, 6)
    
    truncate_percentage = [0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
    results = {}
    results['precentages'] = truncate_percentage
    results['areas'] = []
    # Iterate over truncate_percentage and the above cross validation sets
    for percentage in truncate_percentage:
        area = 0
        for i in range(6):
            sub_positive_training_pos = cross_positive_training_pos[i]
            sub_positive_testing_pos = cross_positive_testing_pos[i]
            sub_negative_training_pos = cross_negative_training_pos[i]
            sub_negative_testing_pos = cross_negative_testing_pos[i]
            # Extract the training_matrix, positive_testing_matring, negative_testing_matrix
            adjusted_sub_negative_training_pos = (np.array(sub_negative_training_pos) + len(positive_training_set)).tolist()
            adjusted_sub_negative_testing_pos = (np.array(sub_negative_testing_pos) + len(positive_training_set)).tolist()
            
            sub_training_pos = copy.deepcopy(sub_positive_training_pos)
            sub_training_pos.extend(adjusted_sub_negative_training_pos)
            
            sub_training_matrix = big_multiplication_matrix[sub_training_pos, :][:,sub_training_pos]
            sub_positive_testing_matrix = big_multiplication_matrix[sub_positive_testing_pos, :][:, sub_training_pos]
            sub_negative_testing_matrix = big_multiplication_matrix[adjusted_sub_negative_testing_pos,:][:, sub_training_pos]
            
            # Create the training_indicator vector
            p = len(sub_positive_training_pos)
            n = len(sub_negative_training_pos)
            training_indicator = np.vstack((np.ones((p, 1)), np.zeros((n, 1))))

#            # Load the truncate_package
#            truncate_package = {}
#            truncate_package['p'] = len(sub_positive_training_pos)
#            truncate_package['n'] = len(sub_negative_training_pos)
#            truncate_package['training_matrix'] =sub_training_matrix
#            truncate_package['positive_testing_matrix'] = sub_positive_testing_matrix
#            truncate_package['negative_testing_matrix'] = sub_negative_testing_matrix
            # Truncate the model
            number = math.floor((n+p) * percentage)
            truncated_training_sum = Truncated(sub_training_matrix , p, n, number)
            truncated_positive_testing_sum = Truncated(sub_positive_testing_matrix, p, n, number)
            truncated_negative_testing_sum = Truncated(sub_negative_testing_matrix, p, n, number)
            # Train the model
            parameter_BFGS, loss_BFGS = Train_BFGS(truncated_training_sum, training_indicator, 
                                                   reg=1, rho=0.85, c=1e-2, termination=1e-3)
            # Prepare the testing data
            truncated_testing_sum = np.vstack((truncated_positive_testing_sum, truncated_negative_testing_sum))
            positive_indicator = np.ones((len(sub_positive_testing_pos), 1))
            negative_indicator = np.zeros((len(adjusted_sub_negative_testing_pos), 1))
            testing_indicator = np.vstack((positive_indicator, negative_indicator))
            print('Number of testing samples ', len(testing_indicator))
            # Calculate the area under the precision_recall curve
            area += Prediction_area(truncated_testing_sum, testing_indicator, parameter_BFGS) 
            
            
        # Store the percentage and the corresponding area
        results['areas'].append(area)
        print('area ' , area)
    return results
###################################################################################
'''
Order the working files
'''
def Order_files():
    ordered_headers = []
    for diff in range(4):
        for i in range(2,5):
            for j in range(2, 5):
                if abs(i-j) == diff:
                     header = str(i)+'_'+str(j)
                     ordered_headers.append(header)
                     print(header)
    return ordered_headers
#ordered_headers = Order_files()
#ordered_headers
ordered_headers = ['2_2','3_3', '2_3', '3_2', '3_4', '4_3','2_4', '4_2','4_4',\
                   '1_1','1_2', '2_1', '1_3', '3_1', '1_4', '4_1']

#len(ordered_headers)
###################################################################################
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
            results =  main(wd_results, wd_negative_samples, header)
            
             # Save the results
            os.chdir(wd_results)
            with open(header+'_cross_results_LogisticRegression', 'w') as f:
                json.dump(results, f,  cls=NumpyEncoder)
 ##################################################################################
os.chdir('/home/leo/Documents/Database/Pipeline/Results/1_free')
with open('2_2_test_results', 'r') as f:
    test_results = json.load(f)
test_results.keys()
#test_results['0']
#with open('2_2_aa_train', 'r') as f:
#    training_set = json.load(f)      
#len(training_set)
