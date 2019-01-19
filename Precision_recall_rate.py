# Import all the function in RBFN_coverage for latter usage.
import json
import os
import numpy as np
from matplotlib import pyplot as plt
os.chdir('/home/leo/Documents/Database/Pipeline/Codes and testing data')
from RBFN_coverage import Pre_processing
from RBFN_coverage import Cross_coverage_centers
from RBFN_coverage import Get_cross_coeff
from RBFN_coverage import Generate_testing_design_matrix


######################################################################
# We want to calculate the precision and recall rate for different set of data
# Let's try one set first, then, we can do batch processing.
# Try to take look at the 2_2 in the 1_free file
#os.chdir('/home/leo/Documents/Database/Pipeline/Results/1_free')
#with open('2_2_cross_results', 'r') as f:
#    results = json.load(f)
#with open('2_2_aa_train', 'r') as f:
#    positive_training_set = json.load(f)
#with open('2_2_aa_test', 'r') as f:
#    positive_testing_set = json.load(f)
#    
#os.chdir('/home/leo/Documents/Database/Pipeline/Negative samples/1_free')
#with open('2_2_train_negative', 'r') as f:
#    negative_training_set = json.load(f)
#with open('2_2_test_negative','r') as f:
#    negative_testing_set = json.load(f)
#len(positive_training_set)
#results.keys()    
#results['type']
#results['best_hyperparameter']
#results['0']['areas']
###########################################################################
'''
Use the Pre_processing in the RBFN_coverage to create the data_dict
'''
    
#data_dict= Pre_processing(positive_training_set, positive_testing_set,\
#                              negative_training_set, negative_testing_set)

#data_dict.keys()
#len(data_dict['training_set'])
#len(positive_training_set)
#len(positive_testing_set)
#len(data_dict['testing_set'])
#len(data_dict['observed_values'])
#data_dict['observed_values'][:6]
#data_dict['observed_values'][-6:]
#data_dict['observed_values'][2593:2598]
#data_dict['design_matrix'].shape
#data_dict['distance_matrix'].shape
###########################################################################
# We can use the Cross_coverage_ceters function to calculate the centers, the only
# difference is that here we only have the best percentage

'''
We reuse the Cross_coverage_centers function
Input:
    Replace the cross_data_dict with the data dict
    cross_parameter should be replaced with parameter, with 
    parameter['center_percentages'] = best_hyperparameter
Output:
    parameter['centers'] gives the indices of the selected centers
''' 
#results['best_hyperparameter']
#parameter = {}
#parameter['center_percentages'] = [results['best_hyperparameter']]# pay attention to the form of the data
#parameter = Cross_coverage_centers(data_dict, parameter)

#parameter.keys()
#len(parameter['centers'][0])
#parameter['centers'][0][:6]
#parameter['centers'][0][-6:]
###########################################################################
'''
Use the Get_cross_coeff to train the model and get the best coeff
'''

#parameter = Get_cross_coeff(cross_parameter=parameter, cross_data_dict=data_dict)
#len(parameter['all_coeff'][0])
#len(parameter['all_coeff'][0]['coeff'])
############################################################################
'''
Precision_recall_rate:
    Calculate the precision and recall relation
Input:
    parameter
    data_dict
Output:
    precision_recall:
        a dictionary, contains:
            precision_recall['precision']: a list of precision rate
            precision_recall['recall']: a list of recall rate
            *****the above two list should be corresponding to each other*****
'''
#data_dict.keys()
#parameter.keys()
def Precision_recall_rate(data_dict, parameter):
    # take out the values
    training_set = data_dict['training_set']
    testing_set = data_dict['testing_set']
    centers = parameter['centers'][0]
    coeff = np.array(parameter['all_coeff'][0]['coeff']).reshape((-1,1))
    # Copy some code from Hyperparameter_evaluation
    testing_design_matrix = Generate_testing_design_matrix(testing_set, centers, training_set,\
                                       mode = 'Multiplication', basis_function = 'Markov')

    prediction = testing_design_matrix.dot(coeff)
        
    # Related the observed testing values with the predicted values
    pred_observed = []
    for i in range(len(testing_set)):
        pred_observed.append([prediction[i,0], testing_set[i][2]])
    pred_observed.sort(key=lambda x:x[0], reverse=True)
        
    # Calculate the recall and precision
    precision=[]
    recall = []
    denominator = 0.5 * len(testing_set)
    n_positive = 0
    n_negative = 0
    for i in pred_observed:
        if i[1] == 1:
            n_positive += 1
        else:
            n_negative += 1

        precision.append(n_positive/(n_positive+n_negative))
        recall.append(n_positive/denominator)
    # Load the results
    precision_recall = {}
    precision_recall['precision']=precision
    precision_recall['recall']=recall
    
    return precision_recall

#precision_recall = Precision_recall_rate(data_dict, parameter)
#
#len(precision_recall['precision'])
#len(precision_recall['recall'])
###############################################################
'''
Plot_precision_recall:
    This function is to plot the precision recall rate
Input:
    precision_recall
Output:
    
'''
def Plot_precision_recall(precision_recall):
    plt.plot(precision_recall['recall'], precision_recall['precision'])
    return

#Plot_precision_recall(precision_recall)
##############################################################################
'''
Define a main function, with input the working directory
Input:
    wd_results: 
        the working directory of the Results directory
    wd_negative_samples: 
        the working directory of the negative samples
    header: 
        a string in the form of '2_2', '3_2', etc. It can be applied to the beginning
        of the file names
    
'''

def main(wd_results, wd_negative_samples, header):
    
    os.chdir(wd_results)
    with open(header+'_cross_results', 'r') as f:
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
    
    # Generate the data_dict
    data_dict= Pre_processing(positive_training_set, positive_testing_set,\
                              negative_training_set, negative_testing_set)
    # Generate the centers
    parameter = {}
    parameter['center_percentages'] = [results['best_hyperparameter']]# pay attention to the form of the data
    parameter = Cross_coverage_centers(data_dict, parameter)
    # Train the model
    parameter = Get_cross_coeff(cross_parameter=parameter, cross_data_dict=data_dict)
    # Calculate the precision and recall rate
    precision_recall = Precision_recall_rate(data_dict, parameter)
    
    # Load the result to a dictionary called result
    result = {}
    result['centers'] = parameter['centers'][0]
    result['precision'] = precision_recall['precision']
    result['recall'] = precision_recall['recall']
    result['pattern'] = header
    
    return result
###################################################################################
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
##############################################################################
# Do the batch processing as we do in the RBFN_coverage
if __name__ == '__main__':
    # Oder the sequence of the files
    ordered_headers = ['2_2','3_3', '2_3', '3_2', '3_4', '4_3','2_4', '4_2','4_4',\
                   '1_1','1_2', '2_1', '1_3', '3_1', '1_4', '4_1']
    
    # Match up the directories
    directory_paires = [['/home/leo/Documents/Database/Pipeline/Results/1_free',\
      '/home/leo/Documents/Database/Pipeline/Negative samples/1_free'],\
     ['/home/leo/Documents/Database/Pipeline/Results/0_free',\
                           '/home/leo/Documents/Database/Pipeline/Negative samples/0_free']]
    # Do the work
    for i in range(2):
        wd_results = directory_paires[i][0]
        wd_negative_samples = directory_paires[i][1]
        for header in ordered_headers:
            print('Working on '+ header+'  '+ wd_results[-6:])
            result =  main(wd_results, wd_negative_samples, header)            
             # Save the results
            os.chdir(wd_results)
            with open(header+'_test_results', 'w') as f:
                json.dump(result, f,  cls=NumpyEncoder)
    

#directory_paires = [['/home/leo/Documents/Database/Pipeline/Results/1_free',\
#  '/home/leo/Documents/Database/Pipeline/Negative samples/1_free'],\
# ['/home/leo/Documents/Database/Pipeline/Results/0_free',\
#                       '/home/leo/Documents/Database/Pipeline/Negative samples/0_free']]
#wd = directory_paires[0][0]
#wd
#wd[-6:]





















