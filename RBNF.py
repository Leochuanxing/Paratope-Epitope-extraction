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
# Import the data
os.chdir("/home/leo/Documents/Database/Pipeline/Ready_2_2_0_1")
os.listdir()


with open('training_2_2_1_1_RBNF', 'r') as f:
    training = json.load(f)
with open('testing_2_2_1_1', 'r') as f:
    testing = json.load(f)
    
training[:36]
testing[:36]
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
# Define Radial Basis Functions
'''
Gaussian: 
    The Gaussian radial basis function
Inputs:
    x: the sample, in the form of [Ab_seq, Ag_seq, contact_number, total_contact], which is 
       in the same form of the elements in training or testing dataset
       
    center: in the same form as x
    
    radius:
        float, a real positive number, gives the scaling factor of the distance between x and the center
        
    distance_mode:
        a string, takes the value of either 'Multiplication_distance' or 'Addition_distance',
        Gives different way of combining the Ab_Ab distances and the Ag_Ag distances.        
'''

def Gaussian(x, center, radius, distance_mode = 'Multiplication_distance'):
    Ab_seq1 = To_seq(x[0])
    Ag_seq1 = To_seq(x[1])
    Ab_seq2 = To_seq(center[0])
    Ag_seq2 = To_seq(center[1])
    if distance_mode == 'Multiplication_distance':
        distance = Multiplication_distance(Ab_seq1, Ab_seq2, Ag_seq1, Ag_seq2)
    if distance_mode == 'Addition_distance':
        distance = Addition_distance(Ab_seq1, Ab_seq2, Ag_seq1, Ag_seq2)
    return math.exp(-distance**2/radius**2)
    
def Mrakov(x, center, radius, distance_mode = 'Multiplication_distance'):
    Ab_seq1 = To_seq(x[0])
    Ag_seq1 = To_seq(x[1])
    Ab_seq2 = To_seq(center[0])
    Ag_seq2 = To_seq(center[1])
    if distance_mode == 'Multiplication_distance':
        distance = Multiplication_distance(Ab_seq1, Ab_seq2, Ag_seq1, Ag_seq2)
    if distance_mode == 'Addition_distance':
        distance = Addition_distance(Ab_seq1, Ab_seq2, Ag_seq1, Ag_seq2)
    return math.exp(-distance/radius)
'''
beta is a positive number
'''
def Inverse_Multi_Quadric(x, center, c, beta, distance_mode = 'Multiplication_distance'):
    Ab_seq1 = To_seq(x[0])
    Ag_seq1 = To_seq(x[1])
    Ab_seq2 = To_seq(center[0])
    Ag_seq2 = To_seq(center[1])
    if distance_mode == 'Multiplication_distance':
        distance = Multiplication_distance(Ab_seq1, Ab_seq2, Ag_seq1, Ag_seq2)
    if distance_mode == 'Addition_distance':
        distance = Addition_distance(Ab_seq1, Ab_seq2, Ag_seq1, Ag_seq2)
    return (distance**2 + c**2)**(-beta)

def Design_matrix(Ab_Ag_paires, centers, basis_function = 'Markov', distance_mode = 'Multiplication_distance'):

    training = Ab_Ag_paires

    design_matrix = np.zeros((len(training), len(centers)))
    for (i, j), value in np.ndenumerate(design_matrix):
        if basis_function == 'Gaussian':
            design_matrix[i, j] = Gaussian(training[i], centers[j], radius=1,  distance_mode=distance_mode)
        if basis_function == 'Markov':
            design_matrix[i, j] = Mrakov(training[i], centers[j], radius=1,  distance_mode=distance_mode)
        if basis_function == 'Inverse_Multi_Quadric':
            design_matrix[i, j] = Inverse_Multi_Quadric(training[i], centers[j], c=1, beta=2, distance_mode=distance_mode)
    
    return design_matrix

##############################################################################
###############################################################################
class RBFN(object):
    def __init__(self, training, testing, distance_mode = 'Multiplication_distance'):
        self.training = training
        self.testing = testing
        self.distance_mode = distance_mode
        self.params = {}
        observed_contact = []
        for parepi in training:
            observed_contact.append(parepi[2])
        self.observed_contact = np.asanyarray(observed_contact)
        
    def Centers_from_NewDataAnalysis(self):
        with open('RBFN_centers_Ab_in_Ag', 'r') as f:
            centers = json.load(f)
        self.centers = centers
        self.params['W'] =  0.1 * np.random.randn(len(self.centers), 1)
            


    
    
    '''
    Pruning:
        to return the prunned centers
    Input: 
        training: the training set
        delta: the cover distance. If delta == 0, all the training data are selected as 
        centers after the repeatations has been deleted.
    Output: 
        the centers, in the same form as the training set
    '''
#    def Pruning(self, delta = 0):
#        
#        distance_mode = self.distance_mode
#        if distance_mode == 'Multiplication_distance':
#            distance_function = Multiplication_distance
#            self.distance_function = distance_function
#        elif distance_mode == 'Addition_distance':
#            distance_function = Addition_distance
#            self.distance_function = distance_function
#        else:
#            return print('The distance mode is not proper.')
#        
#        clustered = {}
#        for parepi in self.training:
#            if (tuple(parepi[0]), tuple(parepi[1])) not in clustered:
#                clustered[(tuple(parepi[0]), tuple(parepi[1]))] = [parepi]
#            else:
#                clustered[(tuple(parepi[0]), tuple(parepi[1]))].append(parepi)
#        
#        training_set = []
#        for key, value in clustered.items():
#            training_set.append(value[0])
#            
#
#        if delta == 0:
#            self.centers = training_set
#            
#        else:
#            centers = []
#            while training_set != []:
#                centers.append(training_set [0])
#                new_center = training_set[0]
#                training_set.remove(training_set[0])  
#                              
#                if training_set != []:                    
#                    coverage = []
#                    for parepi in training_set:
#                        if distance_function(To_seq(new_center[0]), To_seq(parepi[0]), To_seq(new_center[1]), To_seq(parepi[1])) <= delta:
#                            coverage.append(parepi)
#                    for covered in coverage:
#                        training_set.remove(covered)                
#                
#            self.centers = centers
#            
#        self.params['W'] =  0.1 * np.random.randn(len(self.centers), 1)
        
        
        


        
        
        

    def Formula_ridge_find_reg(self, basis_function = 'Markov', distance_mode = 'Multiplication_distance'):
        design_matrix = Design_matrix(self.training, self.centers, basis_function , distance_mode)
        centers = self.centers
        REG_vector = [0.1*x for x in range(20)]
        GCV_vector = []
        GCV_lowest = 1000000
        reg_best = 0
        for reg in REG_vector:
            
            A = (design_matrix.T).dot(design_matrix) + np.eye(len(centers))*reg
            inverse_A = np.linalg.inv(A)
            W = inverse_A.dot((design_matrix.T)).dot(self.observed_contact)
            
            self.params['W'] = W
            
            P = np.eye(len(training)) - design_matrix.dot(inverse_A).dot(design_matrix.T)
            
            GCV = len(self.training) * ((self.observed_contact).T).dot(P.dot(P)).dot(self.observed_contact)
            
            GCV /= (np.trace(P))**2
            
            GCV_vector.append(GCV)
            
            if GCV <= GCV_lowest:
                GCV_lowest = GCV
                reg_best = reg
        
        self.GCV_lowest = GCV_lowest
        self.reg_best = reg_best
            
        plt.title('GCV Vs cut REG')
        plt.xlabel('REG')
        plt.ylabel('GCV')
        plt.plot(REG_vector, GCV_vector)
        plt.show()
        plt.close() 
        
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
                # conbine the test_Ag with disfferent positive centers Ab, and the predicted contact 
                # of the combination willcaltulated later
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
                                          basis_function = 'Markov', distance_mode = 'Multiplication_distance')
            # calculate the scores
            scores = design_matrix.dot(W.T)
            # indexed scores
            indexed_scores = []
            for i in range(len(scores)):
                indexed_scores.append([i, scores[i]])
            # sort the indexed scores to find the top n
            indexed_scores.sort(key = lambda x:x[1], reverse=True)
            #find the top_n prediction
            top_n_predictions = []
            for i in range(top_n):
                top_n_predictions.append(indexed_scores[1][0])
            # Check if the prediction is correct
            if predicted_center in top_n_predictions:
                correct_prediction += 1
                
        self.correct_rate =round(correct_prediction/len(testing), 2)


                
                
            

                
                
            
  
        
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

training[:20]
rbfn = RBFN(training, testing, distance_mode = 'Addition_distance')
rbfn.Centers_from_NewDataAnalysis()
rbfn.centers
#rbfn.Pruning(delta= -0.2)
len(rbfn.centers)
for i in [x*0.05 for x in range(21)]:
    rbfn.Pruning(delta = i)
    print(len(rbfn.centers))
rbfn.Design_matrix()
rbfn.design_matrix_train.shape
rbfn.Formula_ridge_find_reg()
rbfn.reg_best
rbfn.GCV_lowest
rbfn.Design_matrix(train_or_test = 'test', basis_function = 'Markov')
rbfn.Formula_ridge_predict_Ab()
keys = list(rbfn.predictions)
len(keys)
len(rbfn.testing)
rbfn.predictions[keys[0]]
keys[0]
clustered = {}
for parepi in rbfn.training:
    if (tuple(parepi[0]), tuple(parepi[1])) not in clustered:
        clustered[(tuple(parepi[0]), tuple(parepi[1]))] = [parepi]
    else:
        clustered[(tuple(parepi[0]), tuple(parepi[1]))].append(parepi)

rbfn.design_matrix_test.shape
rbfn.params['W'].shape
(rbfn.design_matrix_test.dot(rbfn.params['W'].T)).shape
