'''
This file is to find the complexes with measured affinities and different from each other by a few mutations
'''
import json
import os
import copy
########################################################################
'''
Define an aligner
'''        
from Bio import Align
aligner = Align.PairwiseAligner()
aligner.open_gap_score = -1 # it is unlikly that two sequence of the same source could be different by a gap
aligner.extend_gap_score = -1
aligner.match = 1
aligner.mismatch = -1 
aligner.mode = 'global'
####################################################################
'''
Parse_the_summary:
    a function to parse the summary file.
    input:
    output:
        ids_affinity: a list with the pdbids matched up with their affinities
'''
def Parse_the_summary():
    os.chdir('/home/leo/Documents/Database/Pipeline/Complexes with Affinity')
    os.listdir()
    summary = open('summary.tsv', 'r')
    file = summary.readlines()        
    summary.close
    
    parsed = []
    for string in file:
       parsed.append(string.split('\t')) 
     
    for i in range(len(parsed[0])):
        if parsed[0][i] == 'affinity':
            affinity_pos = i
    
    ids_affinity = []
    for variables in parsed[1:]:
        if len(variables)>= affinity_pos:
            if [variables[0], variables[affinity_pos]] not in ids_affinity:    
                ids_affinity.append([variables[0], variables[affinity_pos]])
    
    return ids_affinity

    

'''
Define a function to find all the complexes different with each other by a few mutations
'''   

#################################################################################

           

def To_seq( aa_sequence):
    
    TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
   ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
   ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]
    
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    
    seq_obj = None

    seq_single_letter = ''
    for aa in aa_sequence:
        for TS in TripleSingle:
            if TS[0] == aa:
                seq_single_letter += TS[1]
    seq_obj = Seq(seq_single_letter, IUPAC.protein)
    
    return seq_obj
    '''
    SeqCDR:
        to convert the sequence into sequence object and be prepared for futher alignment
    Input:
        sequence:
            a dictionary with keys pdbids and the values a dictionary with keys the names of the 
            chains.
        good_matched_ids:
    Output:
        in the form of a list [['1abvhH', Seq_obj], ['1abvlL', Seq_obj]]
        ['1abvhH', Seq_obj] means the pdbid is 1abv heavy chain with chain name H, and 
        the Seq_obj is the sequence object of the concatenated CDR1, CDR2, CDR3 of chain H.
    '''   
    #We should align more than the CDR in order to rule out the confounding effect.
def SeqCDR(sequence, good_matched_ids):
#    l_range = [[20, 43], [46, 66], [86, 110]]
#    h_range = [[22, 40], [47, 74], [96, 132]]
    # Concatenate the sequences 
    # All the seqed cdr sequences will be stored in a dictionary
    # like seq_dict = {}
    seq_dict = {}
    for pdbid in good_matched_ids:
        seq_dict[pdbid] = []
        # Take the first combination of the antibody and antigen chains from the complex
        combination = good_matched_ids[pdbid][0]
        if combination[2] != '':
            if combination[0] != '':
                l = len(To_seq(sequence[pdbid][combination[0]]))
                CDRH = [pdbid+'h'+combination[0]]
                CDRH.append(To_seq(sequence[pdbid][combination[0]][22:40]))
                CDRH.append(To_seq(sequence[pdbid][combination[0]][47:74]))
                CDRH.append(To_seq(sequence[pdbid][combination[0]][96:min(132,l)]))
                seq_dict[pdbid].append(CDRH)
            else:
                seq_dict[pdbid].append([])
            if combination[1] != '':
                l = len(sequence[pdbid][combination[1]])
                CDRL = [pdbid+'l'+combination[1]]
                CDRL.append(To_seq(sequence[pdbid][combination[1]][20:43]))
                CDRL.append(To_seq(sequence[pdbid][combination[1]][46:66]))
                CDRL.append(To_seq(sequence[pdbid][combination[1]][86:min(110, l)]))
                seq_dict[pdbid].append(CDRL) 
            else:
                seq_dict[pdbid].append([])
                       
    return seq_dict
        
#seq_dict = SeqCDR(sequence, good_matched_ids)

'''
Hamming_like_dist:
    A function to calculate the difference between the CDRs of two given complex
    
    input:
        u: one value of the seq_dict
        v: the other value of the seq_dict
'''
    
def CDR_diff (u, v):
    # 1 creat the similarity matrix
    # 2 Do hyerachical clustering
    # 3 Chose the representative
    dist = 0 
    l = 0
    from Bio import Align
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -1 # it is unlikly that two sequence of the same source could be different by a gap
    aligner.extend_gap_score = -1
    aligner.match = 1
    aligner.mismatch = -1 
    aligner.mode = 'global'
    
    H_u = u[0]
    L_u = u[1]
    try: 
        H_v = v[0]
    except:
        print(v)
    L_v = v[1]
    for i in range(1,4):
        if H_v != [] and H_u != []:
            dist += aligner.score(H_u[i], H_v[i])
            l += min(len(H_u[i]), len(H_v[i]))
        if L_v != [] and L_u != []:
            dist += aligner.score(L_u[i], L_v[i])
            l  += min(len(L_u[i]), len(L_v[i]))
        
    return l - dist # This valuse reflect how many differences are there in the sequences


#import numpy as np
'''
Paires_CDR:
    a function to get the paires with the mismatches defined in CDR_diff
    equals to cut
Inuput:
    seq_dict:
        a dictionary, the return of SeqCDR
    cut:
        an integer, gives the number of 'mismatches' defined by the Hamming_like_dist
Output:
    paires:
        a list, gives paires with the mismatch of the CDRs equal to the cut.
'''
def Paires_CDR(seq_dict, cut):
    keys = list(seq_dict.keys())
#    dist_matrix = np.zeros((len(keys), len(keys))) 
    paires = []
    for i in range(len(keys)):
        for j in range(i, len(keys)):
            distance = CDR_diff(seq_dict[keys[i]], seq_dict[keys[j]])
#            dist_matrix[i,j] = distance
            if distance == cut and i != j:
                paires.append([keys[i], keys[j]])
    return paires
    
#paires = Get_Paires(seq_dict, cut =0)
#paires

#########################################################################
'''
Chain_diff:
    input: 
        pdbid1, chain1, pdbid2, chain2
        all the inputs are strings, pdbid is in the form of '1adq'
        chain is in the form of 'H'
        
        sequence:
            a dictionary, contains all the amino acid sequence information for any given 
            pdbid and any chain name.
    output:
        the difference score
'''
def Chain_diff(pdbid1, chain1, pdbid2, chain2, sequence):
    
    align_score = aligner.score(To_seq(sequence[pdbid1][chain1]),\
                                To_seq(sequence[pdbid2][chain2]))
    
    possible_max = max(len(sequence[pdbid1][chain1]), len(sequence[pdbid2][chain2]))
    
    return possible_max - align_score

########################################################################

##########################################################################
'''

'''

#########################################################################
'''
Paires_Ab_fixed:
    This function is to get the paires with Ab_fixed and the Ag chains different within
    the given value
    
    Input:
        sequence, good_matched_ids
        
        Ag_cut: an integer, gives the difference vaule: length of Ag - alignment_score
    output:
        paires_Ab_fixed:
            all the pairs with the same Ab but different Ag chains, and the difference of 
            the Ag chain is within the given value Ag_cut
'''
def Paires_Ab_fixed(sequence, good_matched_ids, Ag_cut):
    # Creat an empty list to contain the results.
    paires_Ab_fixed = []
    
    keys = list(good_matched_ids.keys())
    
    for i in range(len(keys)):
        for j in range(i, len(keys)):
            # Take out the name of the pdbid
            pdbid1 = keys[i]
            pdbid2= keys[j]
            # Take out the first element of the value corresponding to 
            # pdbid1 and pdbid2 in the dictionary of good_matched_ids
            combination_1 = good_matched_ids[pdbid1][0]
            combination_2 = good_matched_ids[pdbid2][0]
            # Take out the name of the heavy chain
            H1=combination_1[0]
            H2 = combination_2[0]
            # Take out the name of the light chain 
            L1 = combination_1[1]
            L2 = combination_2[1]
            # Take out the naem of the antigen chain
            A1= combination_1[2]
            A2 = combination_2[2]
            
            # Align the heavy chains
            Ab_diff_H = 0
            if H1 != '' and H2 != '':
                Ab_diff_H =  Chain_diff(pdbid1, H1, pdbid2, H2, sequence)
            # Align the light chains
            Ab_diff_L = 0
            if L1 != '' and L2 != '':
                Ab_diff_L = Chain_diff(pdbid1, L1, pdbid2, L2, sequence)
            # Align the Ag chains
            Ag_diff = 0
            if A1 != '' and A2 != '':
                Ag_diff = Chain_diff(pdbid1, A1, pdbid2, A2, sequence)
            
            # Check the requirement that the Ab chains should be fixed and the Ag
            # chains should be different within certain range.
            if Ab_diff_H + Ab_diff_L == 0 and Ag_diff <= Ag_cut and Ag_diff > 0:
                paires_Ab_fixed.append([keys[i], keys[j]])
            
    return paires_Ab_fixed
############################################################################            
'''
Affinity_control:
    this function is to select from a list of given paires, those paires that
    the fold change of the affinity is no less than a given value 'control'
    
    Input:
        paires:
            a list of pdbid paires
        ids_affinity:
            the return value of Parse_the_summary
        control:
            a float, gives the threshold of the affinity change
    Output:
        workable_list:
            a list gives the pairs satisfy the affinity_control, with the correspnonding
            affinity values taged along. The elements are in the form shown below:
            [['1hh9', '1hh6'], ['1.00E-05', '1.00E-07']]
'''            
def Affinity_control(paires, ids_affinity, control):
    workable_list = []
    # Go into one match
    for match in paires:
        # Creat a list of length 2 to contain the affinities
        # Find the corresponding affinities.
        affinity_match = [0, 0]
        for i in ids_affinity:
            if i[0] == match[0]:
                affinity_match[0] = i[1]
            if i[0] == match[1]:
                affinity_match[1] = i[1]
        # Calculate the corresponding affinity changes      
        change = float(affinity_match[0])/float(affinity_match[1])
        # Select according to the control
        if max(change, 1/change) > control:
            workable_list.append([match, affinity_match])
    
    return workable_list
#######################################################################
'''
Here is the working flow of finding the paires with the Ab fixed and the difference
of the Ag chains are within a given range.
'''
os.chdir('/home/leo/Documents/Database/Pipeline/Complexes with Affinity')
with open('sequence', 'r') as f:
    sequence = json.load(f)
with open('good_matched_ids', 'r') as f:
    good_matched_ids = json.load(f)
with open('contact', 'r') as f:
    contact = json.load(f)

ids_affinity = Parse_the_summary()
paires_Ab_fixed = Paires_Ab_fixed(sequence, good_matched_ids,Ag_cut=6)
workable_list = Affinity_control(paires_Ab_fixed, ids_affinity, control = 1)
paires_Ab_fixed
workable_list
######################################################################

pdbid1 = '3a6b'
pdbid2 = '3a67'
combination_1 = good_matched_ids[pdbid1][0]
combination_2 = good_matched_ids[pdbid2][0]
# Take out the name of the heavy chain
H1=combination_1[0]
H2 = combination_2[0]
# Take out the name of the light chain 
L1 = combination_1[1]
L2 = combination_2[1]
# Take out the naem of the antigen chain
A1= combination_1[2]
A2 = combination_2[2]

# Align the heavy chains
Ab_diff_H = 0
if H1 != '' and H2 != '':
    Ab_diff_H =  Chain_diff(pdbid1, H1, pdbid2, H2, sequence)
# Align the light chains
Ab_diff_L = 0
if L1 != '' and L2 != '':
    Ab_diff_L = Chain_diff(pdbid1, L1, pdbid2, L2, sequence)
# Align the Ag chains
Ag_diff = 0
if A1 != '' and A2 != '':
    Ag_diff = Chain_diff(pdbid1, A1, pdbid2, A2, sequence)
Ab_diff_H
Ab_diff_L
Ag_diff
alignments=aligner.align(To_seq(sequence[pdbid1][H1]), To_seq(sequence[pdbid2][H2]))
for alignment in alignments:
    print(alignment)
good_matched_ids[pdbid1]
good_matched_ids[pdbid2]
########################################################################
'''
Paires_Ag_fixed:
    A function to find the paires with the antigen chain fixed and the difference 
    of the CDRs should be within a given range. Here we have to make sure that all 
    the difference of the antibody occur in the CDR regions
input:
    sequence, good_matched_ids, Ab_cut: they are the same as the input of the function
    Paires_Ab_fixed
    
    seq_dict: the returned value of SeqCDR
ouput:
    paires_Ab_fixed, a list of paires meet the requirement that the Ag chain are the 
    same, and the Ab chains are different with the given range.    
'''
def Paires_Ag_fixed(sequence, good_matched_ids, seq_dict, Ab_cut):
    paires_Ag_fixed = []
    
    keys = list(good_matched_ids.keys())
    
    for i in range(len(keys)):
        for j in range(i, len(keys)):
            # Take out the name of the pdbid
            pdbid1 = keys[i]
            pdbid2= keys[j]
            # Take out the first element of the value corresponding to 
            # pdbid1 and pdbid2 in the dictionary of good_matched_ids
            combination_1 = good_matched_ids[pdbid1][0]
            combination_2 = good_matched_ids[pdbid2][0]
            # Take out the name of the heavy chain
            H1=combination_1[0]
            H2 = combination_2[0]
            # Take out the name of the light chain 
            L1 = combination_1[1]
            L2 = combination_2[1]
            # Take out the naem of the antigen chain
            A1= combination_1[2]
            A2 = combination_2[2]
            
            # Align the heavy chains
            Ab_diff_H = 0
            if H1 != '' and H2 != '':
                Ab_diff_H =  Chain_diff(pdbid1, H1, pdbid2, H2, sequence)
            # Align the light chains
            Ab_diff_L = 0
            if L1 != '' and L2 != '':
                Ab_diff_L = Chain_diff(pdbid1, L1, pdbid2, L2, sequence)
            # Align the Ag chains
            Ag_diff = 0
            if A1 != '' and A2 != '':
                Ag_diff = Chain_diff(pdbid1, A1, pdbid2, A2, sequence)
            # Align the CDRs
            CDRs_diff = 0
            if seq_dict[pdbid1] != [] and seq_dict[pdbid2] != []:
                CDRs_diff = CDR_diff(seq_dict[pdbid1], seq_dict[pdbid2])  
                     
            # We have to make sure the antigen chain is fixed and the differces of the 
            # Ab chains are within the given range and thoes differences are located 
            # in the CDR regions.
            if Ag_diff == 0 and Ab_diff_H + Ab_diff_L <= CDRs_diff and\
            CDRs_diff <= Ab_cut and CDRs_diff != 0:
                paires_Ag_fixed.append([keys[i], keys[j]])
                
    return paires_Ag_fixed
                
#######################################################################
'''
Here is the workflow while the Ag is fixed and the difference of CDRs is 
within the a given threshhold
''' 
os.chdir('/home/leo/Documents/Database/Pipeline/Complexes with Affinity')
with open('sequence', 'r') as f:
    sequence = json.load(f)
with open('good_matched_ids', 'r') as f:
    good_matched_ids = json.load(f)
with open('contact', 'r') as f:
    contact = json.load(f)

ids_affinity = Parse_the_summary()
seq_dict = SeqCDR(sequence, good_matched_ids)
paires_Ag_fixed = Paires_Ag_fixed(sequence, good_matched_ids, seq_dict, Ab_cut=12)
workable_list_Ag_fixed = Affinity_control(paires_Ag_fixed, ids_affinity, control = 1)
paires_Ag_fixed
workable_list_Ag_fixed
#######################################################################
'''
The following is a step by step check
'''
pdbid1 = '2nyy'
pdbid2= '2nz9'
# Take out the first element of the value corresponding to 
# pdbid1 and pdbid2 in the dictionary of good_matched_ids
combination_1 = good_matched_ids[pdbid1][0]
combination_2 = good_matched_ids[pdbid2][0]
# Take out the name of the heavy chain
H1=combination_1[0]
H2 = combination_2[0]
# Take out the name of the light chain 
L1 = combination_1[1]
L2 = combination_2[1]
# Take out the naem of the antigen chain
A1= combination_1[2]
A2 = combination_2[2]

# Align the heavy chains
Ab_diff_H = 0
if H1 != '' and H2 != '':
    Ab_diff_H =  Chain_diff(pdbid1, H1, pdbid2, H2, sequence)
# Align the light chains
Ab_diff_L = 0
if L1 != '' and L2 != '':
    Ab_diff_L = Chain_diff(pdbid1, L1, pdbid2, L2, sequence)
# Align the Ag chains
Ag_diff = 0
if A1 != '' and A2 != '':
    Ag_diff = Chain_diff(pdbid1, A1, pdbid2, A2, sequence)
# Align the CDRs
CDRs_diff = 0
if seq_dict[pdbid1] != [] and seq_dict[pdbid2] != []:
    CDRs_diff = CDR_diff(seq_dict[pdbid1], seq_dict[pdbid2]) 

Ab_diff_H
Ab_diff_L
Ag_diff
CDRs_diff
len(sequence[pdbid1][A1])

if H1 != '':
    alignments_H = aligner.align(To_seq(sequence[pdbid1][H1]), To_seq(sequence[pdbid2][H2]))
    for alignment in alignments_H:
        print(alignment)
 
if L1 != '':
    alignments_L = aligner.align(To_seq(sequence[pdbid1][L1]), To_seq(sequence[pdbid2][L2]))
    for alignment in alignments_L:
        print(alignment)

seq1 = seq_dict[pdbid1]
seq2 = seq_dict[pdbid2]
for i in range(1, 4):    
    if seq1[0] != []:
        for alignment in  aligner.align(seq1[0][i], seq2[0][i]):
            print(alignment, '      Heavy Chain')
    if seq1[1] != []:
        for alignment in  aligner.align(seq1[1][i], seq2[1][i]):
            print(alignment, '      Light Chain')

for i in range(len(sequence[pdbid1][H1])):
    if sequence[pdbid1][H1][i] != sequence[pdbid2][H2][i]:
        print(i)
pdbid1
pdbid2
sequence[pdbid1][H1][28:31]
sequence[pdbid2][H2][28:31]

good_matched_ids[pdbid1]
good_matched_ids[pdbid2]
#####################################################################

'''
The structure of the workable_dicts:
    1. workable_dicts is a list of dictionaries
    2. each element in the workable_dicts is a dictionary with keys mutations and affinities
       the value of mutations: is a list, for example
       [['3eys', 'Q', [3], ['GLU', 'ALA'], ['H', 'L']],
    ['3eyu', 'Q', [4, 5], ['ASP'], ['H', 'L']],
   ['3eyu', 'Q', [0], [''], ['H', 'L']]]
        This list contains three elements, each element is in the form 
        ['3eys', 'Q', [3], ['GLU', 'ALA'], ['H', 'L']]
        where '3eys' gives the pdbid
               'Q' gives the mutated chain
               [3] gives the position of the mutated amino acids in the mutated chain
               ['GLU', 'ALA'] gives the amino acids mutated into
               ['H', 'L'] gives the opposite chains of the mutated chians
               
        the value of 'affinities' is in the form of 
        [['3eys', '3.00E-06'], ['2r0z', '3.40E-06']]
        with two elements represent two conterpart complexes.
        '3eys' is mutated into '2r0z' or '2r0z' is mutated into '3eys'
        ['3eys', '3.00E-06'] says the affinity of '3eys' is 3.00E-06. Here the affinity
        is measured in Kd.
'''

workable_dicts_Ab_fixed  = [{'mutations': [['1hh9', 'C', [7, 9], ['GLY', 'ARG'], ['A', 'B']],
    ['1hh6', 'C', [7, 9], ['ASN', 'LYS'], ['A', 'B']]],
  'affinities': [['1hh9', '1.00E-05'], ['1hh6', '1.00E-07']]},
                 
 {'mutations': [['1nbz', 'C', [95, 96], ['ALA', 'LYS'], ['A', 'B']],
    ['1nby', 'C', [95, 96], ['LYS', 'ALA'], ['A', 'B']]],
  'affinities': [['1nbz', '1.06E-06'], ['1nby', '9.09E-05']]},
                
 {'mutations': [['1nbz', 'C', [96], ['LYS'], ['A', 'B']],
    ['1dqj', 'C', [96], ['ALA'], ['A', 'B']]],
  'affinities': [['1nbz', '1.06E-06'], ['1dqj', '2.86E-09']]},
                
 {'mutations': [['1nby', 'C', [95], ['LYS'], ['A', 'B']],
    ['1dqj', 'C', [95], ['ALA'], ['A', 'B']]],
  'affinities': [['1nby', '9.09E-05'], ['1dqj', '2.86E-09']]},
                
  {'mutations': [['2fx8', 'P', [3], ['LYS'], ['H', 'L']],
                 ['2fx8', 'P', [6], ['TRP'], ['H', 'L']],
                 ['2fx8', 'P', [9,10], ['LYS', 'LYS', 'LYS', 'LYS'], ['H', 'L']],
    ['2fx9', 'P', [3], ['ASP'], ['H', 'L']],
    ['2fx9', 'P', [6], ['ASN'], ['H', 'L']],
    ['2fx9', 'P', [9,10,11,12], ['ARG', 'ARG'], ['H', 'L']]],
  'affinities': [['2fx8', '3.02E-07'], ['2fx9', '1.70E-08']]}]
                 
                   
workable_dicts_Ag_fixed =[ {'mutations': [['3a6b', 'L', [30,31], ['ASP','ASN'], ['Y']],
                   ['3a67', 'L', [30,31], ['ASN', 'ASP'], ['Y']]],             
  'affinities': [['3a6b', '1.08E-07'], ['3a67', '5.62E-09']]},
 
  {'mutations': [['3a6c', 'L', [31], ['ASP'], ['Y']],
                   ['3a6c', 'L', [91], ['ASN'], ['Y']],
                   ['3a6b', 'L', [31], ['ASN'], ['Y']],
                   ['3a6b', 'L', [91], ['ASP'], ['Y']]],             
  'affinities': [['3a6c', '3a6b'], ['7.14E-09', '1.08E-07']]},
                 
  {'mutations': [['3a6c', 'L', [30], ['ASP'], ['Y']],
                   ['3a6c', 'L', [91], ['ASN'], ['Y']],
                   ['3a67', 'L', [30], ['ASN'], ['Y']],
                   ['3a67', 'L', [91], ['ASP'], ['Y']]],             
  'affinities': [['3a6c', '3a67'], ['7.14E-09', '5.62E-09']]},
                 
  {'mutations': [['2nyy', 'D', [28,29,30], ['SER', 'ASP', 'HIS'], ['A']],

                  ['2nz9', 'D', [28,29,30],['LYS', 'TYR', 'ASP'], ['A']]],             
  'affinities': [['3a6c', '3a67'], ['7.14E-09', '5.62E-09']]}]

os.chdir('/home/leo/Documents/Database/Pipeline/Complexes with Affinity')
with open('workable_dicts_Ab_fixed', 'w') as f:
    json.dump(workable_dicts_Ab_fixed, f)
    
with open('workable_dicts_Ag_fixed', 'w') as f:
    json.dump(workable_dicts_Ag_fixed, f)

