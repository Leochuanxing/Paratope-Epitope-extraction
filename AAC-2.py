import os
os.getcwd()
os.chdir("/home/leo/Documents/Database/Pipeline/Homo")
os.listdir()
###############################################
import json
with open('good_combined_ids', 'r') as f:
    good_combined_ids = json.load(f)
with open('good_matched_ids', 'r') as f:
    good_matched_ids = json.load(f)
    
len(good_combined_ids)
len(good_matched_ids)
################################################
# Extract the chain sequences 
'''
Chain_seq is to extract sequences for all the chains in the complex
Imputs 
        file, a pdb file
        combined_chain_id, a list of the form ['BDF', 'ACE', 'GH']
         in the order of heavy chains, light chains, and antigen chains
Returns: seq, a dictionary of sequences, with the chain id as keys
'''
def Chain_seq(file, combined_chain_id):
    # Combine all the ids together
    ids = combined_chain_id[0] + combined_chain_id[1] + combined_chain_id[2]
    # creat an empty dictionary, set a tracker to track whether an aa should be
    # added to the seq of a partitular chain
    seq = {}
    tracker = {}
    for i in ids:
        seq[i] = []
        tracker[i] = ''
    # load the sequences
    for line in file:
        if line[:6] == 'ATOM  ' and line[21] in ids:
            """Make sure only record the aa when the position number changes"""
            if tracker[line[21]] != line[22:27]:
                seq[line[21]].append(line[17:20])
                tracker[line[21]] = line[22:27]
    return seq

'''take a look'''
sequence = {}
for i in good_combined_ids:
    with open(i + '.pdb', 'r') as file:
        sequence[i] = Chain_seq(file, good_combined_ids[i]) 
len(sequence)
sequence.keys()
sequence['1adq'].keys()
len(sequence['1adq']['A'])
###############################################################

'''
Coordinates is to extract the coordinates used to calculate the interactions.
inputs: file, a pdb file
        id_dict, combined_chain_id, a list of the form ['BDF', 'ACE', 'GH']
        in the order of heavy chains, light chains, and antigen chains
return: cdn, a dictionary in the form of with keys ['h1H', 'h2H', 'h3H',
         'l1L', 'l1L', 'l1L', ..Antigen chain ids..]
        and the coordinates are in the form of [15.1, 2.2, 3.2, pos, aa]
        pos is an integer, indicates the position in the corresponding chain.
        aa, is the name of the amino acid.
'''
def Coordinates(file, combined_chain_id):
    # creat an empty dictionary to contain the results
    cdn = {}
    for i in combined_chain_id[0]:
        cdn['h1'+i], cdn['h2'+i], cdn['h3'+i] = [], [], []
    for i in combined_chain_id[1]:
        cdn['l1'+i], cdn['l2'+i], cdn['l3'+i] = [], [], []
    for i in combined_chain_id[2]:
        cdn[i] = []
        
    # creat a tracker dictionary, and a counter dictionary
    tracker = {}
    counter = {}
    ids = combined_chain_id[0] + combined_chain_id[1] + combined_chain_id[2]
    for i in ids:
        tracker[i] = ''
        counter[i] = -1
        
    # creat a dictionary to indicate the types of chains
    chain_type = {}
    for i in combined_chain_id[0]:
        chain_type[i] = 'H'
    for i in combined_chain_id[1]:
        chain_type[i] = 'L'
    for i in combined_chain_id[2]:
        chain_type[i] = 'A'
    
    # set the range of CDRh and CDRl, all the numbers take the same counting system
    # as python, with the firt one numbered 0.
    l_range = [[23, 35], [45, 55], [88, 96]]
    h_range = [[25, 35], [46, 64], [90, 109]]
    
    # extract the coordinates
    for line in file:
        if line[:6] == 'ATOM  ' and line[21] in ids:
            # update the parameters
            if tracker[line[21]]!= line[22:27]:
                counter[line[21]] += 1
                tracker[line[21]] = line[22:27]
            # extract all the parameters corresponding to line[21]
            c_type = chain_type[line[21]]
            count = counter[line[21]]
            # collect the coordinates according to c_type
            if c_type == 'H':
                #Tell the CDR type and load
                if count in range(h_range[0][0], h_range[0][1]+1):
                    cdn['h1'+line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])
                if count in range(h_range[1][0], h_range[1][1]+1):
                    cdn['h2'+line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])
                if count in range(h_range[2][0], h_range[2][1]+1):
                    cdn['h3'+line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])
    
            if c_type == 'L':
                #Tell the CDR type and load
                if count in range(l_range[0][0], l_range[0][1]+1):
                    cdn['l1'+line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])
                if count in range(l_range[1][0], l_range[1][1]+1):
                    cdn['l2'+line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])
                if count in range(l_range[2][0], l_range[2][1]+1):
                    cdn['l3'+line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])
            if c_type == 'A':
                cdn[line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]),
                                       count, line[17:20]])            
    
    return cdn
'''Extract all the coordinates and store them in dictionary coordinates with keys
   pdbid and the elements cdn
'''
coordinates = {}
for i in good_combined_ids:
    with open(i + '.pdb', 'r') as file:
        coordinates[i] = Coordinates(file, good_combined_ids[i])
'''Take a look at the results'''
len(coordinates)
coordinates.keys()
coordinates['1adq'].keys()
coordinates['1adq']['h1H']
##############################################################
# Extract the contact
'''
Get_contact is to calculate the contact number between amino acids. This function
        is time consuming.
inputs, cdn, a dictionary in the form of with keys ['h1H', 'h2H', 'h3H',
         'l1L', 'l1L', 'l1L', ..Antigen chain ids..] 
         and the coordinates are in the form of [15.1, 2.2, 3.2, pos, aa]
        pos is an integer, indicates the position in the corresponding chain
        aa, is the name of the amino acid.
        cutoff, a float, gives the cutoff distance
        matched_ids, a list in the form of [[H,L,A], [L, M, N]], where 
        [H, L, A] means those three are in a contacting group
return: contact, a list, in the form of [[h1HA, 32, 15, 8], ....]
        this means, the amino acid at position 32, which is located at CDRh1 of chain H, 
        contact with amino acid at position 15 of chain 'A'. The contact number 
        under the given cutoff is 8. The contact number is calculated by the following way:
            if atomA1 from aaA contacts with atomB1 from aaB, then the contact number 
            between aaA and aaB increased by 1. The contact between atomA1 and atomB1
            is only counted once.
'''
def Get_contact(cdn, matched_ids, cutoff = 4):
    # Creat an empty list to contain the temporary results
    contact_temp = []
    squared_cutoff = cutoff**2
    # sorting the keys into CDR and antigen groups
    # it is better to use the information of the matched ids
    # the grouped results should be stored in the form of[ [[h1H, h2H,h3H], [A]], ...]
    grouped =[]
    for matched in matched_ids:
        if matched[2] != '':
            if matched[0] != '':
                grouped.append([['h1'+matched[0], 'h2'+matched[0], 'h3'+matched[0]], [matched[2]]])
            if matched[1] != '':
                grouped.append([['l1'+matched[1], 'l2'+matched[1], 'l3'+matched[1]], [matched[2]]])
    #calculate the contact according to the grouped
    for match in grouped: 
    # calculate the distance and iterating through all possible combinations
        for i in match[0]:
            for j in match[1]:
                for atom1 in cdn[i]:
                    for atom2 in cdn[j]:
                        # We can accelerate this process by selecting the max abs first
                        diff = [atom1[0]-atom2[0],atom1[1]-atom2[1],atom1[2]-atom2[2]]                        
                        # is it faster to compare the square than the sequare root?
                        s = 0
                        for component in diff:
                            s += component**2
                            if s > squared_cutoff:# this step can accelerate the calculation by a lot.
                                break                        
                        if s <= squared_cutoff:
                            contact_temp.append([i+j, atom1[3], atom2[3]])         
    # Count method: Creat a dictionary to count\
    contact = []
    count_dict = {}
    for i in contact_temp:
        string = i[0] + '_' + str(i[1]) + '_' + str(i[2])
        if string in count_dict:
            count_dict[string] += 1
        else:
            count_dict[string] = 1
    # change the count_dict to contact
    contact = []
    for i in count_dict:
        element = i.split('_')
        element[1] = int(element[1])
        element[2] = int(element[2])
        element.append(count_dict[i])
        contact.append(element)
            
    return contact
##################################################################
# do the calculation directly, see what happens
# collect all the pdb ids
import time
start =time.clock()
contact = {}
n = 0
for i in good_matched_ids:
    n += 1
    print('Calculating   ' + i + '     ' + str(n))
    contact[i] = Get_contact(coordinates[i], good_matched_ids[i], cutoff = 4)
end = time.clock()
print('Running time: %s Seconds'%(end-start))
# remove the dud before saving
dud_AAC = []
for pdbid in contact:
    if contact[pdbid] == []:
        dud_AAC.append(pdbid)
dud_AAC
for dud in dud_AAC:
    del contact[dud]
    
'''Save the results'''
# save the results: sequence, contact
#with open('sequence', 'w') as f:
#    json.dump(sequence, f)
#with open('contact', 'w') as f:
#    json.dump(contact, f)
#with open('dud_AAC', 'w') as f:
#    json.dump(dud_AAC, f)

#with open('sequence', 'r') as f:
#    old_sequence = json.load(f)
#with open('contact', 'r') as f:
#    old_contact = json.load(f)
#with open('dud_AAC', 'r') as f:
#    old_dud_AAC = json.load(f)
    
#if old_sequence == sequence:
#    print('The sequences are the same!')
#else:
#    print('The sequences are different!')
#    
#if old_contact == contact:
#    print('The contacts are the same!')
#else:
#    print('The contacts are different!')
#    
#if old_dud_AAC == dud_AAC:
#    print('The dud_AACs are the same!')
#else:
#    print('The dud_AACs are different!')