import os
os.getcwd()
os.chdir("C:\\Users\\leo\\Documents\\Research\\Database\\Pipeline")
os.getcwd()
os.listdir()
###############################################
import json
with open('here_conbined_chain_id', 'r') as f:
    here_conbined_chain_id = json.load(f)
with open('here_id_dict', 'r') as f:
    here_id_dict = json.load(f)
    
here_conbined_chain_id['1dee']
################################################
#Eliminate the wierd pdb, and set aside those wierd pdb ids
# by weird, it means a chain id shows up more than once
'''
inputs: conbined_chain_id, a dictionary with keys pdb id, and values in the form
        of ['BDF', 'ACE', 'GH']
return: witches, a list of weird pdb ids
        GoodPeople, a dictionary with keys of good pdb ids, and values in the form
        of ['BDF', 'ACE', 'GH']
'''

def witch_hunt(conbined_chain_id):
    # creat a container to contain the witches and GoodPeople
    witches = []
    GoodPeople = conbined_chain_id
    for i in conbined_chain_id:
        antibody_ids = conbined_chain_id[i][0] + conbined_chain_id[i][1] 
        antigen_ids = conbined_chain_id[i][2]
        witch = False
        for i in antibody_ids:
            if i in antigen_ids:
                witch = True
                break
        if witch:
            witches.append(i)
            del GoodPeople[i]
    return witches, GoodPeople
# take a look at the results           
witches, GoodPeople = witch_hunt(here_conbined_chain_id)
witches # Here there is no witches, but in the large scale calculation, we have to use Goodpeople
GoodPeople.keys()
GoodPeople['1dee'] 
here_conbined_chain_id.keys()  
# creat good ids
good_combined_ids = {}
good_matched_ids = {}
for i in GoodPeople:
    good_combined_ids[i] = here_conbined_chain_id[i]
    good_matched_ids[i] = here_id_dict[i]
good_matched_ids
#################################################
# Extract the chain sequences 
'''
Imputs: file, a pdb file
       : combined_chain_id, a list of the form ['BDF', 'ACE', 'GH']
         in the order of heavy chains, light chains, and antigen chains
Returns: seq, a dictionary of sequences, with keys as the chain id
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
            if tracker[line[21]] != line[22:27]:
                seq[line[21]].append(line[17:20])
                tracker[line[21]] = line[22:27]
    return seq

# take a look
sequence = {}
for i in good_combined_ids:
    with open(i + '.pdb', 'r') as file:
        sequence[i] = Chain_seq(file, good_combined_ids[i])       
###############################################################
# extract the coordinates used to calculate the interactions
'''
inputs: file, a pdb file
        id_dict, combined_chain_id, a list of the form ['BDF', 'ACE', 'GH']
        in the order of heavy chains, light chains, and antigen chains
return: cdn, a dictionary in the form of with keys ['h1H', 'h2H', 'h3H',
         'l1L', 'l1L', 'l1L', ..Antigen chain ids..]
        and the coordinates are in the form of [15.1, 2.2, 3.2, pos, aa]
        pos is an integer, indicates the position in the corresponding chain
        aa, is the name of the amino acid.
'''
def Coordinates(file, combined_chain_id):
    # creat and empty dictionary to contain the results
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
    
    # extracting the coordinates
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
# take a look at the results
coordinates = {}
for i in good_combined_ids:
    with open(i + '.pdb', 'r') as file:
        coordinates[i] = Coordinates(file, good_combined_ids[i])
coordinates['1adq'].keys()
##############################################################
# Extract the contact
'''
inputs, cdn, a dictionary in the form of with keys ['h1H', 'h2H', 'h3H',
         'l1L', 'l1L', 'l1L', ..Antigen chain ids..] 
         and the coordinates are in the form of [15.1, 2.2, 3.2, pos, aa]
        pos is an integer, indicates the position in the corresponding chain
        aa, is the name of the amino acid.
        cutoff, a float, gives the cutoff distance
        matched_ids, a list in the form of [[H,L,A], [L, M, N]], where 
        [H, L, A] means those three are in a contacting group
return: contact, a list, in the form of [[h1HA, 32, 15, 8], ....]
        this means, the amino acid at position 32, which is located at CDRh1, 
        contact with amino acid at position 15 of chain 'A'. The contact number 
        under the given cutoff is 8
'''
def Get_contact(cdn, matched_ids, cutoff = 4):
    # Creat an empty list to contain the temporary results
    contact_temp = []
    squared_cutoff = cutoff**2
    # sorting the keys into CDR and antigen groups
    # it is better to use the information of the matched ids
    # the grouped results should be stored in the form of[ [[h1H, h2H,h3h], [A]], ...]
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
#                        if max(abs(diff)) < cutoff:
                            # is it faster to compare the square than the sequare root?
                        s = 0
                        for component in diff:
                            s += component**2
                            if s > squared_cutoff:# this step can accelerate the calculation by a lot.
                                break                        
                        if s <= squared_cutoff:
                            contact_temp.append([i+j, atom1[3], atom2[3]])
#                            if np.dot(diff,diff) <= squared_cutoff:
#                                contact_temp.append([i+j, atom1[3], atom2[3]])  
    #Count the contact number
    # Count method 1, use while to count
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#    contact = []
#    while contact_temp != []:
#        contact = contact_temp[0]
#        count = 0
#        while contact in contact_temp:
#            count += 1
#            contact_temp.remove(contact)
#        contact.append([contact[0], contact[1], contact[2], count])
            
    # Count method 2 Creat a dictionary to count\
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

# Take a look at the results, caltulate the running time
import time
start =time.clock()
contact = {}
for i in coordinates:
    contact[i] = Get_contact(coordinates[i], good_matched_ids[i], cutoff = 4)
end = time.clock()
print('Running time: %s Seconds'%(end-start))

# compare with the results of the previous method
#with open('contact_older', 'r') as f:
#    contact_older = json.load(f) 
#contact_older['1adq']
## systematic comparison
#different_id = []
#for i in contact:
#    for j in contact_older[i]:
#        if j not in contact[i] and i not in different_id:
#            different_id.append(i)
#different_id 
#contact['1g9m'][:6]  
#contact_older['1g9m'][:6]
## Check if the number of contacting amino acids are the same
#for i in contact:
#    if len(contact[i]) == len(contact_older[i]):
#        print ('The length of '+ i + ' are the same')
#    else:
#        print ('The length of '+ i + ' are different')
## Check if the total number of contacts are the same
#for keys in contact:
#    new_s = 0
#    for new in contact[keys]:
#        new_s += new[3]
#    old_s = 0
#    for old in contact_older[keys]:
#        old_s += old[3]
#    if new_s == old_s:
#        print(keys + '         same')
#    else:
#        print(keys + '         different')
#        
#
#contact['5kel']
#len(contact['5kel'])
#len(contact_older['5kel'])
#GoodPeople['5kel']
#with open('sequence_older', 'r') as f:
#    sequence_older = json.load(f)            
#sequence['1g9m'].keys()
#sequence_older['1g9m'].keys()
#sequence['1g9m']['G'][37]
#sequence_older['1g9m']['G'][37]
