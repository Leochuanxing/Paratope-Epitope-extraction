import os
os.chdir('C:\\Users\\leo\\Documents\\Research\\Database\\Pipeline')

import json
with open('ac_contact', 'r') as f:
    ac_contact = json.load(f)
ac_contact.keys() 
ac_contact['1deeFh']
#########################################################
# define a finction, to find the consetutive aa in the ref frame with a given 
# lenth
'''
inputs: ac_contact, a dictionary returned from AlignmentConstraint module
        length, an integer, gien the consecutive length of the ref frame
        ref_chain, a string, takes the value of "Ab" or "Ag"
returns: max_ctn_aa, a dictonary in the form of {'1adqHh':[['h1HA', 'Ag', [21, 22], num_contact, all_contact],...],...}
         In this case, 'h1HA' means CDRh1, and the name of the antigen chain is 'A', antibody chain is 'H' 
         'Ag' mean antigen, [21, 22], means the 21, 22 position 
         aa in the 'Ag' chain is the consecutive aas with the largest contact number 'num_contact'
         with frame length 'length'. 'all_contact' is the total contact number of CDRh1 of chain 'H'.
'''   
def Max_ctn_aa(ac_contact, length, ref_chain = 'Ag'):
    
    # Creat an empty container to contain the final results
    max_ctn_aa = {}
        
    for Ab_chian in ac_contact:
        #Creat an empty container fot this key
        max_ctn_aa[Ab_chian] = []
        # group the values of ac_contact into different CDR groups
        # with keys 'h1HA'... and values a list of fcdn
        # creat a temporary empty container 
        temp_dict = {}
        for fcdn in ac_contact[Ab_chian]:
            if fcdn[0] not in temp_dict:
                temp_dict[fcdn[0]] = [fcdn]
            else:
                temp_dict[fcdn[0]].append(fcdn) 
        #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
        # DO THE WORK UNDER 'temp_dict'
        # translate the ref_chain information into index 
        if ref_chain == 'Ab':
            ind = 1
        if ref_chain == 'Ag':
            ind = 2        
        # Into the environment of the values of temp_dict
        for cdr in temp_dict:#cdr in the form of 'h1HA'
            # creat a temp container to contain the amino acid positions
            temp_aapos = []
            for fcdn in temp_dict[cdr]:
                if fcdn[ind] not in temp_aapos:
                    temp_aapos.append(fcdn[ind])
            #sort temp_aapos for the convenience to find the consecutive aa
            temp_aapos.sort()
            
            # find the consecutive aa psotions with the given length
            # Creat a temporary container 'temp_aa_ctn', to contain the results
            temp_aa_ctn = []
            if len(temp_aapos) >= length:
                for i in range(len(temp_aapos) - length + 1):
                    '''Can be modified to allow 0-free or 1-free'''
                    if temp_aapos[i + length - 1] -  temp_aapos[i] == length - 1:
                        temp_aa_ctn.append(temp_aapos[i : i + length ])
                        
            #Find the consecutive aa with the largest contact number
            # Creat a contact number recorder 'n', and a temporary empty list to 
            # contain the possible consecutive aas with the largest contact number.
            n = 0
            temp_max_aa_ctn = None
            # Find the one with the largest contact number
            if temp_aa_ctn != []:
                for aactn in temp_aa_ctn:
                    s = 0
                    for pos in aactn:
                        for fcdn in temp_dict[cdr]:
                            if fcdn[ind] == pos:
                                s += fcdn[3]
                    if s > n:
                        n = s
                        temp_max_aa_ctn = aactn  
            # calculate the total contact number
            all_contact = 0
            for fcdn in temp_dict[cdr]:
                all_contact += fcdn[3]                        
            # load the 'temp_max_aa_ctn' to 'max_ctn_aa'
            max_ctn_aa[Ab_chian].append([cdr, ref_chain, temp_max_aa_ctn, n, all_contact])
         #Clear up the temporary variables  
        del temp_aapos , temp_aa_ctn, temp_max_aa_ctn, n 
    del temp_dict
        #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#                    
    return max_ctn_aa
#Check the function
max_ctn_aa = Max_ctn_aa(ac_contact, length = 2, ref_chain = 'Ag')       
max_ctn_aa.keys()
max_ctn_aa['1bvkBh']
with open('ParaEpi_pos_framed', 'r') as f:
    ParaEpi_pos_framed = json.load(f)
with open('ParaEpi_framed', 'r') as f:
    ParaEpi_framed = json.load(f)
ParaEpi_pos_framed['1adq']  
ParaEpi_framed['1adq']        
# It looks good.   
#########################################################################
#Find the interacting paratope-epitope amino acids acording to the return of 
#'Max_ctn_aa'. This function should have a parameter to control the length of the 
# amino acids in the other chain other than the chain specified in the'Max_ctn_aa'.
# and the returned results should be the amino acids swith the largest contact number.
# the amino acids should be consecutive, otherwise indicate the jump.
'''
inputs: max_ctn_aa, a dictionary returned from 'Max_ctn_aa'
       complement_length, an integer, gives the length of the complemetary chain.
       ac_contact, a dictionary returned from the AlignmentConstraint module
       seq, a dictionary, give all the sequences.{'1adq'{'A':['ALA', 'GLU', ...]...}}
       CN_gate, an integer, set a lower limit for an aa to be inluded in the returns
return:parepi_pos, a dictionary in the form of   {'1adqHh': ['h1HA', [30], [16, 17], 12, 17],.. }
          'h1HA' means CDRh1, antibody chain is 'H', antigen chain is 'A';
          [30], gives the positions of the consecutive paratope aas no longer than 'complement_length'
          [16, 17] gives the positions of the epitope aas. 
          14, is the number of contact between [30] and [16,17]
          17, is the total number of contact in 'h1HA'
      parepi_aa, a corresponding dictionary of parepi_pos, the only difference is the pos is translated
          into amino acids.
'''
def ParaEpi(max_ctn_aa, ac_contact, seq, complement_length, CN_gate = 1):
    # Creat empty containers for the final results
    parepi_pos = {}; parepi_aa = {}
    # extract the ref_chain information from max_ctn_aa(Now I feel the advantage of class)
    for i in max_ctn_aa:
        for j in max_ctn_aa[i]:
            ref_chain = j[1]
            break
        break
    # here we have to make sure the 'ind' can indicate the complement position of the
    # ref_chain. 'c_ind' means the complementary index
    if ref_chain == 'Ab':
        ind = 1; c_ind = 2
    if ref_chain == 'Ag':
        ind = 2; c_ind = 1        
    # Into a Ab_chian
    for Ab_chian in max_ctn_aa:
        # Creat empty containers
        parepi_pos[Ab_chian], parepi_aa[Ab_chian] = [], []
        #Load the parepi_pos[Ab_chian] first
        # filter the 'ac_contact[Ab_chian]' through 'CN_gate' , group them into different CDRs
        #and store the results in 'CN_gated_fcdn' in the form of {'h1HA': [fcdn,..],...}.
        CN_gated_dict = {}
        for fcdn in ac_contact[Ab_chian]:
            if fcdn[3] >= CN_gate:
                if fcdn[0] not in CN_gated_dict:
                    CN_gated_dict[fcdn[0]] = [fcdn]
                else:
                    CN_gated_dict[fcdn[0]].append(fcdn) 
        # use the information in max_ctn_aa to index into fcdn in CN_gated_dict
        for max_ctn in max_ctn_aa[Ab_chian]:
            # we should guranttee that max_ctn[2] exists
            if max_ctn[2] != None:
                # Into one CDR of the CN_gated_dict
                ''' we just throw away the least contact ones if the c_aa is too long'''
                # in order to do this , sort the list of fcdn first
                CN_gated_dict[max_ctn[0]].sort(key = lambda x:x[3], reverse = True)
                # extract all the completmentary positions contact with max_ctn
                # and store them in 'temp_c_pos'
                temp_c_pos = []
                for fcdn in CN_gated_dict[max_ctn[0]]:
                    if len(temp_c_pos) < complement_length and fcdn[ind] in max_ctn[2]:
                        if fcdn[c_ind] not in temp_c_pos:
                            temp_c_pos.append(fcdn[c_ind])
                temp_c_pos.sort()
#                calculate the number of contacts between the c_aas and the 
#                 creat a contact recorder 'c'
                c = 0
                for fcdn in CN_gated_dict[max_ctn[0]]:
                    if fcdn[ind] in max_ctn[2] and fcdn[c_ind]  in temp_c_pos :
                        c += fcdn[3]
                # Everything is prepared, load the data
                if ref_chain == 'Ab':
                   parepi_pos[Ab_chian].append([max_ctn[0], max_ctn[2], temp_c_pos, c, max_ctn[4]])
                if ref_chain == 'Ag':
                    parepi_pos[Ab_chian].append([max_ctn[0], temp_c_pos, max_ctn[2], c, max_ctn[4]])
#                    
#                #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
                # load the parepi_aa according to [max_ctn[0], temp_c_pos, max_ctn[2], c, max_ctn[4]]
                # when loading, if there are jumps in 'temp_c_pos', be sure to indicace
                # first, find out the chain ids
                if ref_chain == 'Ab':
                    cid = max_ctn[0][2]; c_cid = max_ctn[0][3]
                if ref_chain == 'Ag':
                    cid = max_ctn[0][3]; c_cid = max_ctn[0][2]
                # map to the seq
                # creat temporary containers
                temp_c_aa, temp_max_ctn_aa  = [], []
                # load the the 'temp_c_aa', but we have to calculate whether there are jumps
                temp_c_pos.sort()# can not sort this way, has to sort according to 
                # creat an indicator to show if ther are jumps
                indicator = temp_c_pos[0]
                temp_c_aa.append(seq[Ab_chian[:4]][c_cid][indicator])
                if len(temp_c_pos) >= 2:
                    for pos in range(1, len(temp_c_pos)):
                        if temp_c_pos[pos] - indicator >= 2:
                            temp_c_aa.append(str(temp_c_pos[pos] - indicator-1))
                            temp_c_aa.append(seq[Ab_chian[:4]][c_cid][temp_c_pos[pos]])
                            indicator = temp_c_pos[pos]
                        else:
                            temp_c_aa.append(seq[Ab_chian[:4]][c_cid][temp_c_pos[pos]])
                            indicator = temp_c_pos[pos]
                # load the the 'max_ctn_aa', don't have to calculate jumps
                for pos in max_ctn[2]:
                    temp_max_ctn_aa .append(seq[Ab_chian[:4]][cid][pos])
                # load to parepi_aa
                if ref_chain == 'Ab':
                   parepi_aa[Ab_chian].append([max_ctn[0], temp_max_ctn_aa , temp_c_aa, c, max_ctn[4]])
                if ref_chain == 'Ag':
                    parepi_aa[Ab_chian].append([max_ctn[0], temp_c_aa, temp_max_ctn_aa , c, max_ctn[4]])
    # Clean the house before left               
    del CN_gated_dict, temp_c_pos, c, temp_c_aa, temp_max_ctn_aa   
#    return indicator
    return parepi_pos, parepi_aa
# Check if it works
with open('sequence_older', 'r') as f:
    sequence = json.load(f)
parepi_pos, parepi_aa = ParaEpi(max_ctn_aa, ac_contact, seq =sequence, complement_length = 10, CN_gate = 1)
parepi_aa.keys()           
parepi_aa['5vebAh'] 
ParaEpi_framed['5veb']['h1AX']
ParaEpi_framed['5veb']['h2AX']
ParaEpi_framed['5veb']['h3AX']
parepi_pos['5vebAh']['h1AX']

ParaEpi_pos_framed['5veb']
sequence['5veb']['A'][102]         
#for AB_c in max_ctn_aa:
#    for i in max_ctn_aa[AB_c]:
#        print(i)
#        break
#    break
print(26/71)         

                    
                    


    

        
        
  
        
        
        
        