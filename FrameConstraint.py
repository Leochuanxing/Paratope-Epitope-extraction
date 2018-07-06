import os
os.chdir('C:\\Users\\leo\\Documents\\Research\\Database\\Pipeline')

import json
with open('ac_contact', 'r') as f:
    ac_contact = json.load(f)
ac_contact.keys() 
ac_contact['1bvkBh']
#########################################################
# define a finction, to find the consetutive aa in the ref frame with a given 
# lenth
'''
inputs: ac_contact, a dictionary returned from AlignmentConstraint module
        length, an integer, gien the consecutive length of the ref frame
        ref_chain, a string, takes the value of "Ab" or "Ag"
returns: max_ctn_aa, a dictonary in the form of {'1adqHh':[['h1', 'Ag', [21, 22], num_contact, all_contact],...],...}
         In this case, 'h1' means CDRh1, 'Ag' mean antigen, [21, 22], means the 21, 22 position 
         aa in the 'Ag' chain is the consecutive aas with the largest contact number 'num_contact'
         with frame length 'length'. 'all_contact' is the total contact number of CDRh1 of chain 'H', 
         pdb id '1adq'
         2, in contact region CDRh1
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
            max_ctn_aa[Ab_chian].append([cdr[:2], ref_chain, temp_max_aa_ctn, n, all_contact])
            
        del temp_aapos , temp_aa_ctn, temp_max_aa_ctn, n 
    del temp_dict
        #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#                    
    return max_ctn_aa
#Check the function
max_ctn_aa = Max_ctn_aa(ac_contact, length = 2, ref_chain = 'Ag')       
max_ctn_aa.keys()
max_ctn_aa['1adqHh']
with open('ParaEpi_pos_framed', 'r') as f:
    ParaEpi_pos_framed = json.load(f)
ParaEpi_pos_framed['1adq']           
# It looks good.       

        
        
  
        
        
        
        