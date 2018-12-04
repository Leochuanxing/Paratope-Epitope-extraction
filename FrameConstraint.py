
import os
os.chdir("/home/leo/Documents/Database/Pipeline/Homo")
os.listdir()

import json
with open('ready_for_FrameConstraint', 'r') as f:
    ac_contact = json.load(f)
with open('sequence', 'r') as f:
    sequence = json.load(f)
len(ac_contact.keys())
ac_contact['6b5sHh'] 
#########################################################

##################################################################################

####################################################################################
'''
The above functions are two long, it is better to write them into class object
CN_gate:
    an integer
cn_gated:
    stores all the contacting amino acids with the contact number no less than CN_gate
Ag_all_ctn:
    continuous amino acids of Atigen epitopes satisfies the Ag_free_type requirement

'''
class FrameConstraint(object):

    def __init__(self, ac_contact, sequence,
                  CN_gate = 1, cn_gated = None, Ag_frame_length = None, Ag_free_type = None,
                 Ag_all_ctn = None, Ab_frame_length = None, Ab_free_type = None,Ab_all_ctn = None, 
                 parepi_pos = None, parepi_aa = None):
        self.ac_contact = ac_contact
        self.sequence = sequence
        self.CN_gate = CN_gate  
    '''
    CN_gated: To get rid of the amino acid contact with contact number lower than a given number
    Output:self.cn_gated, a dictionary in the same form as the self.ac_contact, the only difference
          is that the amino acids with contact number smaller than the given number is deleted
    '''             
        
    def CN_gated(self):
        # this function select the fcdns meet the lower limit of the CN_gate
        cn_gated = {}
        for Ab_chain in self.ac_contact:
            cn_gated[Ab_chain] = []
            for fcdn in self.ac_contact[Ab_chain]:
                if fcdn[3] >= self.CN_gate:
                    cn_gated[Ab_chain].append(fcdn )
        self.cn_gated = cn_gated  
    '''
    All_ctn: gives the position index satisfying the given length and free type constraint
    return: 
    '''
             
    def All_ctn(self, length, free_type = 0, chain_id = 'Ag'):
        ac_contact = self.cn_gated
        # get the posion of the chain id
        if chain_id == 'Ag':
            ind = 2
        if chain_id == 'Ab':
            ind = 1
            
        all_ctn = {}
        for Ab_chian in ac_contact:
            all_ctn[Ab_chian] = {}# will be returned
            # group the values of ac_contact into different CDR groups
            # with keys 'h1HA'... and values a list of fcdn
            # creat a temporary empty container 
            temp_dict = {}
            for fcdn in ac_contact[Ab_chian]:
                if fcdn[0] not in temp_dict:
                    temp_dict[fcdn[0]] = [fcdn]
                else:
                    temp_dict[fcdn[0]].append(fcdn) 
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////#
            # DO THE WORK UNDER 'temp_dict'      
            # Go Into the environment of the values of temp_dict
            for cdr in temp_dict:#cdr in the form of 'h1HA'
                all_ctn[Ab_chian][cdr] = []# will be returned
                # creat a temp container temp_aapos to contain the amino acid positions
                pos = []
                for fcdn in temp_dict[cdr]:
                    if fcdn[ind] not in pos:
                        pos.append(fcdn[ind])
                #sort Ab_pos, and Ag_pos for the convenience to find the consecutive aa
                pos.sort() 
                # find the consecutive aa psotions with the given length
                # here we have to consider the free_type
                # Lets do it this way, first we partition the sequences into many contious
                # length according to the free_type, then we select the continuous amino acids
                # meet our length requirement
                # Creat empty contianers 'Ab_partitions' and 'Ag_partitions' to contain those partitioned continuous sequences
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
                partitions =[]
                #Do the partition
                while pos != []:
                    # Use temp to contain one partition
                    temp = []
                    # the beginning position is start
                    start = pos[0]
                    # load this value to temp
                    temp.append(start)
                    # if the len(temp_aapos) >= 2, pull out the longest sequence with start
                    # as the beginning
                    if len(pos) >= 2:
                        for i in range(1, len(pos)):
                            if pos[i] - pos[i - 1] <= free_type + 1:# free_type here
                                temp.append(pos[i])
                            # if the two neighbors has a difference largee than free_type + 1
                            # this longest sequence is end.
                            else:
                                break
                    # don't forget to remove the selected, so that this while loop can end
                    for d in temp:
                        pos.remove(d)
                    # we are done with the temp, load this temp to partitions
                    partitions.append(temp)
                # Ok, now lets find all the consecutive amino acids with the given length
                # from partitions, which the longest consecutive sequenences.
                for longest in partitions:
                    if len(longest) >= length:
                        for i in range(len(longest) - length+1):
                            all_ctn[Ab_chian][cdr].append(longest[i: i+length])
                            # the resluts of slicing a list is a list as well
        if chain_id == 'Ag':
            self.Ag_frame_length = length
            self.Ag_free_type = free_type
            self.Ag_all_ctn = all_ctn
        if chain_id == 'Ab':
            self.Ab_frame_length = length
            self.Ab_free_type = free_type
            self.Ab_all_ctn = all_ctn

    def Parepi_pos(self):
        # find the match with the largest contact
        Ab_all_ctn = self.Ab_all_ctn
        Ag_all_ctn = self.Ag_all_ctn
        ac_contact = self.ac_contact
        parepi_pos = {}
        for Ab_chain in ac_contact:
            parepi_pos[Ab_chain] = []                   
            for cdr in Ab_all_ctn[Ab_chain]:
                # the max math is taken over all cdrs
                max_c = 0
                max_match = [cdr, [], [], max_c]            
                for Ab_ctn in Ab_all_ctn[Ab_chain][cdr]:
                #Ab_ctn in the form of [359, 360]
                    for Ag_ctn in Ag_all_ctn[Ab_chain][cdr]:
                         #Ag_ctn in the form of [359, 360]
                        # find the max contact and corresponding mathed pos                    
                        c = 0
                        Ab_ctn_select = []
                        Ag_ctn_select = []
                        for fcdn in ac_contact[Ab_chain]:
                            if fcdn[1] in Ab_ctn and fcdn[2] in Ag_ctn:
                                c += fcdn[3]
                                # in case there is the case like[1,2][3,4], the contact
                                # between [1] and [3,4] is the largest among all 2-2 pattern
                                # while there is no contact between [2] and [3,4], we have to 
                                # select the fcdn again
                                Ab_ctn_select.append(fcdn[1])
                                Ag_ctn_select.append(fcdn[2])
                        if set(Ab_ctn_select) == set(Ab_ctn) and set(Ag_ctn_select) == set(Ag_ctn):
                            # Under this condition, the above problem can be avoided
                            if c > max_c:
                                max_c = c
                                max_match = [cdr, Ab_ctn, Ag_ctn, max_c]        
                # calculat all the contact in this cdr and add it to the returned results
                # adjust the orders of the Ab_ctn. Take the first element and the last element 
                #of Ab_ctn and calculate and their contacts with the first and the last elements
                #of the Ag_ctn, say f1 and f2, and l1, l2
                # if f1 + l2 > f2 + l1, keep the order of Ab_ctn, otherwise, reverse its order
                #
                max_match[1].sort()
                max_match[2].sort()
                if len(max_match[1]) >= 2 and len(max_match[2]) >= 2:
                    first_Ab = max_match[1][0]
                    last_Ab = max_match[1][-1]
                    first_Ag = max_match[2][0]
                    last_Ag = max_match[2][-1]
    #                print(first_Ag, last_Ag)
                    f1 = 0; f2 = 0; l1 = 0; l2 = 0 
                    for fcdn in ac_contact[Ab_chain]:
                        if fcdn[1] == first_Ab and fcdn[2] == first_Ag:
                            f1 = fcdn[3]
                        if fcdn[1] == first_Ab and fcdn[2] == last_Ag:
                            f2 = fcdn[3]
                        if fcdn[1] == last_Ab and fcdn[2] == first_Ag:
                            l1 = fcdn[3]
                        if fcdn[1] == last_Ab and fcdn[2] == last_Ag:
                            l2 = fcdn[3]
                    if f1 + l2 < f2 + l1:
                        max_match[1].sort(reverse = True)
                 #calculate the total number of contact in this cdr
                n = 0
                for fcdn in ac_contact[Ab_chain]:
                    if fcdn[0] == cdr:
                        n += fcdn[3]
                max_match.append(n)
                # load to the returned result                    
                parepi_pos[Ab_chain].append(max_match)
        self.parepi_pos = parepi_pos
        
    def Parepi_aa(self):
        parepi_pos = self.parepi_pos
        sequence = self.sequence
        parepi_aa = []
        # mapping the amino acids from the sequence, while the returned results should
        # be in the form of a list such as [[['TYR', 'VAL'], ['ASN'], 13, 41]...]
        # ['TYR', 'VAL'] is from the antibody, and ['ASN'] is from the antigen
        for Ab_chain in parepi_pos:
            for parepi in parepi_pos[Ab_chain]:
                # set a condition to make sure there is no empty Ab_aas or Ag_aas
                if parepi[1] != [] and parepi[2] != []:
                    # now continue the mapping process under this condition
                    Ab_aas = []
                    for Ab_pos in parepi[1]:
                        Ab_aas.append(sequence[Ab_chain[:4]][parepi[0][2]][Ab_pos])
                    Ag_aas = []
                    for Ag_pos in parepi[2]:
                        Ag_aas.append(sequence[Ab_chain[:4]][parepi[0][3]][Ag_pos])
                    parepi_aa.append([Ab_aas, Ag_aas, parepi[3], parepi[4]])
        self.parepi_aa = parepi_aa                
                        
        
testing = FrameConstraint(ac_contact, sequence, CN_gate = 1)     
testing.CN_gated() 
len(testing.cn_gated)           
testing.cn_gated.keys()       
testing.All_ctn(length = 2, free_type = 1, chain_id = 'Ag')        
testing.Ag_all_ctn.keys()
testing.Ag_all_ctn['5yaxBl']
len(testing.Ag_all_ctn)
#testing.Ag_all_ctn['6eluKh'] 
#testing.ac_contact['6eluKh']       
testing.Ag_free_type    
testing.Ag_frame_length
testing.All_ctn(length = 2, free_type = 0, chain_id = 'Ab')    
testing.Ab_all_ctn.keys()
len(testing.Ab_all_ctn)  
#testing.Ab_all_ctn['6eluKh']   
testing.Ab_frame_length    
testing.Ab_free_type    
testing.Parepi_pos()  
testing.parepi_pos.keys()
#testing.parepi_pos['6eluKh']   
testing.Parepi_aa()    
len(testing.parepi_aa)
testing.parepi_aa.sort()
testing.parepi_aa[:10]
#with open('parepi_aa_2_2', 'r') as f:
#    parepi_aa_2_2 = json.load(f)    
#parepi_aa_2_2['6eluKh']
#sequence['6elu']['K'][102]
#sequence['6elu']['J'][164]      
#    
os.chdir("/home/leo/Documents/Database/Pipeline/Ready_2_2_0_1")
os.listdir()    
with open('ready_2_2_0_1', 'w') as f:
    json.dump(testing.parepi_aa, f)
    
    
    
    
    
    
    

                    


    

        
        
  
        
        
        
        