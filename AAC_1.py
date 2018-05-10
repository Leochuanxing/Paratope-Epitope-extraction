#import os
#os.getcwd()
#os.chdir('./Documents/Research/Database/PDB Learning')
#os.listdir()
#   
class Parsed_antibody_antigen_complex:
    def __init__(self, PDBseq, Coordinates):
        self.PDBseq = PDBseq
        self.Coordinates = Coordinates

def Find_Chain_Coordinates(pdb, combined_chain_id): #combined_chain_id = [heavy,light, antigen]
    PDBseq = {}
    DBREF = {}
    Coordinates = {}
    counter = 0
    normal_tracker = 0
    insersion_tracker = ' '
    CDRLindex = [list(range(23, 36)), list(range(45, 56)), list(range(88, 97))]
    CDRHindex = [list(range(25, 36)), list(range(46, 65)), list(range(90, 110))]
    for line in pdb:
        
        if line[:5] == 'DBREF' and line[12] not in DBREF:
            DBREF[line[12]] = [int(line[14:18])] # DBREF{chainid:[begin_int]}
            
#         Find the sequences for chains            
        elif line[:6] == 'SEQRES':
            if line[11] not in PDBseq:
                PDBseq[line[11]] = line[19:70].strip().split(' ')
            else: 
                PDBseq[line[11]].extend(line[19:70].strip().split(' '))
           
#              Find the coordinates of CDRHs 
        elif (line[:4] == 'ATOM' ) and line[21] in combined_chain_id[0]:
            (DBREF, Coordinates) = get_cdr_coordinates(line, DBREF, Coordinates, 'h1', CDRHindex[0])
            (DBREF, Coordinates) = get_cdr_coordinates(line, DBREF, Coordinates, 'h2', CDRHindex[1])
            (DBREF, Coordinates) = get_cdr_coordinates(line, DBREF, Coordinates, 'h3', CDRHindex[2])
#              Find the coordinates of CDRLs 
        elif (line[:4] == 'ATOM' ) and line[21] in combined_chain_id[1]:
            (DBREF, Coordinates) = get_cdr_coordinates(line, DBREF, Coordinates, 'l1', CDRLindex[0])
            (DBREF, Coordinates) = get_cdr_coordinates(line, DBREF, Coordinates, 'l2', CDRLindex[1])
            (DBREF, Coordinates) = get_cdr_coordinates(line, DBREF, Coordinates, 'l3', CDRLindex[2])
#              Find the coordinates of Antigen
        elif (line[:4] == 'ATOM' ) and line[21] in combined_chain_id[2]:
#            (DBREF, Coordinates) = get_cdr_coordinates(line, DBREF, Coordinates, '', list(range(10000)))                
            (counter, normal_tracker, PDBseq, Coordinates, insersion_tracker) = Get_antigen_chain_and_coordinates(line, line[21], counter, normal_tracker, 
            PDBseq, Coordinates, insersion_tracker) 
    
    return  PDBseq, Coordinates


def get_cdr_coordinates(line, DBREF, Coordinates, cdr, index):
    insersion_tracker = 0
    if len(DBREF[line[21]]) == 1:
        DBREF[line[21]]=[int(line[22:26])-DBREF[line[21]][0], int(line[22:26]), line[26]] # DBREF{chainid:[ tracker_int, counter_str, insersion]}           
    else:        
        if DBREF[line[21]][2] != line[26] and line[26] != ' ':
            DBREF[line[21]][2] = line[26]
            insersion_tracker = 1 
        elif DBREF[line[21]][2] != line[26] and line[26] == ' ':
            DBREF[line[21]][2] = line[26]
            
        DBREF[line[21]][0] += insersion_tracker + (int(line[22:26])-DBREF[line[21]][1])
        DBREF[line[21]][1] = int(line[22:26])
        
    if (DBREF[line[21]][0] in index) and ((cdr+line[21]) not in Coordinates):
        Coordinates[cdr+line[21]] = [[float(line[30:38]), float(line[38:46]), float(line[46:54]), DBREF[line[21]][0]]]
    elif(DBREF[line[21]][0] in index) and ((cdr+line[21]) in Coordinates):    
        Coordinates[cdr+line[21]].append([float(line[30:38]), float(line[38:46]), float(line[46:54]), DBREF[line[21]][0]])
    return  DBREF, Coordinates    

def Get_antigen_chain_and_coordinates(line, chain_id, counter = 0, normal_tracker = 0, 
    antigen_chain = {}, antigen_coordinates = {},insersion_tracker = ' '):
    #line begins with"ATOM", chain_id is for antigen
    break_indexer = 0 # track the distance between two adjacent  ATOM data in terms of amino acids
#     track the normal sequence, to see whether it jumps    
    if  normal_tracker != 0:
        break_indexer += int(line[22:26]) - normal_tracker
        normal_tracker = int(line[22:26])
    elif normal_tracker == 0:
        normal_tracker = int(line[22:26])
        antigen_chain[chain_id] = [line[17:20]]
# track whether there is any insersion        
    if line[26] != insersion_tracker and line[26] != ' ':
        break_indexer += 1
        insersion_tracker = line[26]
    elif line[26] != insersion_tracker and line[26] == ' ':
        insersion_tracker = line[26]
#     extract coordinates and sequences, if there is a break, insert 'BRK' in the sequence
    if break_indexer == 0:
        if chain_id in antigen_coordinates:
            antigen_coordinates[chain_id].append([float(line[30:38]), float(line[38:46]), float(line[46:54]), counter])
        else:
            antigen_coordinates[chain_id] = [[float(line[30:38]), float(line[38:46]), float(line[46:54]), counter]]
    elif break_indexer == 1:
        counter += 1
        antigen_chain[chain_id].append(line[17:20])
        if chain_id in antigen_coordinates:
            antigen_coordinates[chain_id].append([float(line[30:38]), float(line[38:46]), float(line[46:54]), counter])
        else:
            antigen_coordinates[chain_id] = [[float(line[30:38]), float(line[38:46]), float(line[46:54]), counter]]
    elif break_indexer >= 2:
        counter += 2
        antigen_chain[chain_id].extend(['BRK', line[17:20]])
        if chain_id in antigen_coordinates:
            antigen_coordinates[chain_id].append([float(line[30:38]), float(line[38:46]), float(line[46:54]), counter])
        else:
            antigen_coordinates[chain_id] = [[float(line[30:38]), float(line[38:46]), float(line[46:54]), counter]]
    
    return counter, normal_tracker, antigen_chain, antigen_coordinates, insersion_tracker
        
#with open('1g9m.pdb', 'r') as f:
#    data_1g9m = f.readlines()
#counter = 0
#normal_tracker = 0
#antigen_chain = {}
#antigen_coordinates = {}
#insersion_tracker = ' '
#for line in data_1g9m:
#    if line[:4] == 'ATOM':
#        (counter, normal_tracker, antigen_chain, antigen_coordinates, insersion_tracker) = Get_antigen_chain_and_coordinates(line, 'G', counter, normal_tracker, 
#        antigen_chain, antigen_coordinates,insersion_tracker)
#type(antigen_chain)
#antigen_chain.keys()
#j = 0       
#for i in antigen_chain['G']:
#    j += 1
#    if i == 'BRK':
#        print(j)
#antigen_chain[154]   
#antigen_chain[155]  
#antigen_coordinates.keys()
#for i in antigen_coordinates['G']: 
#    if i[3] == 48:
#        print(i)
    

def Distance(coordinate1, coordinate2):
    distance_square = 0
    for i in range(0,3):
        distance_square += (coordinate1[i]-coordinate2[i])**2
    distance = distance_square**0.5
    return distance

def findContact(coordinates_for_one_pdb, id_dict_for_one_pdb, cutoff):
    ct = cutoff
    contact_count = []
    for i in id_dict_for_one_pdb:#id_dict_for_one_pdb
#find contact between between CDRHs and the Antigen
        if i[2] != '' and i[0] != '':
            CDRHs = ['h1'+i[0], 'h2'+i[0], 'h3'+i[0]]
            contact_count.extend(findContact_sub_function(CDRHs, i[2], coordinates_for_one_pdb, ct))
#find contact between between CDRLs and the Antigen                   
        if i[2] != '' and i[1] != '':           
            CDRLs = ['l1'+i[1], 'l2'+i[1], 'l3'+i[1]]
            contact_count.extend(findContact_sub_function(CDRLs, i[2], coordinates_for_one_pdb, ct))           
    return contact_count          

def findContact_sub_function(CDRs, achain, coordinates_for_one_pdb, cutoff):
    contact_sub_all = []
    contact_sub_count =[]
    temp_dict = {}
    for j in CDRs:
        for k in coordinates_for_one_pdb[j]:
            for l in coordinates_for_one_pdb[achain]:
                if Distance(k[:3],l[:3]) <= cutoff:
                    contact_sub_all.append(j+achain+'_'+str(k[3])+'_'+str(l[3]))
    for m in contact_sub_all:
        if m in temp_dict:
            temp_dict[m] += 1
        else:
            temp_dict[m] = 1
    for n in temp_dict:
        temp_list = n.split('_')
        contact_sub_count.append([temp_list[0], int(temp_list[1]), int(temp_list[2]), temp_dict[n]])
    return contact_sub_count    
    
with open('summary.TSV', 'r') as summary:
    file = summary.readlines()
    
import IDhere
id_dict = IDhere.Id_dict(file)
combined_chain_id_dict = IDhere.Combined_chain_id_dict(id_dict)
here_iddict_combineddict = IDhere.Here_iddict_combineddict(id_dict, combined_chain_id_dict)

res = {}
for i in here_iddict_combineddict[1]:
    with open(i+'.pdb', 'r') as f:
        res[i] = Find_Chain_Coordinates(f, here_iddict_combineddict[1][i])

res_contact = {}
for i in res:
    res_contact[i] = findContact(res[i][1], here_iddict_combineddict[0][i], 5)


       

import json
with open ('Contact_dict_here_5A', 'w') as f:
    json.dump(res_contact, f)
with open ('Chain_and_coordinates_dict_here_5A', 'w') as f:
    json.dump(res, f)

with open('Contact_dict_here_5A', 'r') as data:
    data_check = json.load(data)

#def Contact_checking(coordinates1, coordinates2):
#    for i in coordinates1:
#        for j in coordinates2:
#            d = ((i[0]-j[0])**2+(i[1]-j[1])**2+(i[2]-j[2])**2)**0.5
#            if d <= 5:
#                print([i,j,d,'Contact'])
##            if d >= 5:
##                print([i,j,d,'No Contact'])
##            else:
##                print([i,j,'None Contact'])
#    return()
#
#def Temp_coordinate(coordinates_for_one_pdb, index, order_of_amino):
#    temp_coordinate = []
#    for i in coordinates_for_one_pdb[index]:
#        if i[3] == order_of_amino:
#            temp_coordinate.append(i)
#    return temp_coordinate
#
#res_contact['1g9m']
#res['1g9m'][0]['G'][61]
#res['1g9m'][0]['L'][93]
#len(res['1g9m'][0]['G'])
#for i in res['1g9m'][1]['G']:
#    if i[3] == 340:
#        print(i)
#print(res['1g9m'][1]['A'][:10])
#len(res['1g9m'][0]['G'])
#
#data_check['1g9m']
#res_contact   
#len(res['1g9m'][1]['A'])
#here_iddict_combineddict[0]['1g9m']
#
#keys = list(res.keys())
#keys[3]  
#temp_coordinate1 =  Temp_coordinate(res[keys[3]][1], 'l3L', 93) 
#temp_coordinate2 =  Temp_coordinate(res[keys[3]][1], 'G', 59) 
#print(temp_coordinate1)
#print(temp_coordinate2)
##
#Contact_checking(temp_coordinate1, temp_coordinate2)
#
#with open('1g9m', 'r') as file:
#    kel = file.readlines()
#for line in kel:
#    if line[:4] == 'ATOM' and float(line[30:38]) == 155.484:
#        print(line)
##    





