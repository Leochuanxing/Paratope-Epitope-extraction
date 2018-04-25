# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 15:18:42 2018

@author: leo
"""

import os
import json
os.getcwd()
os.chdir("C:/Users/leo/Documents/Research/Database/Humnized antibody/PDB DATA/Homo Sapiens with paptide 5+/try100")
os.getcwd()

'''Find the Complex Ids and Chain Ids'''
# the first element is pdbid, the second is heavy chain, the third is lightchain and 
# the last one is antigen chain[[pdbid, heavychian, lightchain, achain]...]


def FindIds(file):
    Ids = []
    for l in file:
        b = []
#        Make sure the line is long enough
        if len(l) >= 16:
            a = l.split('\t')            
#  separete the achains connected by |
            a4 = a[4].split('|')            
            for i in a4:
                c = ''
                c = c + i.strip()
#                Combine the requied ids in line l and store them in b
                b = [a[0].strip(), a[1].strip(), a[2].strip(), c]         
#                Remove NA in b
                for i in range(1, 4):
                    if b[i] == 'NA':
                        b[i] = '' 
                Ids.append(b) 
    return Ids

 
# creat a dictionary with all the chainid combined, the input is the output of FindIds 
def Ids_dict (ids):
    Ids = {}
    for s in ids:              
        if Ids.__contains__(s[0]):
            for i in range(0, 3):
                Ids[s[0]][i] = Ids[s[0]][i] + s[i+1]
     # Add this term  if there is no such term         
        else:           
            Ids[s[0]] = [s[1], s[2], s[3]]
    return Ids

'''Find PDB files in current directory'''
def PDB_here ():
    names = os.listdir()
    PDBfiles = []
    for f in names:
        if len(f) == 8 and f[5:8] == 'pdb':
            PDBfiles.append(f[:4])
    return PDBfiles


'''Find the ids for pdb files in current directory'''
# Input: PDBfiles is the PDB files in current directory
# Input: ids is the ids from the summary which can be given by FindIds
def Id_here(PDBfiles, ids):
    Id_in_current_file = []
    for i in PDBfiles:
        for j in ids:
            if j[0] == i:
                Id_in_current_file.append(j)
    return Id_in_current_file
 


def chain_ids(file):
    Ids = FindIds(file)
    PDBfiles = PDB_here()
    Id_in_current_file = Id_here(PDBfiles, Ids)
    Id_dictionary = Ids_dict(Id_in_current_file)
    return Id_dictionary
        #Chain_ids is a dictionary as{'2vyr': ['KEGIHJLF', '', 'ACDAB']}


def findChain(pdb, chain_id):
    chain = []
    for line in pdb:
        if line[:6] == 'SEQRES' and line[11] == chain_id:
            chain.extend(line[19:70].strip().split(' '))
        elif len(chain) > 0:
            break
    return chain

# creat a dictionary of chains
# the input chain_ids should be given by function chain_ids
def Chain_dict(chain_ids):  
    Chain_dict = {}        
    for i in list(chain_ids.keys()):
        
        file = open(i+'.pdb', 'r')
        pdb = file.readlines()
        file.close
        
        k = ''
        for j in chain_ids[i]:
            k = k + j
    #   iterate over k to get the chain sequence
        Chain_dict_sub = {}
        for l in k:
            chain = findChain(pdb, l)
            Chain_dict_sub[l] = chain
        Chain_dict[i] = Chain_dict_sub
    return Chain_dict
#The output {pdbid:{chainid:[sequence], chainid: [sequence]......}, ......}              

# the input dictionary <- Chain_dict['pdbid'], keylist <- Chain_dict['pdbid'].keys(), 
#        relation_list <- []
def recur_compare(dictionary, keylist, relation_list):
    if keylist != []:
        j = keylist[0]
        keylist.remove(keylist[0])
        lbox = [j]        
        
        for i in keylist:            
            if dictionary[i] == dictionary[j]:
                lbox.append(i)
                keylist.remove(i)
            print(lbox)
        relation_list.append(lbox) 
        
        return recur_compare(dictionary, keylist, relation_list)
    
    else:
        return (dictionary, keylist, relation_list)
# Ids of equal chains are given in the relation_list in the form of [['A', 'D'],...]
#        which means chian A and chain D are the same


def Chain_dict_with_relation(Chain_dict):
    for i in Chain_dict.keys():
        keylist = list(Chain_dict[i].keys())
        a = recur_compare(Chain_dict[i], keylist, []) 
        Chain_dict[i]['Relation'] = a[2]
    return Chain_dict


def main():
    summary = open('summary.TSV', 'r')
    file = summary.readlines()        
    summary.close

    chain_id = chain_ids(file)
    Chain_dic = Chain_dict(chain_id)
    Chain_dict_relation = Chain_dict_with_relation(Chain_dic)
    
    return Chain_dict_relation
#    return Chain_dict


if __name__ == '__main__':

    Chain_dict_relation = main()

with open('Chain_dict_relation.json', 'w') as f:        
    json.dump(Chain_dict_relation, f)
        
#with open('Chain_dict_relation.json', 'r') as f:       
#    data = json.load(f)        
        
   
    
        
        
        
        


    
    
