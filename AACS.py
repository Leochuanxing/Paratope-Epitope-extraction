# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 19:59:59 2018

@author: leo
"""
from math import sqrt
import os
os.getcwd()
os.chdir("C:/Users/leo/Documents/Research/Database/PDB Learning")
os.getcwd()

def findChain(pdb, chain_id):
    chain = []
    for line in pdb:
        if line[:6] == 'SEQRES' and line[11] == chain_id:
            chain.extend(line[19:70].strip().split(' '))
        elif len(chain) > 0:
            break
    return chain

#def getCoordinates(pdb, chainId, Interval):
#    CDNT = []
#    m = 0
#    for line in pdb:
#        if line[:4] == 'ATOM' and line[21] == chainId:
#            n = int(line[22:26])
#            if m != n and m == 0:
#                m = n                
#            if (n-m) in Interval:
#                CDNT.append([float(line[30:38]),float(line[38:46]), float(line[46:54]), n-m])
#    return CDNT
#def getCoordinates(pdb, chain_id, residue_index):
#    CDNT = []
#    for line in pdb:
#        if line[:4] == 'ATOM' and line[21] == chain_id:
#            resSeq = int(line[22:26])
#            if resSeq - 1 in residue_index:
#                CDNT.append((float(line[30:38]), float(line[38:46]), float(line[46:54]), resSeq))
#        elif len(CDNT) > 0:
#            break
#    return CDNT
def getCoordinates(pdb, chain_id, residue_index):
    CDNT = []
    start_index = 0
    for line in pdb:
        if line[:4] == 'ATOM' and line[21] == chain_id:
            resSeq = int(line[22:26])
            if start_index == 0:
                start_index = resSeq
            if resSeq - start_index in residue_index:
                CDNT.append((float(line[30:38]), float(line[38:46]), float(line[46:54]), resSeq - start_index))
        elif len(CDNT) > 0:
            break
    return CDNT

# The returned format: [(HLchain index, Antigen index, contact number), ....]
def findContact(CDRcdnt, acdnt, cutoff):
    contactA = []
    for i in CDRcdnt:
        for j in acdnt:
            if sqrt((i[0] - j[0]) ** 2 + (i[1] - j[1]) ** 2 + (i[2] - j[2]) ** 2) <= cutoff:
                contactA.append((i[3], j[3]))
    contactB = list(set(contactA))
    contactB.sort(key=contactA.index)
    contactN = []
    if contactB != []:
        for i in contactB:
            j = i.__add__((contactA.count(i),))
            contactN.append(j)
    return   contactN

# Give the four coordinates contact for a given CDRid
#    In the form ('CDRid HLchain_id Achain_id', HLchain_index, Achain_index, contact_number )
#    CDRid takes the value of l1, l2, l3, h1, h2, h3.

def four_coordinates (CDR_id, HLchain_id, Achain_id, contact):
    four_coordinates = []
    for i in contact:
        if i != ():
            j = (CDR_id+HLchain_id+ Achain_id,).__add__(i)
            four_coordinates.append(j)
    return four_coordinates

# Assigns the ranges of CDRs, which are the largest among Kabat, Chothia, Contact and Abm, as
#    given by http://www.biochem.ucl.ac.uk/~martin/abs/GeneralInfo.html
CDRLindex = [list(range(23, 36)), list(range(45, 56)), list(range(88, 97))]
CDRHindex = [list(range(25, 36)), list(range(46, 65)), list(range(90, 110))]

# the input ID is in the form [pdbid, Hchainid, Lchainid, Achainid]
# the input pdb is the pdbid
# the return values are all contacts in four coordinates in four coordinates form
def contact(pdb, ID, cutoff=6):
    LAcontact1 = []
    LAcontact2 = []
    LAcontact3 = []
    HAcontact1 = []
    HAcontact2 = []
    HAcontact3 = []
    if ID[3] != '': 
        chainA = findChain(pdb, ID[3])
        
        cdnt_a = getCoordinates(pdb, ID[3], list(range(len(chainA))))

        if ID[2] != '':
#            chainL = findChain(pdb, ID[2])
            
            cdnt_l1 = getCoordinates(pdb, ID[2], CDRLindex[0])
            cdnt_l2 = getCoordinates(pdb, ID[2], CDRLindex[1])
            cdnt_l3 = getCoordinates(pdb, ID[2], CDRLindex[2])
            
            cnt_l1 = findContact(cdnt_l1, cdnt_a, cutoff)
            cnt_l2 = findContact(cdnt_l2, cdnt_a, cutoff)
            cnt_l3 = findContact(cdnt_l3, cdnt_a, cutoff)
            
            LAcontact1 = four_coordinates('l1', ID[2], ID[3], cnt_l1)
            LAcontact2 = four_coordinates('l2', ID[2], ID[3], cnt_l2)
            LAcontact3 = four_coordinates('l3', ID[2], ID[3], cnt_l3)  
               
        if ID[1] != '':    
#            chainH = findChain(pdb, ID[1])
            
            cdnt_h1 = getCoordinates(pdb, ID[1], CDRHindex[0])
            cdnt_h2 = getCoordinates(pdb, ID[1], CDRHindex[1])
            cdnt_h3 = getCoordinates(pdb, ID[1], CDRHindex[2])
  
            cnt_h1 = findContact(cdnt_h1, cdnt_a, cutoff)
            cnt_h2 = findContact(cdnt_h2, cdnt_a, cutoff)
            cnt_h3 = findContact(cdnt_h3, cdnt_a, cutoff)
            
            HAcontact1 = four_coordinates('h1', ID[1], ID[3], cnt_h1)
            HAcontact2 = four_coordinates('h2', ID[1], ID[3], cnt_h2)
            HAcontact3 = four_coordinates('h3', ID[1], ID[3], cnt_h3)             
  
    return LAcontact1, LAcontact2, LAcontact3, HAcontact1, HAcontact2, HAcontact3 

    
# open the summary file of the ids and find the corresponding ids in current directory
# the returned values are ids in this file in the form of list as 
# [['5kel', 'Q', 'U', 'I'], ['5kel', 'J', 'N', 'E'],...]

def get_ids_in_summary():    
    summary = open('summary.TSV', 'r')
    file = summary.readlines()        
    summary.close
    
    import IDhere 
    idhere = IDhere.main(file)
    
    return idhere


# the input is the id in this file, the out put is a dictionary with key pdbid
# and values a list of ids in the form as [['5kel', 'Q', 'U', 'I'], ['5kel', 'J', 'N', 'E'],...]
def clustered_id_by_pdbid (idhere):       
    dictid = {}
    for i in idhere:
        if dictid.__contains__(i[0]):
            dictid[i[0]].append(i)
        else:
            dictid[i[0]] = [i]
    return dictid

# the input is the ids here in the formats of dictionary and cutoff value
# the returned value are all the contacts in the format of dictionary with keys 
#    pdbids and values lists of four coordinates contact
def contact_dict(dictid, cutoff):      
    contact_dict = {}     
    pdbids_here_key = dictid.keys() 
    for k in pdbids_here_key:        
        
        f = open(k + '.pdb', 'r') 
        pdb = f.readlines()
        f.close
        
        res = []           
        for i in dictid[k]: 
            temp = contact(pdb, i, cutoff)
            for j in temp :
                res.extend(j)
        contact_dict[k] = res
    return contact_dict


idhere = get_ids_in_summary()
dictid = clustered_id_by_pdbid(idhere)
print(dictid) 
contactdict = contact_dict(dictid, 6)
print(contactdict)
print(dictid)


import csv
with open('contactdict1.csv', 'w') as csv_file:
    writer = csv.writer(csv_file)
    for key, value in contactdict.items():
       writer.writerow([key, value])


print(contactdict['5kel'])
len(contactdict['5kel'])

print(contactdict['1dee'])
len(contactdict['1dee'])

#
#
#idhere = get_ids_in_summary()
#print(idhere)
#dictid = clustered_id_by_pdbid(idhere)
#dee = dictid['1dee']
#dee_dict = {}
#dee_dict['1dee'] = dee
#dee_contact = contact_dict(dee_dict, 6)
#print(dee_contact)
#print(dee_dict)
#print(dee)
#
#f = open('1dee.pdb', 'r')
#pdb = f.readlines()
#ID = dee[2]
#print(ID)
#cutoff = 6
#chainA = findChain(pdb, ID[3])
#chainH = findChain(pdb, ID[1])
#chainL = findChain(pdb, ID[2])
#
#cdnt_a = getCoordinates(pdb, ID[3], list(range(len(chainA))))
#        
#cdnt_l1 = getCoordinates(pdb, ID[2], CDRLindex[0])
#cdnt_l2 = getCoordinates(pdb, ID[2], CDRLindex[1])
#cdnt_l3 = getCoordinates(pdb, ID[2], CDRLindex[2])
#
#cnt_l1 = findContact(cdnt_l1, cdnt_a, cutoff)
#cnt_l2 = findContact(cdnt_l2, cdnt_a, cutoff)
#cnt_l3 = findContact(cdnt_l3, cdnt_a, cutoff)
#
#LAcontact1 = four_coordinates('l1', ID[2], ID[3], cnt_l1)
#LAcontact2 = four_coordinates('l2', ID[2], ID[3], cnt_l2)
#LAcontact3 = four_coordinates('l3', ID[2], ID[3], cnt_l3)
#
#print(chainA)
#print(chainH)
#print(chainL)
#
#print(cdnt_a)
#print(cdnt_l1)
#print(cdnt_l2)
#print(cdnt_l3)
#
#print(cnt_l1)
#
#print(LAcontact1)
#print(LAcontact2)
#print(LAcontact3)
#import csv
with open('contactdict1.csv', 'w') as csv_file:
    keys = contactdict.keys()
    for i in contactdict.keys():
        j = [i]
        k = []
        for l in contactdict[i]:
            k.extend(list(l))
        j.extend(k)
        wr = csv.writer(csv_file, dialect='excel')
        wr.writerow(j)
#
with open('contactdict1.csv', 'r') as csv_file:        
    reader = csv.reader(csv_file)
    contactdict1 = {}
    for line in reader:
        if line != []:
            contactdict1[line[0]] = line[1]









