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

def getCoordinates(pdb, chain_id, residue_index):
    CDNT = []
    for line in pdb:
        if line[:4] == 'ATOM' and line[21] == chain_id:
            resSeq = int(line[22:26])
            if resSeq - 1 in residue_index:
                CDNT.append((float(line[30:38]), float(line[38:46]), float(line[46:54]), resSeq))
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

# Give the six coordinates contact for a given CDRid
#    In the form (CDRid, HLchain_id, Achain_id, HLchain_index, Achain_index, contact_number )
#    CDRid takes the value of l1, l2, l3, h1, h2, h3.

def six_coordinates (CDR_id, HLchain_id, Achain_id, contact):
    six_coordinates = []
    for i in contact:
        if i != ():
            j = (CDR_id, HLchain_id, Achain_id).__add__(i)
            six_coordinates.append(j)
    return six_coordinates


# Structure of the output
#class ParEpi:
#    def __init__(self, PdbId, chainL, chainH, chainA, L1CNT,
#                 L2CNT, L3CNT, H1CNT, H2CNT, H3CNT):
#        self.PdbId = PdbId
#        self.chainL = chainL
#        self.chainH = chainH
#        self.chainA = chainA
#        self.L1CNT = L1CNT
#        self.L2CNT = L2CNT
#        self.L3CNT = L3CNT
#        self.H1CNT = H1CNT
#        self.H2CNT = H2CNT
#        self.H3CNT = H3CNT
#
#    def __str__(self):
#        return (('PDBid: {} \n  chainL: {}\n  chainH: {}\n chainA: {} \n L1CNT: {} \n' +
#                 'L2CNT:{} \n L3CNT:{} \n H1CNT: {} \n H2CNT: {} \n' +
#                 'H3CNT: {} \n')
#                .format(self.PdbId, self.chainL, self.chainH, self.chainA, self.L1CNT, self.L2CNT,
#                        self.L3CNT,self.H1CNT, self.H2CNT, self.H3CNT))

# See http://www.biochem.ucl.ac.uk/~martin/abs/GeneralInfo.html
# We use the minimum of the starting position and the maximum of the ending position
# of Kabat, Chothia, AbM and Contact. For H3, we creat a range with length 20
CDRLindex = [list(range(23, 36)), list(range(45, 56)), list(range(88, 97))]
CDRHindex = [list(range(25, 36)), list(range(46, 65)), list(range(90, 110))]

def contact(ID, cutoff=6):
    with open(ID[0] + '.pdb', 'r') as f:
        pdb = f.readlines()
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
            
            LAcontact1 = six_coordinates('l1', ID[2], ID[3], cnt_l1)
            LAcontact2 = six_coordinates('l2', ID[2], ID[3], cnt_l2)
            LAcontact3 = six_coordinates('l3', ID[2], ID[3], cnt_l3)  
               
        if ID[1] != '':    
#            chainH = findChain(pdb, ID[1])
            
            cdnt_h1 = getCoordinates(pdb, ID[1], CDRHindex[0])
            cdnt_h2 = getCoordinates(pdb, ID[1], CDRHindex[1])
            cdnt_h3 = getCoordinates(pdb, ID[1], CDRHindex[2])
  
            cnt_h1 = findContact(cdnt_h1, cdnt_a, cutoff)
            cnt_h2 = findContact(cdnt_h2, cdnt_a, cutoff)
            cnt_h3 = findContact(cdnt_h3, cdnt_a, cutoff)
            
            HAcontact1 = six_coordinates('h1', ID[1], ID[3], cnt_h1)
            HAcontact2 = six_coordinates('h2', ID[1], ID[3], cnt_h2)
            HAcontact3 = six_coordinates('h3', ID[1], ID[3], cnt_h3)             
  
        return LAcontact1, LAcontact2, LAcontact3, HAcontact1, HAcontact2, HAcontact3 
  
    
summary = open('summary.TSV', 'r')
file = summary.readlines()        
summary.close
        
import IDhere 
l = IDhere.main(file)
print(l)

kel = l[(len(l)-6):len(l)]
print(kel)

dict_5kel = {}
for s in kel:
    if dict_5kel.__contains__(s[0]):
            for i in range(0, 3):
                dict_5kel[s[0]][i] = dict_5kel[s[0]][i] + s[i+1] 
    else:
        dict_5kel[s[0]] = [s[1], s[2], s[3]]
print (dict_5kel)


'''Run a test to see if the classification is correct'''
ckel = dict_5kel['5kel']
print(ckel)

LAID = []
for i in ckel[1]:
    for j in ckel[2]:
        LAID.append(['5kel','', i, j])
print(LAID)
len(LAID)  

f = open('5kel.pdb', 'r') 
pdb = f.readlines()
f.close
 
ID = ['5kel', 'M', 'O', 'F'] 
for i in kel:
    print(contact(i, 6))
    




def LACNT (ID, cutoff = 6):
    cdnt_l1 = getCoordinates(pdb, ID[2], CDRLindex[0])
    cdnt_l2 = getCoordinates(pdb, ID[2], CDRLindex[1])
    cdnt_l3 = getCoordinates(pdb, ID[2], CDRLindex[2]) 
    
    chainA = findChain(pdb, ID[3])
    cdnt_a = getCoordinates(pdb, ID[3], list(range(len(chainA))))

    cnt_l1 = findContact(cdnt_l1, cdnt_a, cutoff)
    cnt_l2 = findContact(cdnt_l2, cdnt_a, cutoff)
    cnt_l3 = findContact(cdnt_l3, cdnt_a, cutoff)
   
    return [ID[2], ID[3]],[cnt_l1, cnt_l2, cnt_l3]     

print(LACNT(ID,  6))

for i in LAID:
    print(LACNT(i, 6))
 
print(kel)
problem = []
LACNT1 = []
for ID in kel:
    if ID[3] != '':
        try:
            Contact = main(ID, 6)                
            LACNT1.append([[ID[2],ID[3]], [Contact.L1CNT, Contact.L2CNT, Contact.L3CNT]])
        except:
            problem.append('Something is wrong with '+ ID[0]+'.pdb')
print(problem)
print (LACNT1)  


# the input is the output of findContact


