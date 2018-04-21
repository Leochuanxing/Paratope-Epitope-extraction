# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 19:59:59 2018

@author: leo
"""
from math import sqrt

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

def findContact(CDRcdnt, acdnt, cutoff):
    contactA = []
    for i in CDRcdnt:
        for j in acdnt:
            if sqrt((i[0] - j[0]) ** 2 + (i[1] - j[1]) ** 2 + (i[2] - j[2]) ** 2) <= cutoff:
                contactA.append((i[3], j[3]))
    contactB = list(set(contactA))
    contactB.sort(key=contactA.index)
    contactN = []
    for i in contactB:
        contactN.append(contactA.count(i))
    return contactB, contactN
# Structure of the output
class ParEpi:
    def __init__(self, PdbId, chainL, chainH, chainA, L1CNT,
                 L2CNT, L3CNT, H1CNT, H2CNT, H3CNT):
        self.PdbId = PdbId
        self.chainL = chainL
        self.chainH = chainH
        self.chainA = chainA
        self.L1CNT = L1CNT
        self.L2CNT = L2CNT
        self.L3CNT = L3CNT
        self.H1CNT = H1CNT
        self.H2CNT = H2CNT
        self.H3CNT = H3CNT

    def __str__(self):
        return (('PDBid: {} \n  chainL: {}\n  chainH: {}\n chainA: {} \n L1CNT: {} \n' +
                 'L2CNT:{} \n L3CNT:{} \n H1CNT: {} \n H2CNT: {} \n' +
                 'H3CNT: {} \n')
                .format(self.PdbId, self.chainL, self.chainH, self.chainA, self.L1CNT, self.L2CNT,
                        self.L3CNT,self.H1CNT, self.H2CNT, self.H3CNT))

# See http://www.biochem.ucl.ac.uk/~martin/abs/GeneralInfo.html
# We use the minimum of the starting position and the maximum of the ending position
# of Kabat, Chothia, AbM and Contact. For H3, we creat a range with length 20
CDRLindex = [list(range(23, 36)), list(range(45, 56)), list(range(88, 97))]
CDRHindex = [list(range(25, 36)), list(range(46, 65)), list(range(90, 110))]

def main(ID, cutoff=6):
    with open(ID[0] + '.pdb', 'r') as f:
        pdb = f.readlines()

    chainL = findChain(pdb, ID[2])
    chainH = findChain(pdb, ID[1])
    chainA = findChain(pdb, ID[3])

    cdnt_l1 = getCoordinates(pdb, ID[2], CDRLindex[0])
    cdnt_l2 = getCoordinates(pdb, ID[2], CDRLindex[1])
    cdnt_l3 = getCoordinates(pdb, ID[2], CDRLindex[2])

    cdnt_h1 = getCoordinates(pdb, ID[1], CDRHindex[0])
    cdnt_h2 = getCoordinates(pdb, ID[1], CDRHindex[1])
    cdnt_h3 = getCoordinates(pdb, ID[1], CDRHindex[2])

    cdnt_a = getCoordinates(pdb, ID[3], list(range(len(chainA))))

    cnt_l1 = findContact(cdnt_l1, cdnt_a, cutoff)
    cnt_l2 = findContact(cdnt_l2, cdnt_a, cutoff)
    cnt_l3 = findContact(cdnt_l3, cdnt_a, cutoff)

    cnt_h1 = findContact(cdnt_h1, cdnt_a, cutoff)
    cnt_h2 = findContact(cdnt_h2, cdnt_a, cutoff)
    cnt_h3 = findContact(cdnt_h3, cdnt_a, cutoff)
  
    Output = ParEpi(ID[0], chainL, chainH, chainA, cnt_l1, cnt_l2, cnt_l3, cnt_h1, cnt_h2, cnt_h3)
    return Output
  
    
summary = open('summary.TSV', 'r')
file = summary.readlines()        
summary.close
        
from IDhere import main
l = main (PDBfiles, ids)

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

for i in LAID:
    print(LACNT(i, 6))
 
print(kel)
problem = []
LACNT1 = []
for ID in kel:
    if ID[3] != '':
        try:
#    if __name__ == '__main__':
            Contact = main(ID, 6)
#            print([ID[1],ID[3]], [Contact.H1CNT, Contact.H2CNT, Contact.H3CNT],
            LACNT1.append([[ID[2],ID[3]], [Contact.L1CNT, Contact.L2CNT, Contact.L3CNT]])
        except:
            problem.append('Something is wrong with '+ ID[0]+'.pdb')
print(problem)
print (LACNT1)    












