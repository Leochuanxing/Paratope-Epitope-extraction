# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 15:15:42 2018

@author: leo
"""


import os
os.getcwd()
os.chdir("C:/Users/leo/Documents/Database/PDB Learning")
os.getcwd()

'''Find the Complex Ids and Chain Ids'''
def FindIds (file):
    Pdbid = []
    Ids = []
    for l in file:
        i = 0
        j = 0
        k = 0
        Id = []
        if len(l) <= 10:
            print(l)
        else:
            for a in l:
                j = j + 1
                if a == '\t' and k <= 4:
                    k = k + 1
                    Id.append(l[i : j - 1])
                    i = j     
            if len(Id) == 5:
                Id[3] = ''.join(Id[4].split('|'))            
                if Id[0] not in Pdbid:
                    Pdbid.append(Id[0])
                    Ids.append([Id[0], Id[1], Id[2], Id[3]])
    return Ids

'''Find the matched files，in the return values, the first one is the pdbid, the se
    second one is the heavy chain, the third one is the light chain and the last 
    one is the antigen. There may be more than on chain id for the antigen.'''
def MatchedFiles (ListFiles, Summary):
    MatchedFiles = []
    for f in ListFiles:
        if len(f) == 8 and f[5:8] == 'pdb':
            for Id in Summary:
                if f[:4] == Id[0]:
                   MatchedFiles.append([f, Id[1], Id[2], Id[3]])
    return MatchedFiles

ListFiles = os.listdir()
file = open('summary.tsv', 'r')
summary = file.readlines()
file.close
Ids = MatchedFiles(ListFiles,FindIds(summary))
print(Ids)

'''Find the CDRs of the lightchains'''
def CDRL (lightchain):
#    motif_cdr_l1_b = ['CYS']
#    motif_cdr_l1_a = (['TRP', ['TRP', 'TYR', 'GLN'], ['TRP', 'LEU', 'GLN'], ['TRP', 'PHE', 'GLN'], ['TRP', 'TYR', 'LEU']])
#    motif_cdr_l2_b = [['ILE', 'TYR'], ['VAL', 'TYR'], ['ILE', 'LYS'], ['ILE', 'PHE']]
#    motif_cdr_l3_b = ['CYS']
#    motif_cdr_l3_A = ['PHE', 'GLY', 'XXX', 'GLY']
    motif_cdr_l3_A = ['PHE', 'GLY',  'GLY']  
#    find cdrl1
    a = []
    i = 18
    while i <= 29:
        if lightchain[i] == 'CYS':
            a.append(i)
        i = i+ 1
    b = a[0]
    for i in a:
        if abs(i - 23) < abs(b - 23):
            b = i
    A = b + 10
    while lightchain[A] != 'TRP' and A <= b + 17:
        A = A + 1
    cdrl1 = lightchain[b+1 : A]   
    l1 = list(range(b+2, A+1))
    """find cdrl2"""
    cdrl2 = lightchain[A+15 : A+22] 
    l2 = list(range(A + 16, A + 23))
    """find cdrl3"""
    B = A + 22 + 30 
    j = B
    while j <= B + 22 + 36:
        if lightchain[j] != 'CYS':
            j = j + 1
        else:
            break
    i = j + 6    
    while i <= len(lightchain) - 4:
        if [lightchain[i], lightchain[i+1], lightchain[i+3] ] != motif_cdr_l3_A:
            i = i + 1
        else:
            break
    cdrl3 = lightchain[j+1:i]
    l3 = list(range(j+2, i+1))
    return cdrl1, cdrl2, cdrl3, l1, l2, l3
   
'''Find the CDRs of the heavychains'''
def CDRH(heavychain):
    """find cdrh1"""
    i = 21
    while i <= 31:
        if heavychain[i] != 'CYS':
            i = i + 1
        else:
            break
    j = i + 3 + 10
    while j <= i + 3 + 12:
        if heavychain[j] != 'TRP':
            j = j + 1
        else:
            break
    cdrh1 = heavychain[i+3+1: j]
    h1 = list(range(i+5, j+1))
    """find cdrh2"""
    end1 = [['LYS', 'ARG'],['LEU','ILE', 'VAL', 'PHE', 'THR', 'ALA'],['THR', 'SER', 'ILE', 'ALA']]
    end2 = []

    for a in end1[0]:
        for b in end1[1]:
            for c in end1[2]:
                end = []
                end.append(a)
                end.append(b)
                end.append(c)
                end2.append(end)
    j = j -1 + 15
    k = j + 10
    while k <= j + 24:
        flag = 1
        for a in end2:
            if [heavychain[k], heavychain[k+1], heavychain[k+2]] == a:
                flag = 0
                break
        if flag == 0:
            break
        else:
            k = k+1
    cdrh2 = heavychain[j:k]
    h2 = list(range(j+1, k+1))
    """find cdrh3"""
    k = k - 1 + 33
    l = k
    while l <= len(heavychain)-4:
        if [heavychain[l],heavychain[l+1],heavychain[l+3]] != ['TRP', 'GLY', 'GLY']:
            l = l + 1
        else:
            break
    cdrh3 = heavychain[k:l]
    h3 = list(range(k+1, l+1))
    return cdrh1, cdrh2, cdrh3, h1, h2, h3


"""Find the chain sequence with given chainId"""
def Chain_finder(pdb, chainId):
    chain = []
    for line in pdb:
        if line[:6] == 'SEQRES' and line[11] == chainId: 
            i = 19
            while i >= 19 and i <= 69:
                if line[i] != ' ':
                    chain.append(line[i] + line[i+1] + line[i+2]) 
                    i = i + 3
                else:
                    i = i + 1 
    return chain

"""Find the coordinates for interval with chainID"""
def Coordinates(pdb, chainId, Interval):
    CDNT = []
    m = 0
    i = 0
    for line in pdb:
        if line[:4] == 'ATOM' and line[21] == chainId:
            n = int(line[22:26])
            if m != n:
                m = n
                i = i + 1
            if i in Interval:
                CDNT.append([float(line[30:38]),float(line[38:46]), float(line[46:54]), i])
    return CDNT
          
"""Find the contact sequence in chain  for a given seq with cutoff """
def FindContact(CDR, acdnt, cutoff): 
    contactA = []
    for c1 in CDR:
        for c2 in acdnt:
            d = ((c1[0]-c2[0])**2+(c1[1]-c2[1])**2+(c1[2]-c2[2])**2)**0.5
            if d <= cutoff:
                contactA.append([c1[3], c2[3]])
    '''reduce redundancy and calculate contact numbers'''
    if contactA !=[]:
        contactB = [contactA[0],]
        contactN = []
        m = 0
        for i in contactA:
            if i  not in contactB:
                contactB.append(i) 
                
        for j in contactB:
            for k in contactA:
                if k == j:
                    m = m + 1
            contactN.append(m)
            m = 0  
    else:
        contactB = []
        contactN = []
    return contactB, contactN
                
    
class ParEpi:
    def __init__(self, PdbId, Lchain, Hchain, Achain, CDRL1, CDRL2, CDRL3,
                 CDRH1, CDRH2, CDRH3, L1CNT, L2CNT, L3CNT, H1CNT, H2CNT, H3CNT ):
        self.PdbId = PdbId
        self.Lchain = Lchain
        self.Hchain = Hchain
        self.Achain = Achain
        self.CDRL1 = CDRL1
        self.CDRL2 = CDRL2
        self.CDRL3 = CDRL3
        self.CDRH1 = CDRH1
        self.CDRH2 = CDRH2
        self.CDRH3 = CDRH3
        self.L1CNT = L1CNT
        self.L2CNT = L2CNT
        self.L3CNT = L3CNT
        self.H1CNT = H1CNT
        self.H2CNT = H2CNT
        self.H3CNT = H3CNT
    
    def __str__(self):
        return (('PDBid: {} \n Lchain: {} \n Hchain: {} \n Achian: {} \n' +
                 'CDRL1: {} \n CDRL2: {} \n CDRL3: {} \n' +
                 'CDRH1: {} \n CDRH2: {} \n CDRH3: {} \n L1CNT: {} \n' +
                 'L2CNT:{} \n L3CNT:{} \n H1CNT: {} \n H2CNT: {} \n' +
                 'H3CNT: {} \n')
                .format(self.PdbId, self.Lchain, self.Hchain, self.Achain, 
                        self.CDRL1, self.CDRL2,self.CDRL3, self.CDRH1,
                       self.CDRH2, self.CDRH3, self.L1CNT, self.L2CNT, self.L3CNT,
                       self.H1CNT, self.H2CNT, self.H3CNT ))             
        

'''Read the file'''
file = open(Ids[0][0], 'r')
pdb = file.readlines()
file.close

'''Get the chains and CDRs'''
Hchain = Chain_finder(pdb,Ids[0][1])
Lchain = Chain_finder(pdb, Ids[0][2])
Achain = Chain_finder(pdb, Ids[0][3])

cdrl = CDRL(Lchain)
CDRL1 = cdrl[0]
CDRL2 = cdrl[1]
CDRL3 = cdrl[2]

cdrh = CDRH(Hchain)
CDRH1 = cdrh[0]
CDRH2 = cdrh[1]
CDRH3 = cdrh[2]

'''Get the coordinates'''

h1cdnt = Coordinates(pdb, Ids[0][1], cdrh[3])
h2cdnt = Coordinates(pdb, Ids[0][1], cdrh[4])
h3cdnt = Coordinates(pdb, Ids[0][1], cdrh[5])

l1cdnt = Coordinates(pdb, Ids[0][2], cdrl[3])
l2cdnt = Coordinates(pdb, Ids[0][2], cdrl[4])
l3cdnt = Coordinates(pdb, Ids[0][2], cdrl[5])

acdnt = Coordinates(pdb,'A', list(range(1, len(Achain)+1)))

'''Get the contact'''
cutoff = 6
L1CNT = FindContact(l1cdnt, acdnt, cutoff)
L2CNT = FindContact(l2cdnt, acdnt, cutoff)
L3CNT = FindContact(l3cdnt, acdnt, cutoff)

H1CNT = FindContact(h1cdnt, acdnt, cutoff)
H2CNT = FindContact(h2cdnt, acdnt, cutoff)
H3CNT = FindContact(h3cdnt, acdnt, cutoff)

'''Output'''

output = ParEpi(Ids[0][0], Lchain, Hchain, Achain, CDRL1, CDRL2, CDRL3, CDRH1, CDRH2, CDRH3, 
                L1CNT,L2CNT, L3CNT, H1CNT, H2CNT, H3CNT)
print(output)


