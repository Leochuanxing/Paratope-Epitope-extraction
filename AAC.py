# -*- coding:utf-8 -*-

from math import sqrt


def main(pdb_id='1bog', cutoff=6):
    with open(pdb_id + '.pdb', 'r') as f:
        pdb = f.readlines()

    chainL = findChain(pdb, 'A')
    chainH = findChain(pdb, 'B')
    chainA = findChain(pdb, 'C')

    cdr_l1, cdr_l2, cdr_l3 = findCDRL(chainL)
    cdr_h1, cdr_h2, cdr_h3 = findCDRH(chainH)

    cdr_l1_res = index2Residue(chainL, cdr_l1)
    cdr_l2_res = index2Residue(chainL, cdr_l2)
    cdr_l3_res = index2Residue(chainL, cdr_l3)

    cdr_h1_res = index2Residue(chainH, cdr_h1)
    cdr_h2_res = index2Residue(chainH, cdr_h2)
    cdr_h3_res = index2Residue(chainH, cdr_h3)

    cdnt_l1 = getCoordinates(pdb, 'A', cdr_l1)
    cdnt_l2 = getCoordinates(pdb, 'A', cdr_l2)
    cdnt_l3 = getCoordinates(pdb, 'A', cdr_l3)

    cdnt_h1 = getCoordinates(pdb, 'B', cdr_h1)
    cdnt_h2 = getCoordinates(pdb, 'B', cdr_h2)
    cdnt_h3 = getCoordinates(pdb, 'B', cdr_h3)

    cdnt_a = getCoordinates(pdb, 'C', list(range(len(chainA))))

    cnt_l1 = findContact(cdnt_l1, cdnt_a, cutoff)
    cnt_l2 = findContact(cdnt_l2, cdnt_a, cutoff)
    cnt_l3 = findContact(cdnt_l3, cdnt_a, cutoff)

    cnt_h1 = findContact(cdnt_h1, cdnt_a, cutoff)
    cnt_h2 = findContact(cdnt_h2, cdnt_a, cutoff)
    cnt_h3 = findContact(cdnt_h3, cdnt_a, cutoff)

    output = ParEpi(pdb_id, cdr_l1_res, cdr_l2_res, cdr_l3_res, cdr_h1_res, cdr_h2_res, cdr_h3_res, cnt_l1, cnt_l2,
                    cnt_l3, cnt_h1, cnt_h2, cnt_h3)
    return output


def findChain(pdb, chain_id):
    chain = []
    for line in pdb:
        if line[:6] == 'SEQRES' and line[11] == chain_id:
            chain.extend(line[19:70].strip().split(' '))
        elif len(chain) > 0:
            break
    return chain


def findCDRL(chain):
    # CDR-L1
    bias = 6  # position bias for "start: approx residue 24"
    cys_index = 0
    for i in range(23 - bias, 24 + bias):
        if chain[i] == 'CYS' and abs(i - 23) < abs(cys_index - 23):
            cys_index = i

    trp_index = 0
    tmp_trp_index = 0  # 低优先级
    for i in range(cys_index + 11, cys_index + 18):
        if chain[i] + chain[i + 1] + chain[i + 2] in ('TRPTYRGLN', 'TRPLEUGLN', 'TRPPHEGLN', 'TRPTYRLEU'):
            trp_index = i
        elif chain[i] == 'TRP' and chain[i + 13] + chain[i + 14] in ('ILETYR', 'VALTYR', 'ILELYS', 'ILEPHE'):
            tmp_trp_index = i
        elif tmp_trp_index == 0 and chain[i] == 'TRP':
            tmp_trp_index = i
    if trp_index == 0:
        trp_index = tmp_trp_index

    cdr_l1 = list(range(cys_index + 1, trp_index))

    # CDR-L2
    l2_start_index = trp_index + 15
    l2_end_index = l2_start_index + 6
    cdr_l2 = list(range(l2_start_index, l2_end_index + 1))

    # CDR-L3
    assert chain[l2_end_index + 32] == 'CYS', 'CDR-L3: Residue before should be CYS but not ' + chain[l2_end_index + 32]
    l3_start_index = l2_end_index + 33
    l3_end_index = 0
    for i in range(l3_start_index + 6, l3_start_index + 11):
        if chain[i + 1] + chain[i + 2] + chain[i + 4] == 'PHEGLYGLY':
            l3_end_index = i
    cdr_l3 = list(range(l3_start_index, l3_end_index + 1))

    return cdr_l1, cdr_l2, cdr_l3


def findCDRH(chain):
    # CDR-H1
    bias = 6  # position bias for "start: approx residue 26"
    cys_index = 0
    for i in range(21 - bias, 22 + bias):
        if chain[i] == 'CYS' and abs(i - 21) < abs(cys_index - 21):
            cys_index = i

    trp_index = 0
    tmp_trp_index = 0  # 低优先级
    for i in range(cys_index + 14, cys_index + 16):
        # print(chain[i] + chain[i + 1])
        if chain[i] + chain[i + 1] in ('TRPVAL', 'TRPILE', 'TRPALA'):
            trp_index = i
        elif chain[i] == 'TRP' and chain[i + 9] + chain[i + 10] + chain[i + 11] + chain[i + 12] + chain[
                    i + 13] == 'LEUGLUTRPILEGLY':
            tmp_trp_index = i
        elif tmp_trp_index == 0 and chain[i] == 'TRP':
            tmp_trp_index = i
    if trp_index == 0:
        trp_index = tmp_trp_index

    cdr_h1 = list(range(cys_index + 4, trp_index))

    # CDR-H2
    h2_start_index = trp_index + 14
    h2_end_index = 0
    for i in range(h2_start_index + 15, h2_start_index + 19):
        if chain[i + 1] in ('LYS', 'ARG') and chain[i + 2] in ('LEU', 'ILE', 'VAL', 'PHE', 'THR', 'ALA') and chain[
                    i + 3] in ('THR', 'SER', 'ILE', 'ALA'):
            h2_end_index = i
    cdr_h2 = list(range(h2_start_index, h2_end_index + 1))

    # CDR-H3
    assert chain[h2_end_index + 30] == 'CYS', 'CDR-H3: Residue before should be CYS-XXX-XXX but not ' + chain[
        h2_end_index + 30] + '-XXX-XXX'
    h3_start_index = h2_end_index + 33
    h3_end_index = 0
    for i in range(h3_start_index + 2, h3_start_index + 25):
        if chain[i + 1] + chain[i + 2] + chain[i + 4] == 'TRPGLYGLY':
            h3_end_index = i
    cdr_h3 = list(range(h3_start_index, h3_end_index + 1))

    return cdr_h1, cdr_h2, cdr_h3


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


def index2Residue(chain, residue_index):
    return [chain[i] for i in residue_index]


class ParEpi:
    def __init__(self, PdbId, CDRL1, CDRL2, CDRL3, CDRH1, CDRH2, CDRH3, L1CNT,
                 L2CNT, L3CNT, H1CNT, H2CNT, H3CNT):
        self.PdbId = PdbId
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
        return (('PDBid: {} \n CDRL1: {} \n CDRL2: {} \n CDRL3: {} \n' +
                 'CDRH1: {} \n CDRH2: {} \n CDRH3: {} \n L1CNT: {} \n' +
                 'L2CNT:{} \n L3CNT:{} \n H1CNT: {} \n H2CNT: {} \n' +
                 'H3CNT: {} \n')
                .format(self.PdbId, self.CDRL1, self.CDRL2, self.CDRL3, self.CDRH1,
                        self.CDRH2, self.CDRH3, self.L1CNT, self.L2CNT, self.L3CNT,
                        self.H1CNT, self.H2CNT, self.H3CNT))


if __name__ == '__main__':
    res = main('1bog', 6)
    print(res)
    # print(res.H3CNT)
