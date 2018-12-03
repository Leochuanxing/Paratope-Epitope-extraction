import os
os.getcwd()
os.chdir("/home/leo/Documents/Database/Pipeline/Homo")

import json
with open('sequence', 'r') as f:
    sequence = json.load(f)
with open('good_matched_ids', 'r') as f:
    good_matched_ids = json.load(f)
with open('dud_AAC', 'r') as f:
    dud_AAC = json.load(f)
with open('contact', 'r') as f:
    contact = json.load(f)
# Get rid of the dud
for dud in dud_AAC:
    del good_matched_ids[dud]
good_matched_ids.keys()
good_matched_ids['4om1']
contact.keys()
contact['1adq']
sequence['1adq'].keys()
len(sequence['1adq']['H'])
#################################################################################



##########################################################################
'''
Write the above functions into a class, named 'AlignmentConstraint'

'''



class AlignmentConstraint(object):
    
    def __init__(self, iddict, sequence, contact):
        self.iddict = iddict
        self.sequence = sequence
        self.contact = contact
        self.TripleSingle =  [['TYR', 'Y'], ['LYS', 'K'],['ASP', 'D'], ['ASN', 'N'], ['TRP', 'W'], ['PHE', 'F'], ['GLN', 'Q'],
       ['GLU', 'E'], ['PRO', 'P'], ['GLY', 'G'], ['THR', 'T'],['SER', 'S'], ['ARG', 'R'], ['HIS', 'H'],
       ['LEU', 'L'], ['ILE', 'I'], ['CYS', 'C'], ['ALA', 'A'], ['MET', 'M'], ['VAL', 'V']]
               
    
    def To_seq(self, aa_sequence):
        from Bio.Seq import Seq
        from Bio.Alphabet import IUPAC
        
        seq_obj = None

        seq_single_letter = ''
        for aa in aa_sequence:
            for TS in self.TripleSingle:
                if TS[0] == aa:
                    seq_single_letter += TS[1]
        seq_obj = Seq(seq_single_letter, IUPAC.protein)
        
        return seq_obj
    '''
    SeqCDR:
        to convert the sequence into sequence object and be prepared for futher alignment
    Input:
        The self object
    Output:
        in the form of a list [['1abvhH', Seq_obj], ['1abvlL', Seq_obj]]
        ['1abvhH', Seq_obj] means the pdbid is 1abv heavy chain with chain name H, and 
        the Seq_obj is the sequence object of the concatenated CDR1, CDR2, CDR3 of chain H.
    '''   
    def SeqCDR(self):
        # Concatenate the sequences 
        seqCDRH = []
        seqCDRL = []
        for pdbid in self.iddict:
            for combination in self.iddict[pdbid]:
                if combination[2] != '':
                    if combination[0] != '':
                        CDRH = [pdbid+'h'+combination[0]]
                        CDRH.append(self.To_seq(self.sequence[pdbid][combination[0]][25:36]))
                        CDRH.append(self.To_seq(self.sequence[pdbid][combination[0]][46:65]))
                        CDRH.append(self.To_seq(self.sequence[pdbid][combination[0]][90:110]))
                        seqCDRH.append(CDRH)
                    if combination[1] != '':
                        CDRL = [pdbid+'l'+combination[1]]
                        CDRL.append(self.To_seq(self.sequence[pdbid][combination[1]][23:36]))
                        CDRL.append(self.To_seq(self.sequence[pdbid][combination[1]][45:56]))
                        CDRL.append(self.To_seq(self.sequence[pdbid][combination[1]][88:97]))
                        seqCDRL.append(CDRL)                        
        self.seqCDRH = seqCDRH
        self.seqCDRL = seqCDRL
        
    '''
    Align according to the given criteria and choose the none redundent one
    '''
    
    def Hamming_like_dist (u, v):
        # 1 creat the similarity matrix
        # 2 Do hyerachical clustering
        # 3 Chose the representative
        dist = 0 
        from Bio import Align
        aligner = Align.PairwiseAligner()
        aligner.open_gap_score = -1 # it is unlikly that two sequence of the same source could be different by a gap
        aligner.extend_gap_score = -1
        aligner.match = 1
        aligner.mismatch = 0 
        aligner.mode = 'local'
        
        dist += aligner.score(u[1], v[1])
        dist += aligner.score(u[2], v[2])
        dist += aligner.score(u[3], v[3])
        l1 = len(u[1]) + len(u[2]) + len(u[3])
        l2 = len(v[1]) + len(v[2]) + len(v[3])
        l = min(l1, l2)
        dist =  1 - dist /l
        return dist
    
    '''
    Input: self
    Output: hcluster_CDRL, hcluster_CDRH. They are given in the form of
             [[0, 1, 0.5, 2],
              [2, 3, 0.5, 2],
              [4, 5, 1, 4]]
        each row is in the form of [idx1, idx2, dist, sample_count]
        In the above example, there are for samples with index 0, 1, 2, 3.
        4 means the cluster formed in row number 4 - len(samples) = 0
        5 means the cluster formed in row number 5 - len(samples) = 1.
        Therefore, the last row means combine the clusters formed in the above 
        two rows, and the sample count is 4, all the samples.
    
    '''
    
    def Hcluster(self):

        from Bio import Align
        aligner = Align.PairwiseAligner()
        aligner.open_gap_score = -1 # it is unlikly that two sequence of the same source could be different by a gap
        aligner.extend_gap_score = -1
        aligner.match = 1
        aligner.mismatch = 0 
        aligner.mode = 'local'
        
        from scipy.cluster.hierarchy import linkage, optimal_leaf_ordering
        from scipy.spatial.distance import squareform
        import numpy as np
        
        CDR = [self.seqCDRH, self.seqCDRL]
        
        for idx in range(2):
            
            toydata = CDR[idx]
            distmatrix = np.zeros((len(toydata), len(toydata)))
            for i in range(len(toydata)):
                for j in range(len(toydata)):
                    distmatrix[i,j] += aligner.score(toydata[i][1], toydata[j][1])
                    distmatrix[i,j] += aligner.score(toydata[i][2], toydata[j][2])
                    distmatrix[i,j] += aligner.score(toydata[i][3], toydata[j][3])
                    l1 = len(toydata[i][1]) + len(toydata[i][2]) + len(toydata[i][3])
                    l2 = len(toydata[j][1]) + len(toydata[j][2]) + len(toydata[j][3])
                    l = min(l1, l2)
                    distmatrix[i,j] =  1 - distmatrix[i,j]/l
                    
            if idx == 0:
                self.distmatrix_CDRH = distmatrix
                dm_CDRH = squareform(distmatrix)
                hcluster_CDRH_unordered = linkage(dm_CDRH, method = 'complete')
                self.hcluster_CDRH = optimal_leaf_ordering(hcluster_CDRH_unordered, dm_CDRH )
            else:
                self.distmatrix_CDRL = distmatrix
                dm_CDRL = squareform(distmatrix)
                hcluster_CDRL_unordered = linkage(dm_CDRL, method = 'complete')
                self.hcluster_CDRL = optimal_leaf_ordering(hcluster_CDRL_unordered, dm_CDRL )
                
    def Show_elbow(self):
        
        from scipy.cluster.hierarchy import cut_tree
        
        n_clusters_CDRH = []
        n_clusters_CDRL = []
        heights = [x * 0.02 for x in range(51)]
        for h in heights:
            cut_CDRH = cut_tree(self.hcluster_CDRH, height = h)
            cut_CDRL = cut_tree(self.hcluster_CDRL, height = h)
            maxH = -1
            maxL = -1
            for i in range(len(cut_CDRH)):
                if cut_CDRH[i][0] >= maxH:
                    maxH = cut_CDRH[i][0]
            n_clusters_CDRH.append(maxH)
            
            for i in range(len(cut_CDRL)):
                if cut_CDRL[i][0] >= maxL:
                    maxL = cut_CDRL[i][0]
            n_clusters_CDRL.append(maxL)
            
        from matplotlib import pyplot as plt
        
        plt.title('CDRH_clusters Vs cut distance')
        plt.xlabel('cut distance')
        plt.ylabel('CDRH_clusters')
        plt.plot(heights, n_clusters_CDRH)
        plt.show()
        plt.close()
       
        
        plt.title('CDRL_clusters Vs cut distance')
        plt.xlabel('cut distance')
        plt.ylabel('CDRL_clusters')
        plt.plot( heights, n_clusters_CDRL)
        plt.show()
        plt.close()
        
        
        
    '''
    Choose a representative from each cluster with the highest contact number to 
    represent the cluster
    '''              
    def Cut_Select(self, CDRH_cut_dist, CDRL_cut_dist):
        # First we cut the tree, then change the number to ids
        # Second, find the representative with the highest contact numbers 
        from scipy.cluster.hierarchy import cut_tree
        
        samples = [self.seqCDRH, self.seqCDRL]
        for indicator in range(2):
            if indicator == 0:
                cut1 = cut_tree(self.hcluster_CDRH, height = CDRH_cut_dist)
            else:
                cut1 = cut_tree(self.hcluster_CDRL, height = CDRL_cut_dist) 
                
        #find the corresponding ids
            n = 0
            for i in range(len(cut1)):
                if cut1[i][0] >= n:
                    n = cut1[i][0]
    
            clusters = []
            for i in range(n+1):
                subcluster = []
                for j in range(len(cut1)):
                    if cut1[j][0] == i:
                        subcluster.append(j)
                clusters.append(subcluster)
            # select the representatives 
            rpts = []
            m_ctns = []
            for cluster in clusters:
                representative = -1
                max_contact = -1
                for i in cluster:
                    contact_number = 0
                    pdbid = samples[indicator][i][0][:4]
                    hl = samples[indicator][i][0][4]
                    HL = samples[indicator][i][0][5]
                    for fcdn in self.contact[pdbid]:
                        if fcdn[0][0] == hl and fcdn[0][2] == HL:
                            contact_number += fcdn[3]
                            
                    if contact_number >= max_contact:
                        max_contact = contact_number
                        representative = i
                rpts.append(samples[indicator][representative][0])
                m_ctns.append(max_contact)
            #load
            if indicator == 0:
                self.CDRH_rpts = rpts
                self.CDRH_max_ctns = m_ctns
            else:
                self.CDRL_rpts = rpts
                self.CDRL_max_ctns = m_ctns
                
    def Prepare_for_FrameConstraint(self):
        ac_contact = {}
        all_rep = self.CDRH_rpts[:]
        all_rep.extend(self.CDRL_rpts)
        for rpt in all_rep:
            container = []
            pdbid = rpt[:4]
            hl = rpt[4]
            HL = rpt[5]
            for fcdn in self.contact[pdbid]:
                if fcdn[0][0] == hl and fcdn[0][2] == HL:
                    container.append(fcdn)
            ac_contact[pdbid + HL + hl] = container
        self.ready_for_FrameConstraint = ac_contact
        
    def Save(self):
        import json
        with open('ready_for_FrameConstraint', 'w') as f:
            json.dump(self.ready_for_FrameConstraint, f)
                

        
        


                
        
 



alignment =  AlignmentConstraint(good_matched_ids, sequence, contact)
alignment.SeqCDR()
alignment.Hcluster()
alignment.Show_elbow()
alignment.Cut_Select(CDRH_cut_dist = 0.1, CDRL_cut_dist = 0.1)
alignment.CDRH_max_ctns
alignment.CDRH_rpts
len(alignment.CDRH_rpts)
len(alignment.CDRH_max_ctns)
len(alignment.CDRL_rpts)
len(alignment.CDRL_max_ctns)
alignment.Prepare_for_FrameConstraint()
alignment.ready_for_FrameConstraint.keys()
alignment.ready_for_FrameConstraint['5tzuLl']
len(alignment.ready_for_FrameConstraint)
alignment.CDRH_rpts[:6]
alignment.CDRL_rpts[:6]
alignment.Save()



#from scipy.cluster.hierarchy import cut_tree
#samples = alignment.seqCDRH
#contacts = alignment.contact
#cut1 = cut_tree(alignment.hcluster_CDRH, height = 0.1)
#n = 0
#for i in range(len(cut1)):
#    if cut1[i][0] >= n:
#        n = cut1[i][0]
#
#clusters = []
#for i in range(n+1):
#    subcluster = []
#    for j in range(len(cut1)):
#        if cut1[j][0] == i:
#            subcluster.append(j)
#    clusters.append(subcluster)






'''
Input: 
     clusters: the clusters, in the form of [[0], [1, 2], [3, 4, 9, 10]], the numbers are 
     the indices of the samples
     distmatrix: the distance matrix used to generate the hierarchical clustering tree
Output:
     max_inner: a float number, gives the smallest value among all the inner distances
     of all the clusters
     min_inter: gives the smallest distances among all the distances between clusters
'''    
def Proof_read_clusters(clusters):
    import numpy as np
    
# to check whether the clusters are correct
    innerdistance = [] 
    for cluster in clusters:
        score = 0 
        for i in cluster:
            for j in cluster:
                if alignment.distmatrix_CDRH[i, j] >= score:
                    score = alignment.distmatrix_CDRH[i, j]
        innerdistance.append(score)
    
    max_inner = max(innerdistance)
   
    min_inter = 2
    interdistances = np.zeros((len(clusters), len(clusters)))
    for m in range(len(clusters)):
        for n in range(len(clusters)):
            score = 0
            for i in clusters[m]:
                for j in clusters[n]:
                    if alignment.distmatrix_CDRH[i, j] >= score:
                        score = alignment.distmatrix_CDRH[i, j]
            interdistances[m, n] = score
            if m != n and score <= min_inter:
                min_inter= score
    
    return max_inner, min_inter






          



#import numpy as np
#z_test = np.array([[0, 1, 0.2, 2],
#              [2, 3, 0.1, 2],
#              [5, 6, 0.4, 4],
#              [7, 4, 0.6, 5]])
#        




'''
Hcluster: A function to perform hierarchical clustering
Input: distance_matrix, a symmetric matrix, with the main diagal 0
Output: tree, in the form of an array, which could be used to show the 
        dendrogram.
             [[0, 1, 0.5, 2],
              [2, 3, 0.5, 2],
              [4, 5, 1, 4]]
        each row is in the form of [idx1, idx2, dist, sample_count]
        In the above example, there are for samples with index 0, 1, 2, 3.
        4 means the cluster formed in row number 4 - len(samples) = 0
        5 means the cluster formed in row number 5 - len(samples) = 1.
        Therefore, the last row means combine the clusters formed in the above 
        two rows, and the sample count is 4, all the samples.
'''

    



