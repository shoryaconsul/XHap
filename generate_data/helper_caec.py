# Helper file for CAECSeq_Haplo

import numpy as np
import sys
from itertools import permutations
import random

def import_config(configname):
    with open(sys.argv[1], 'r') as f:
        df = f.readlines()

    data = []
    for item in df:
        data.append(item.strip('\n').split(' : '))
        
    return data

# import SNV fragment matrix
def import_SNV(SNVmatrix_name):
    with open(SNVmatrix_name, 'r') as f:
        SNVmatrix_list = f.readlines()

    SNVmatrix = np.zeros((len(SNVmatrix_list), len(SNVmatrix_list[0][0:-1:2])))
    SNVonehot = np.zeros((SNVmatrix.shape[1], 4, SNVmatrix.shape[0]))

    for i in range(len(SNVmatrix_list)):
        SNVmatrix[i, :] = np.fromstring(SNVmatrix_list[i][0:-1:1], dtype = int, sep = ' ')
        for j in range(SNVmatrix.shape[1]):
            if SNVmatrix[i, j]:
                SNVonehot[j, int(SNVmatrix[i, j] - 1), i] = 1

    SNVonehot = np.moveaxis(SNVonehot, -1, 0)
    SNVonehot = np.moveaxis(SNVonehot, -1, 1)

    return SNVmatrix.astype(int), SNVonehot.reshape((SNVonehot.shape[0], SNVonehot.shape[1], SNVonehot.shape[2], 1))

from scipy.sparse import coo_matrix, csr_matrix, dok_matrix
def import_sparseSNV(mat_file):
    """
    Read read-SNP matrix into sparse matrix. Useful for large matrices.

    mat_file: str
        Path to file containing read-SNP matrix
    
    Returns
        NDArray: Sparse read-SNP matrix
    """

    with open(mat_file, "r") as f:
        vals, rows, cols = [], [], []
        idx_val_dict = {}
        nSNV = int(f.readline().strip())
        nReads = 0
        for line in f:
            line = line.strip()
            ind_vals = line.split()
            for iv in ind_vals:
                snv, val = iv.split(",")
                vals.append(int(val))
                rows.append(nReads)
                cols.append(int(snv))
                idx_val_dict[(nReads, int(snv))] = int(val)
            nReads = nReads + 1
    
    SNV_matrix = coo_matrix((vals, (rows, cols)), shape=(nReads, nSNV)).todok()
    
    SNV_onehot = np.zeros((SNV_matrix.shape[0], 4, SNV_matrix.shape[1], 1))
    nnz_idx = SNV_matrix.nonzero()
    for i, j in zip(nnz_idx[0], nnz_idx[1]):
        SNV_onehot[i, int(SNV_matrix[i, j] - 1), j, 0] = 1

    return SNV_matrix.toarray(), SNV_onehot


# get the ACGT statistics of a read matrix
def ACGT_count(submatrix):
    out = np.zeros((submatrix.shape[1], 4))
    for i in range(4):
        out[:, i] = (submatrix == (i + 1)).sum(axis = 0)

    return out 

# 
def origin2haplotype(origins, SNVmatrix, n_cluster):
    V_major = np.zeros((n_cluster, SNVmatrix.shape[1])) # majority voting result
    ACGTcount = ACGT_count(SNVmatrix)
    
    for i in range(n_cluster):         
        reads_single = SNVmatrix[origins == i, :] # all reads from one haplotypes
        single_sta = np.zeros((SNVmatrix.shape[1], 4))
        
        if len(reads_single) != 0:
            single_sta = ACGT_count(reads_single) # ACGT statistics of a single nucleotide position
        V_major[i, :] = np.argmax(single_sta, axis = 1) + 1            

        uncov_pos = np.where(np.sum(single_sta, axis = 1) == 0)[0]

        for j in range(len(uncov_pos)):
            if len(np.where(ACGTcount[uncov_pos[j], :] == max(ACGTcount[uncov_pos[j], :]))[0]) != 1: # if not covered, select the most doninant one based on 'ACGTcount'     
                tem = np.where(ACGTcount[uncov_pos[j], :] == max(ACGTcount[uncov_pos[j], :]))[0]
                V_major[i, uncov_pos[j]] = tem[int(np.floor(random.random() * len(tem)))] + 1
            else:
                V_major[i, uncov_pos[j]] = np.argmax(ACGTcount[uncov_pos[j], :]) + 1

    return V_major

# calculate hamming distance between two sequences
def hamming_distance(read, haplo):
    return sum((haplo - read)[np.where(read != 0)] != 0)

# calculate MEC
def MEC(SNVmatrix, Recovered_Haplo):
    res = 0
    
    for i in range(len(SNVmatrix)):
        dis = [hamming_distance(SNVmatrix[i, :], Recovered_Haplo[j, :]) for j in range(len(Recovered_Haplo))]
        res += min(dis)
        
    return res

def target_distribution_haplo(q, SNVmatrix, n_clusters):
    q = np.argmax(q, axis = 1)
    haplotypes = origin2haplotype(q, SNVmatrix, n_clusters)
    
    index = []
    for i in range(SNVmatrix.shape[0]):
        dis = np.zeros((haplotypes.shape[0]))
        for j in range(haplotypes.shape[0]):
            dis[j] = hamming_distance(SNVmatrix[i, :], haplotypes[j, :])
        index.append(np.argmin(dis))

    p = np.zeros((SNVmatrix.shape[0], n_clusters))
    for i in range(p.shape[0]):
        p[i, index[i]] = 1

    return p

# convert list of sequence to numpy array
def list2array(ViralSeq_list):
    ViralSeq = np.zeros((len(ViralSeq_list), len(ViralSeq_list[0])))

    for i in range(len(ViralSeq_list)):
        for j in range(len(ViralSeq_list[0])):
            if ViralSeq_list[i][j] == 'A':
                ViralSeq[i, j] = 1
            elif ViralSeq_list[i][j] == 'C':
                ViralSeq[i, j] = 2
            elif ViralSeq_list[i][j] == 'G':
                ViralSeq[i, j] = 3
            elif ViralSeq_list[i][j] == 'T':
                ViralSeq[i, j] = 4
            else:
            	ViralSeq[i, j] = 0
    
    return ViralSeq

# evaluate the recall rate and reconstruction rate
def recall_reconstruction_rate(Recovered_Haplo, SNVHaplo):
    distance_table = np.zeros((len(Recovered_Haplo), len(SNVHaplo)))
    for i in range(len(Recovered_Haplo)):
        for j in range(len(SNVHaplo)):
            distance_table[i, j] = hamming_distance(SNVHaplo[j, :], Recovered_Haplo[i, :])
    
    index = list(permutations(list(range(SNVHaplo.shape[0]))))
    distance = []
    for item in index:
        count = 0
        for i in range(len(item)):
            count += distance_table[i, item[i]]
        distance.append(count)
    index = index[np.argmin(np.array(distance))]
    
    reconstruction_rate = []
    for i in range(len(index)):
        reconstruction_rate.append(1 - distance_table[i, index[i]] / SNVHaplo.shape[1])
    
    recall_rate = reconstruction_rate.count(1) / len(reconstruction_rate)
    
    CPR = 1 - min(distance) / (len(distance_table) * SNVHaplo.shape[1])
    
    return distance_table, reconstruction_rate, recall_rate, index, CPR


def best_match(rec_hap, new_hap):
	if np.shape(rec_hap) != np.shape(new_hap):
		raise ValueError("Input arguments should have the same shape.")

	distance_table = np.zeros((len(rec_hap), len(new_hap)))
	for i, rh in enumerate(rec_hap):
		for j, nh in enumerate(new_hap):
			distance_table[i, j] = hamming_distance(rh, nh)

	index = permutations(range(new_hap.shape[0]))
	min_distance = np.inf 
	distance = []
	for matching in index:
		count = 0
		for i, match_idx in enumerate(matching):
			count += distance_table[i, match_idx]
		distance.append(count)
		if count < min_distance:
			best_matching = matching
			min_distance = count

	return best_matching