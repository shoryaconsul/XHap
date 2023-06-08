import numpy as np
from numpy import random as rn
from scipy.spatial.distance import pdist, squareform
from os import path
from copy import deepcopy
from tqdm import tqdm
import time
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

import logging
from typing import Any, List, Tuple, Optional
from nptyping import NDArray

import torch
from torch import nn
from torch import optim
from torch.utils.data import Dataset, DataLoader
from torch.nn import functional as F

from matplotlib import pyplot as plt

from read_embeddings import save_ckp, load_ckp, MyFilter, \
    ReadAE, SNVMatrixDataset, learn_embed
from xhap_parallel import CorrTransformer
from helper import * # read_hap, read_true_hap, compute_cpr


def readW(SNVdataset: SNVMatrixDataset, 
          ae: ReadAE, xformer: CorrTransformer,
          device: torch.cuda.device=torch.device("cpu")
         ):
    
    m = len(SNVdataset)
    dataloader_full = DataLoader(SNVdataset, batch_size=m,
                                num_workers=0)
    for i, (data, idx) in enumerate(dataloader_full):
        SNV_onehot = data.to(device)
    
    with torch.no_grad():  # Do not compute grads for speedup    
        ae.eval()  # Set eval flags
        xformer.eval()  
    
    embed, _ = ae(SNV_onehot)  # Read embeddings
    W_full = xformer(embed[None,:]).cpu().detach().numpy()  # Converting to numpy
    return W_full[0]

def order_reads(SNV_matrix, read_origin):
    SNP_cover = SNV_matrix > 0
    nreads, nSNP = SNV_matrix.shape  # Number of reads
    read_idx = np.arange(nreads)
    
    if nreads == 1:
        return read_idx

    return np.lexsort(SNP_cover.T[::-1])
    # return np.lexsort(tuple(
    #     [SNP_cover[:, i] for i in range(nSNP-1,0,-1)] + [read_origin]
    #     ))
    # return np.lexsort(SNP_cover)

# Parameters
outhead = 'hap_n2_cov30_soltub'
d_model = 128  # Size of embeddings
device = torch.device("cuda")
# -----------------------------------

datapath = 'generate_data/' + outhead + '/' + outhead + '_SNV_matrix.txt'  
gt_file = 'generate_data/' + outhead + '/combined.fa'
pos_file = 'generate_data/' + outhead + '/' + outhead + '_SNV_pos.txt'
model_file = 'generate_data/' + outhead + '/xhap_model'

SNVdata = SNVMatrixDataset(datapath)
SNV_matrix = np.loadtxt(datapath, dtype=int)
SNV_matrix = SNV_matrix[np.sum(SNV_matrix != 0, axis=1) > 1]  # Removing uninformative reads
print('SNP matrix: ', SNV_matrix.shape)

nSNP = SNV_matrix.shape[1] # Number of SNVs
num_read = len(SNVdata)  # Number of reads
batch_size = int(np.ceil(num_read/5))

dataloader = DataLoader(SNVdata, batch_size=batch_size,
                        shuffle=True, num_workers=0)    

# Instantiate read embedding encoder
savefile="read_AE"
embedAE = learn_embed(SNVdata, num_epoch=0, embed_dim=d_model, savefile=savefile)
corr_xformer = CorrTransformer(d_model, d_model//2)  # Transformer

ckp = torch.load(model_file)
embedAE.load_state_dict(ckp['embed_ae'])
corr_xformer.load_state_dict(ckp['xformer'])

embedAE.to(device)
corr_xformer.to(device)

# Read transformer-inferred W matrix
W_best = readW(SNVdata, embedAE, corr_xformer, device)

# Measured/observed read correlation matrix
num_reads = np.shape(SNV_matrix)[0]
W_sim = np.zeros((num_reads, num_reads))
W_dissim = np.zeros((num_reads, num_reads))
W_mask = np.zeros((num_reads, num_reads), dtype=bool)

# Computing similarity matrix for supervision
for i, read_i in enumerate(tqdm(SNV_matrix)):
    for j, read_j in enumerate(SNV_matrix):
        if np.any((read_i != 0) & (read_j != 0)):  # Only if reads overlap
            W_mask[i, j] = True
            W_sim[i, j] = np.sum((read_i == read_j)[(read_i != 0) & (read_j != 0)])
            W_dissim[i, j] = np.sum((read_i != read_j)[(read_i != 0) & (read_j != 0)])

W_mask = torch.from_numpy(W_mask)
W_sup = (W_sim - W_dissim)/(W_sim + W_dissim + 1e-10)
np.fill_diagonal(W_sup, 1.)

# Read origins from W_Sup
true_hap = read_true_hap(gt_file, pos_file)
read_origin = []
for SNV_read in SNV_matrix:
    dist_list = [hamming_distance(SNV_read, hap) for j, hap in enumerate(true_hap)]
    read_origin.append(np.argmin(dist_list))
read_origin = np.array(read_origin)

# Sort correlation matrices
# idx = order_reads(W_sup)
# idx = order_reads(SNV_matrix, read_origin)
linkage_matrix = hierarchy.linkage(squareform(1 - W_sup), method='single')
idx = hierarchy.leaves_list(linkage_matrix)

W_sup_sort = W_sup[idx, :]
W_sup_sort = W_sup_sort[:, idx]

W_best_sort = W_best[idx, :]
W_best_sort = W_best_sort[:, idx]

# Create heatmap
fig, ax = plt.subplots(1,2, figsize=(16,8))
im1 = ax[0].imshow(W_sup_sort, cmap='coolwarm')
im2 = ax[1].imshow(W_best_sort, cmap='coolwarm')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.83, 0.15, 0.05, 0.7])
fig.colorbar(im1, cax=cbar_ax)
plt.savefig('test_plots_1.jpg', bbox_inches='tight')
# plt.savefig('generate_data/' + outhead + '/corr_plots.jpg', bbox_inches='tight')

corr_thresh = 0.1
r_nz, c_nz = np.nonzero(np.abs(W_best_sort - W_sup_sort)*(W_sup_sort == 0)
                        > corr_thresh)
nz_sort = np.argsort(np.abs(c_nz - r_nz))[::-1]
print(np.sort(np.abs(c_nz - r_nz))[::-1])
print(W_best_sort[r_nz[nz_sort], c_nz[nz_sort]])
print("Number of discovered correlations %d, Discovered fraction of correlations %.3f" 
      %(len(r_nz), len(r_nz)/np.sum(W_sup_sort == 0))
      )
# idx_nz = np.argsort(np.abs(c_nz - r_nz))[::-1][:100]
# print(list(zip(r_nz[idx_nz], c_nz[idx_nz])))

# print('CHECKING READS')
# with open(pos_file, "r") as f:
# 		pos_str = f.readline().split()
# 		pos = np.array([int(ps) - 1 for ps in pos_str])

# read_gap_list = []
# read_dist_list = []
# for r, c in zip(r_nz, c_nz):
#     r1_snp = np.nonzero(SNV_matrix[idx][r])[0]
#     r2_snp = np.nonzero(SNV_matrix[idx][c])[0]
#     if r1_snp[0] > r2_snp[0]:
#         r1_snp, r2_snp = r2_snp, r1_snp
#     read_gap_list.append(pos[r2_snp[0]] - pos[r1_snp[-1]])
#     read_dist_list.append(pos[r2_snp[0]] - pos[r1_snp[0]])
# print('MAX READ GAP:')
# print(np.amax(read_gap_list))
# print('MAX READ DISTANCE:')
# print(np.amax(read_dist_list))