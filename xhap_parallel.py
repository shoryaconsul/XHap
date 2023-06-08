import numpy as np
from numpy import random as rn
from scipy.spatial.distance import pdist, squareform
from os import path
from copy import deepcopy
from tqdm import tqdm
import time

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
# from ScheduleOptim import ScheduledOptim
from helper import * # read_hap, read_true_hap, compute_cpr
# from kernel_kmeans import KernelKMeans
from kernel_kmeans_torch import KernelKMeans

class CorrTransformer(nn.Module):
    def __init__(self, d_model: int=512, d_qk: int=None):
        super().__init__()
        d_ff = 4*d_model  # Dimension of intermediate feedforward layer
        n_head = 4  # Number of attention heads
        if d_qk is None:
            d_qk = int(d_model//n_head)  # Dim of query and key vectors in last attn layer

        self.d_model = d_model
        self.d_qk = d_qk
        encoderLayer = nn.TransformerEncoderLayer(d_model, n_head, d_ff,
                                                  batch_first=True)
        self.encoder = nn.TransformerEncoder(encoderLayer, num_layers=3)
        self.attn_lastq = nn.Linear(d_model, d_qk)
#         self.attn_lastk = nn.Linear(d_model, d_qk)
#         self.activation = nn.Tanh()
#         self.activation = nn.Softmax(dim=-1)

    def forward(self, x):
        x_enc = self.encoder(x)
        lastq = self.attn_lastq(x_enc)
        Z = F.normalize(lastq, 2, dim=2)
        
        return torch.matmul(Z, Z.transpose(1,2))


def det_memhap(SNVdataset: SNVMatrixDataset,
               ae: ReadAE, xformer: CorrTransformer,
               device: torch.cuda.device=torch.device("cpu"),
               num_clusters: int=2
               ):
    """
    SNVdataset:
        Dataset object containing read-SNP data
    ae:
        Autoencoder to generated read embeddings
    xformer:
        Transformer used to learn similarity matrix across reads
    num_clusters:
        Number of haplotypes
    device:  torch.cuda.device
        Device on which models are trained
    
    Returns
        haplotype labels
    """

    m = len(SNVdataset)
    dataloader_full = DataLoader(SNVdataset, batch_size=m,
                                num_workers=0)
    for i, (data, idx) in enumerate(dataloader_full):
        SNV_onehot = data.to(device)
        
    ae.eval()  # Set eval flags
    xformer.eval()  
    
    embed, _ = ae(SNV_onehot)  # Read embeddings
    # W_full = xformer(embed[None,:]).cpu().detach().numpy()  # Converting to numpy
    W_full = xformer(embed[None,:]).detach()
    W_full = W_full[0]  # Removing dummy dimension

    kmeans = KernelKMeans(n_clusters=num_clusters, max_iter=1000, tol=1e-3, device=device)
    kmeans.fit(W_full)
    
    if num_clusters == 2:  # Alphabet is {-1, +1}
        return 2*kmeans.labels_ - 1
    else:  # Alphabet is {0, 1, ,,,}
        return kmeans.labels_


# def origin2hap(SNV_matrix: torch.Tensor, origin: torch.Tensor,
#               num_hap: int=2, device: torch.cuda.device=torch.device("cpu")):    
#     """
#     SNV_matrix:
#         Full read-SNV matrix
#     origin: 
#         Specifies origin of each read by an int from (0, 1, ..., num_hap-1)
        
#     Returns
#         matrix of haplotypes (haplotypes x SNPs)
#     """
    
#     origin_val = torch.unique(origin)
#     accepted_val = torch.arange(num_hap, device=device)
#     if torch.any(np.intersect1d(origin_val, accepted_val) != origin_val):
#         raise ValueError("Invalid origin values passed as argument.")

#     hap_matrix = np.zeros((num_hap, SNV_matrix.shape[1]), dtype=int)
#     ACGTcount = ACGT_count(SNV_matrix)  # Stats of entire read matrix
#     for h in range(num_hap):
#         reads_h = SNV_matrix[origin == h]  # Reads attributed to haplotype i
#         h_stats = torch.zeros((SNV_matrix.shape[1], 4), device=device)
        
#         if len(reads_h) != 0:
#             h_stats = ACGT_count(reads_h) # ACGT statistics of a single nucleotide position
#         hap_matrix[h, :] = torch.argmax(h_stats, axis = 1) + 1  # Most commonly occuring base at each pos  
        
#         uncov_pos = torch.where(torch.sum(h_stats, axis = 1) == 0)[0]  # Positions uncovered by reads
#         for j in range(len(uncov_pos)):  # if not covered, select the most doninant one based on 'ACGTcount'  
#             base_max = np.flatnonzero(ACGTcount[uncov_pos[j], :] == np.amax(ACGTcount[uncov_pos[j], :])) + 1
#             if len(base_max) == 1:  # Single dominant base
#                 hap_matrix[h, uncov_pos[j]] == base_max[0]
#             else:  # Choose one of the dominant bases at random
#                 hap_matrix[h, uncov_pos[j]] = np.random.choice(base_max)

#     return hap_matrix

def origin2hap(SNV_matrix: NDArray[(Any, Any), int], origin: NDArray[int],
              num_hap: int=2) -> NDArray[(Any, Any), int]:    
    """
    SNV_matrix:
        Full read-SNV matrix
    origin: 
        Specifies origin of each read by an int from (0, 1, ..., num_hap-1)
        
    Returns
        matrix of haplotypes (haplotypes x SNPs)
    """
    
    origin_val = np.unique(origin)
    accepted_val = np.arange(num_hap)
    if np.any(np.intersect1d(origin_val, accepted_val) != origin_val):
        raise ValueError("Invalid origin values passed as argument.")

    hap_matrix = np.zeros((num_hap, SNV_matrix.shape[1]), dtype=int)
    ACGTcount = ACGT_count(SNV_matrix)  # Stats of entire read matrix
    for h in range(num_hap):
        reads_h = SNV_matrix[origin == h]  # Reads attributed to haplotype i
        h_stats = np.zeros((SNV_matrix.shape[1], 4))
        
        if len(reads_h) != 0:
            h_stats = ACGT_count(reads_h) # ACGT statistics of a single nucleotide position
        hap_matrix[h, :] = np.argmax(h_stats, axis = 1) + 1  # Most commonly occuring base at each pos  
        
        uncov_pos = np.where(np.sum(h_stats, axis = 1) == 0)[0]  # Positions uncovered by reads
        for j in range(len(uncov_pos)):  # if not covered, select the most doninant one based on 'ACGTcount'  
            base_max = np.flatnonzero(ACGTcount[uncov_pos[j], :] == np.amax(ACGTcount[uncov_pos[j], :])) + 1
            if len(base_max) == 1:  # Single dominant base
                hap_matrix[h, uncov_pos[j]] == base_max[0]
            else:  # Choose one of the dominant bases at random
                hap_matrix[h, uncov_pos[j]] = np.random.choice(base_max)

    return hap_matrix


def xformer_loss(xformer_output: torch.Tensor, origin: torch.Tensor,
                 Ws: torch.Tensor, Wmask: torch.Tensor,
                 lambda_reg: float = 0.1) -> float:
    
    origin_main = origin.float()  # origin.type(torch.FloatTensor)
    obj_main = -1*torch.matmul(origin_main, torch.matmul(xformer_output, origin_main))
    
    origin_onehot = F.one_hot(origin + 1).float()
    pair_mem = torch.matmul(origin_onehot, origin_onehot.transpose(0,1))
    obj_cont = torch.sum((1-pair_mem)*(1 + xformer_output)**2 + pair_mem*(1 - xformer_output)**2)
    
    obj_reg = lambda_reg*torch.linalg.norm(xformer_output*Wmask - Ws)**2 
    obj_sparse = torch.sum(torch.abs(xformer_output) - torch.diag_embed(torch.diag(xformer_output)))
    
    res = obj_cont + obj_reg + (1e1)*obj_sparse
#     res = obj_sparse
#     print('Clustering loss: %.3f, regularization: %.3f, sparsity: %.3f' %(obj_cont, obj_reg, obj_sparse))
    return res


def train_xhap(outhead: str, d_model: int = 128, num_hap: int = 2, num_epoch: int = 3000,
               check_cpr:bool = True, verbose: bool=False, gpu: int=-1
               )  -> Tuple[float, Optional[float]]:
    """
    Function to train XHap network using data indicated by outhead with the specified hyperparameters
    
    Parameters
    ----------
    outhead: str
        The prefix string for data (and consequently results)
    d_model: int
        Embedding size for each read
    n_hap: int
        Number of haplotypes to be assembled
    num_epoch: int
        Number of epochs to train XHap
    check_cpr: bool
        True if CPR is to be computed every epoch
    verbose: bool
        True if plots are to be printed
    gpu: int
        -1 for CPU, >=0 to spcify GPU to use
    
    Returns
    --------
    float:
        Lowest MEC score
    float:
        Corresponding CPR score (None if chec_cpr is False)
    """
    
    if gpu >= 0:
        if torch.cuda.is_available():
            device = torch.device("cuda:" + str(gpu))
        else:
            device = torch.device("cpu")
            print("GPUs not available, so CPU being used instead.")
    else:
        device = torch.device("cpu")
    print('DEVICE: ', device)

    datapath = 'generate_data/' + outhead + '/' + outhead + '_SNV_matrix.txt'  
    gt_file = 'generate_data/' + outhead + '/combined.fa'
    pos_file = 'generate_data/' + outhead + '/' + outhead + '_SNV_pos.txt'
    if check_cpr:
        true_haplo = read_true_hap(gt_file, pos_file)
    
    SNVdata = SNVMatrixDataset(datapath)
    SNV_matrix = np.loadtxt(datapath, dtype=int)
    SNV_matrix = SNV_matrix[np.sum(SNV_matrix != 0, axis=1) > 1]  # Removing uninformative reads
    print('SNP matrix: ', SNV_matrix.shape)

    
    nSNP = SNV_matrix.shape[1] # Number of SNVs
    num_read = len(SNVdata)  # Number of reads
    # batch_size = int(np.ceil(num_read/20))
    batch_size = int(np.ceil(num_read/5))
    d_model = 128  # Size of embeddings

    dataloader = DataLoader(SNVdata, batch_size=batch_size,
                            shuffle=True, num_workers=0)    
    
    # Instantiate read embedding encoder
    savefile="read_AE"
    embedAE = learn_embed(SNVdata, num_epoch=0, embed_dim=d_model, savefile=savefile)
    corr_xformer = CorrTransformer(d_model, d_model//2)  # Transformer
    
    # if use_gpu and torch.cuda.device_count() > 1:  # Using multiple GPUs
    #     embedAE = nn.DataParallel(embedAE)
    #     corr_xformer = nn.DataParallel(corr_xformer)
    embedAE.to(device)
    corr_xformer.to(device)

    xform_optimizer = optim.AdamW(list(corr_xformer.parameters()) + list(embedAE.parameters()),
                              lr=1e-5)
    
    # Set up logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    # Create handlers
    c_handler = logging.StreamHandler()
    f_handler = logging.FileHandler('xformer_train.log')
    c_handler.setLevel(logging.WARNING)
    f_handler.setLevel(logging.INFO)
    f_handler.addFilter(MyFilter(logging.INFO))

    # Create formatters and add it to handlers
    c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    f_format = logging.Formatter('%(message)s')
    c_handler.setFormatter(c_format)
    f_handler.setFormatter(f_format)

    # Add handlers to the logger
    logger.addHandler(c_handler)
    logger.addHandler(f_handler)
    
    # Setting up training parameters
    xform_train_loss_arr = []
    mec = []
    cpr = []
    xformer_savefile = 'generate_data/' + outhead + '/xhap_ckp'


    hap_origin = det_memhap(SNVdata, embedAE, corr_xformer, num_clusters=num_hap, device=device)  # Initial haplotype memberships
    if num_hap == 2:
        hap_matrix = origin2hap(SNV_matrix, (hap_origin.cpu().detach().numpy().astype(int) + 1)/2) 
    else:
        hap_matrix = origin2hap(SNV_matrix, hap_origin.cpu().detach().numpy().astype(int), num_hap) 
    mec.append(MEC(SNV_matrix, hap_matrix))
    if check_cpr:
        cpr.append(compute_cpr(hap_matrix, true_haplo))
    mec_min = np.inf
    cpr_max = 0

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
    W_tmp = (W_sim - W_dissim)/(W_sim + W_dissim + 1e-10)
    np.fill_diagonal(W_tmp, 1.)
    W_sup = torch.from_numpy(W_tmp)

    W_sup_dev = W_sup.to(device)
    W_mask_dev = W_mask.to(device)
    # hap_origin_dev = torch.from_numpy(hap_origin).to(device)

    time_xformer = []
    time_kmeans = []
    time_cpu = []

    for epoch in range(num_epoch):
        xform_train_loss = 0
        embedAE.train()  # Set train flags
        corr_xformer.train()

        t0 = time.time()
        for batch_data, batch_idx in dataloader:
            xform_optimizer.zero_grad()
            input_data = batch_data.to(device)

            embed, recon = embedAE(input_data)
            Z_batch = corr_xformer(embed[None,:])

            xform_loss = xformer_loss(Z_batch[0],
                                      hap_origin[batch_idx],
                                      W_sup_dev[batch_idx][:,batch_idx],
                                      W_mask_dev[batch_idx][:,batch_idx],
                                      lambda_reg=1e2) # + 0.1*AE_loss
            xform_loss.backward()
            xform_optimizer.step()
            xform_train_loss += xform_loss.item()

        xform_train_loss = xform_train_loss / len(dataloader)
        xform_train_loss_arr.append(xform_train_loss)
        t1 = time.time()

        hap_origin = det_memhap(SNVdata, embedAE, corr_xformer, num_clusters=num_hap,
                            device=device)  # Initial haplotype memberships
        t2 = time.time()
        
        if num_hap == 2:
            hap_matrix = origin2hap(SNV_matrix, (hap_origin.cpu().detach().numpy().astype(int) + 1)/2)
        else:
            hap_matrix = origin2hap(SNV_matrix, hap_origin.cpu().detach().numpy().astype(int), num_hap)
        t3 = time.time()

        time_xformer.append(t1-t0)
        time_kmeans.append(t2-t1)
        time_cpu.append(t3-t2)

        mec_curr = MEC(SNV_matrix, hap_matrix)
        if mec_curr <= mec_min:
            mec_min = mec_curr
            W_best = readW(SNVdata, embedAE, corr_xformer, device=device)
            hap_origin_best = 1*hap_origin
            hap_matrix_best = 1*hap_matrix
            epoch_best = 1*epoch
            print('Epoch = %d, MEC = %d' %(epoch, mec_curr))

            xhap_best = {
                'embed_ae': embedAE.state_dict(),
                'xformer': corr_xformer.state_dict()
                        }
            torch.save(xhap_best, 'generate_data/' + outhead + '/xhap_model')

        mec.append(mec_curr)

        if check_cpr:
            cpr_curr = compute_cpr(hap_matrix, true_haplo)
            if cpr_curr > cpr_max:
                cpr_max = cpr_curr
            cpr.append(cpr_curr)

            if (mec_curr == mec_min and 
                cpr_curr > compute_cpr(hap_matrix_best, true_haplo)):  # Store better solution
                W_best = readW(SNVdata, embedAE, corr_xformer, device=device)
                hap_origin_best = 1*hap_origin
                hap_matrix_best = 1*hap_matrix
                epoch_best = 1*epoch
                print('Epoch = %d, MEC = %d' %(epoch, mec_curr))

                xhap_best = {
                    'embed_ae': embedAE.state_dict(),
                    'xformer': corr_xformer.state_dict()
                            }
                torch.save(xhap_best, 'generate_data/' + outhead + '/xhap_model')

        # Display epoch training loss
        logger.info("epoch : {}/{}, loss = {:.2f}".format(epoch + 1, num_epoch, xform_train_loss))
        if epoch % 100 == 0:
            print("epoch : {}/{}, loss = {:.2f}".format(epoch + 1, num_epoch, xform_train_loss))
            # print("Time profiles: %.3f on transformer, %.3f on clustering, %.3f on origin2hap" 
            #         %(np.mean(time_xformer), np.mean(time_kmeans), np.mean(time_cpu)))
        if xformer_savefile and (epoch % 10 == 0):
            checkpoint = {
            'epoch': epoch + 1,
            'embed_ae': embedAE.state_dict(),
            'xformer': corr_xformer.state_dict(),
            'optimizer': xform_optimizer.state_dict()
            }
            save_ckp(checkpoint, xformer_savefile)

    if verbose:
    # Plot results
        if check_cpr:
            plt.figure(figsize=(10,18))
            plt.subplot(3,1,1)
            plt.plot(mec)
            plt.ylabel('MEC')
            plt.xlabel('Epoch')

            plt.subplot(3,1,2)
            plt.plot(cpr)
            plt.ylabel('CPR')
            plt.xlabel('Epoch')

            plt.subplot(3,1,3)
            plt.scatter(mec, cpr)
            plt.xlabel('MEC')
            plt.ylabel('CPR')

            idx_best = np.argsort(mec)[:10]
            plt.scatter(np.array(mec)[idx_best], np.array(cpr)[idx_best], color='red')

            print('Pearson correlation coefficient of MEC and CPR: %.3f' %(np.corrcoef(mec, cpr)[0, 1]))
            print(list(zip(np.array(mec)[idx_best], np.array(cpr)[idx_best])))

            # idx_same = np.where(np.diff(mec) == 0)
            # print(np.array(mec)[idx_same], np.array(cpr)[idx_same])
        else:  # CPR not computed
            plt.figure(figsize=(10,6))
            plt.plot(mec)
            plt.ylabel('MEC')
            plt.xlabel('Epoch')
    
    # Save parameters for best result
    hap_origin_best = hap_origin_best.cpu().numpy()
    if check_cpr:
        np.savez('generate_data/' + outhead + '/haptest_xformer_res', rec_hap=hap_matrix_best,
                 rec_hap_origin=hap_origin_best, true_hap=true_haplo)
    else:
        np.savez('generate_data/' + outhead + '/haptest_xformer_res', rec_hap=hap_matrix_best,
                 rec_hap_origin=hap_origin_best)


    if check_cpr:
        print('MAX CPR = %.3f, CORRESPONDING MEC = %d' %(cpr_max, mec[np.argmax(cpr)]))
        cpr_best = np.amax(np.array(cpr)[np.array(mec) == mec_min])
        print(np.array(cpr)[np.array(mec) == mec_min], np.argmin(mec))
        return mec[epoch_best + 1], cpr_best
    else:
        return mec[epoch_best + 1], None

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



