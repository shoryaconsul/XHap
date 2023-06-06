# File to read reconstructed haplotypes and compute MEC/CPR

import numpy as np
from os import path
import sys

current = path.dirname(path.realpath(__file__))
sys.path.append(path.dirname(current))
import re
from helper import compute_cpr, MEC, compute_ver, read_true_hap
import argparse

from typing import Any, Dict, Tuple
from nptyping import NDArray

def read_hap_blocks(hap_file: str, num_hap: int
                ) -> Tuple[Dict[int, NDArray], int]:
    """
    Read reconstructed haplotypes and save them as a dictionary of blocks
    
    Parameters
    ----------
    hap_file:
        File containing reconstructed haplotypes
    hap__num:
        Number of reconstructed haplotypes
        
    Returns
        Dictionary of haplotype blocks keyed by index of starting SNP for each haplotype block
        Number of missing haplotype entries
    """
    
    if not path.isfile(hap_file):
            raise OSError("Haplotype file not found.")
    
    blocks_dict = {}  # Key is first SNP number of block
    num_gaps = 0  # Number of missing haplotype values
    with open(hap_file, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            elif line[0] in ["#", "B"]:
                snp_idx = int(re.split(r"[ :]+", line)[2])
                block_len = int(re.split(r"[ :]+", line)[4])
                hap_block = np.zeros((num_hap, block_len), dtype=np.int32)
            elif line[0] != "*":
                line_split = line.split()
                if int(line_split[0]) - snp_idx >= block_len:
                    hap_block = np.append(hap_block, np.zeros((num_hap, int(line_split[0]) - snp_idx - block_len + 1), dtype=np.int32), axis=1)
                    block_len = int(line_split[0]) - snp_idx + 1
                hap_block[:, int(line_split[0]) - snp_idx] = np.array([-1 if h == '-' else int(h) for h in line_split[1:1+num_hap]])
                num_gaps = num_gaps + np.sum(hap_block[:, int(line_split[0]) - snp_idx] == -1)
            else:
                blocks_dict[snp_idx] = hap_block
    
    # Reached end of file
    if ~np.all(hap_block == 0):
        blocks_dict[snp_idx] = hap_block
    
    return blocks_dict, num_gaps


def compute_cpr_blocks(blocks_dict: Dict[int, NDArray], vcf_file: str, true_haplo: NDArray[(Any, Any), int]):
    """
    This function computes the CPR of each haplotype block and returns the average CPR
    across the blocks.

    Parameters
    ----------
    blocks_dict: Dict[int, NDArray]
        Dictionary of haplotype blocks keyed by starting SNP index
    vcf_file: str
        Variants file
    true_haplo: NDArray[(Any, Any), int]
        True haplotype matrix
    
    Returns
    -------
    float: Average CPR across haplotype blocks 

    """
    if not path.isfile(vcf_file):
            raise OSError("Variant VCF file not found.")
            
    base2int = {'A': 1, 'C': 2, 'G': 3, 'T': 4, '-': 0}  # mapping of base to int
    
    variant_alleles = []
    with open(vcf_file, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            elif line[0] != "#":
                line_split = line.split()
                alleles = [line_split[3]] + line_split[4].split(',')
                variant_alleles.append([base2int[b] for b in alleles])
    
    cpr = []
    for snp_idx, hap_block in blocks_dict.items():
        block_len = np.shape(hap_block)[1]
        hap_block_bases = np.zeros(np.shape(hap_block))
        for (i, j), h in np.ndenumerate(hap_block):  # Determine base at each SNP for inferred haplotypes
            if h >= 0:
                hap_block_bases[i, j] = variant_alleles[snp_idx - 1 + j][h]
            else:
                hap_block_bases[i, j] = -1  # So that missing bases are penalized in CPR
        
        cpr.append(compute_cpr(hap_block_bases, true_haplo[:, snp_idx-1:snp_idx-1+block_len])*block_len)
    
    return sum(cpr)/np.shape(true_haplo)[1]


def compute_MEC_blocks(blocks_dict: Dict[int, NDArray], vcf_file: str, SNV_matrix: NDArray[(Any, Any), int]):
    """
    This function computes the MEC of each haplotype block and returns the total MEC
    across the blocks.

    Parameters
    ----------
    blocks_dict: Dict[int, NDArray]
        Dictionary of haplotype blocks keyed by starting SNP index
    vcf_file: str
        Variants file
    tSNV_matrix: NDArray[(Any, Any), int]
        read=SNP matrix
    
    Returns
    -------
    int: Total MEC over haplotype blocks

    """   
    if not path.isfile(vcf_file):
            raise OSError("Variant VCF file not found.")
            
    base2int = {'A': 1, 'C': 2, 'G': 3, 'T': 4, '-': 0}  # mapping of base to int
    
    variant_alleles = []
    with open(vcf_file, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            elif line[0] != "#":
                line_split = line.split()
                alleles = [line_split[3]] + line_split[4].split(',')
                variant_alleles.append([base2int[b] for b in alleles])
    
    mec = []
    for snp_idx, hap_block in blocks_dict.items():
        block_len = np.shape(hap_block)[1]
        hap_block_bases = np.zeros(np.shape(hap_block))
        for (i, j), h in np.ndenumerate(hap_block):  # Determine base at each SNP for inferred haplotypes
            if h >= 0:
                hap_block_bases[i, j] = variant_alleles[snp_idx - 1 + j][h]
            else:
                hap_block_bases[i, j] = 0
        
        mec.append(MEC(SNV_matrix[:, snp_idx-1:snp_idx-1+block_len], hap_block_bases))
    
    return sum(mec)

def compute_ver_blocks(blocks_dict: Dict[int, NDArray], vcf_file: str, true_haplo: NDArray[(Any, Any), int]):
    """
    This function computes the switch/vector errors in each haplotype block and returns
    the vector/switch error rate over the full haplotypes. This ignores the vector/switch
    errors between haplotype blocks. Switch error rate (SWER) is returned for diploids,
    otherwise the vector error rate (VER) is returned.

    Parameters
    ----------
    blocks_dict: Dict[int, NDArray]
        Dictionary of haplotype blocks keyed by starting SNP index
    vcf_file: str
        Variants file
    true_haplo: NDArray[(Any, Any), int]
        True haplotype matrix
    
    Returns
    -------
    float: SWER/VER

    """
    if not path.isfile(vcf_file):
            raise OSError("Variant VCF file not found.")
            
    base2int = {'A': 1, 'C': 2, 'G': 3, 'T': 4, '-': 0}  # mapping of base to int
    
    variant_alleles = []
    with open(vcf_file, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            elif line[0] != "#":
                line_split = line.split()
                alleles = [line_split[3]] + line_split[4].split(',')
                variant_alleles.append([base2int[b] for b in alleles])
    
    ver = []
    for snp_idx, hap_block in blocks_dict.items():
        block_len = np.shape(hap_block)[1]
        hap_block_bases = np.zeros(np.shape(hap_block))
        for (i, j), h in np.ndenumerate(hap_block):  # Determine base at each SNP for inferred haplotypes
            if h >= 0:
                hap_block_bases[i, j] = variant_alleles[snp_idx - 1 + j][h]
            else:
                hap_block_bases[i, j] = -1  # So that missing bases are penalized in CPR
        
        ver.append(compute_ver(hap_block_bases, true_haplo[:, snp_idx-1:snp_idx-1+block_len])*(block_len-1))
    
    return sum(ver)/(np.shape(true_haplo)[1] - 1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Evaluate quality of reconstructed haplotyes.')
    parser.add_argument('-p', '--pos', required=True, help='File containing 1-indexed SNP positions') 
    parser.add_argument('-g', '--gt', help='FASTA file containing true genome sequences')
    parser.add_argument('-v', '--vcf', required=True, help='Variants VCF file')
    parser.add_argument('-f', '--file', required=True, help='File containing reconstructed haplotypes')
    parser.add_argument('-m', '--mat', required=True, help='Integer-valued read-SNP matrix file')
    # parser.add_argument('-o', '--out', required=True, help='Output file')

    args = parser.parse_args()

    true_hap = read_true_hap(args.gt, args.pos)
    SNV_matrix = np.loadtxt(args.mat, dtype=int)
    mec_base = MEC(SNV_matrix, true_hap)  # Minimum possible MEC
    
    num_hap = np.shape(true_hap)[0]  # Number of haplotypes
    hap_blocks, hap_gaps = read_hap_blocks(args.file, num_hap)   
    #hap_blocks = read_hap_blocks(args.file, 4)   
    recon_hap = {list(hap_blocks.keys())[0]: np.hstack(tuple(hap_blocks.values()))}
    
    cpr_blocks = compute_cpr_blocks(recon_hap, args.vcf, true_hap)
    mec_blocks = compute_MEC_blocks(recon_hap, args.vcf, SNV_matrix)
    ver_blocks = compute_ver_blocks(recon_hap, args.vcf, true_hap)
    print("CPR: ", cpr_blocks,
         " MEC: ", mec_blocks,
         "Diff MEC: ", mec_blocks - mec_base,
         "VER: ", ver_blocks)

    # return cpr_blocks, mec_blocks, len(hap_blocks), hap_gaps
    
    #print(" MEC: ", compute_MEC_blocks(recon_hap, args.vcf, SNV_matrix))
