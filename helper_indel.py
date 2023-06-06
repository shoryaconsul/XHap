import pysam
import re
import numpy as np
from scipy import sparse
from scipy.spatial.distance import cdist
from itertools import zip_longest, chain, product
import subprocess
from os import remove

from matplotlib import pyplot as plt

from typing import Any, Union, List, Tuple, Dict, Optional
from nptyping import NDArray

from helper import hamming_distance


def write_fasta(seq_list: List[str], fname: str) -> None:
    """
    This function takes in a list of sequences and writes them to a FASTA file
    with unique dummy headers.
    
    Parameters
    ----------
    seq_list: List[str]
        List of sequences to be stored in FASTA file
    fname: str
        Name of FASTA file to be written into
        Note that the contents of this file will be overwritten.
        
    This function has no return value.
    
    """
    
    with open(fname, 'w') as f:
        for i, seq in enumerate(seq_list):
            f.write('>header%d\n'%i)
            f.write(seq)
            f.write('\n')


def read_msa(msa_out: str) -> List[str]:
    """
    This function reads the output of an MSA algorithm. Each alignment in such a 
    file is assumed to be preceded by a header of the form ">[header]" and a single
    alignment may be broken into multiple lines. The output is a list of these
    alignments without the headers.
    
    Parameters
    ----------
    msa_out: str
        This is the MSA algorithm read as a string (as would happen if read from stdout).
        msa_out can be replaced with msa_list if the desired input form is a list of strings.
    
    Returns
    -------
    List[str]
        List of the alignments without the headers and any of the preceding and trailing
        whitespaces.
    
    """
    
    msa_list = msa_out.split()
    aln_seq = []

    for i, line in enumerate(msa_list):
        if line[0] == '>':  # New header
            if i > 0:
                aln_seq.append(seq)
            seq = ''
        else:
            seq = seq + line.strip()
    aln_seq.append(seq)

    return aln_seq


def parse_genomes(filehead: str, ref: str) -> Tuple[NDArray[int], NDArray[(Any, Any), int]]:
    
    """
    Parse the reference genome and assemble the recovered genomes from the reconstructed haplotypes.
    
    Parameters
    ----------
    filehead: str
        SAM file is found as 'generate_data/filehead/filehead.sam'
    ref: str
        Reference FASTA is 'generate_data/ref_file'

    Returns
    -------
    NDArray[int]:
        Reference genome
    NDArray[(Any, Any), int]:
        Reconstructed genomes

    """
    
    ref_file = "generate_data/" + ref
    
    xhap_res = np.load("generate_data/" + filehead + "/haptest_xformer_res.npz")
    rec_hap = xhap_res['rec_hap']

    pos_file = "generate_data/" + filehead + "/" + filehead + "_SNV_pos.txt"

    base_int = {'A': 1, 'C': 2, 'G': 3, 'T': 4, '-': 0}  # mapping of base to int
    int_base = {1: 'A', 2: 'C', 3: 'G', 4: 'T', 0: '-'}  # mapping of base to int

    with open(ref_file, "r") as f:  # Read in reference sequence
        while True:
            line = f.readline().strip()
            if not line:
                break
            elif line[0] != '>':
                ref_seq = line

    ref_genome = np.zeros(len(ref_seq), dtype=np.int8)
    for i, b in enumerate(ref_seq):
        ref_genome[i] = base_int[b]

    with open(pos_file, "r") as f:  # Read in SNP positions
        pos_str = f.readline().split()
        pos = np.array([int(ps) - 1 for ps in pos_str]) # Convert to int, assuming pos_file is 1-indexed

    rec_genomes = np.tile(ref_genome, (np.shape(rec_hap)[0], 1))
    for i, r in enumerate(rec_hap):
        for j, p in enumerate(pos):
            rec_genomes[i][p] = r[j]
    
    return ref_genome, rec_genomes


def global_align(seq1: str, seq2: str) -> Tuple[str, str, int]:
    
    def match_score(a, b):
        return 2*(a == b) - 1  # 1 if match, -1 if mismatch
    
    gap_score = -1  # Score for gap
    
    M = len(seq1)
    N = len(seq2)
    
    score = np.zeros((M+1, N+1))
    dirn = np.zeros((M+1, N+1))  # 0 means diagonal, 1 means right, 2 means down
    
    # Initialize first row and column of score matrix
    score[0, :] = -1*np.arange(N+1)
    score[:, 0] = -1*np.arange(M+1)
    dirn[0, :] = 1
    dirn[:, 0] = 2
    
    
    # DP to compute alignment scores
    for i, j in product(range(1, M+1), range(1, N+1)):
        score_ij = [score[i-1, j-1] + match_score(seq1[i-1], seq2[j-1]),
                    score[i, j-1] + gap_score,
                    score[i-1, j] + gap_score]
        dirn[i, j] = np.argmax(score_ij)
        score[i, j] = np.amax(score_ij)
    
    # Backtracking to produce alignment
    aln1, aln2 = '', ''
    i, j = M, N
    while i > 0 or j > 0:
        dirn_ij = dirn[i, j]
        if dirn_ij == 0:  # Match or mismatch, i.e. diagonal
            aln1 = seq1[i-1] + aln1
            aln2 = seq2[j-1] + aln2
            i = i - 1
            j = j - 1
        elif dirn_ij == 1:  # Gap in seq2, i.e., right
            aln1 = seq1[i-1] + aln1
            aln2 = '-' + aln2
            j = j - 1
        else:  # Gap in seq2, i.e., down
            aln1 = '-' + aln1
            aln2 = seq2[j-1] + aln2
            i = i -1
    
    return aln1, aln2, score[-1, -1]