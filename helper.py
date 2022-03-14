import numpy as np
from os import path
from itertools import permutations

import logging
from typing import Any
from nptyping import NDArray


def save_ckp(state, checkpoint_path):
    """
    state: checkpoint we want to save
    checkpoint_path: path to save checkpoint
    """
    f_path = checkpoint_path  # Save path
    torch.save(state, f_path)


def load_ckp(checkpoint_path, model, optimizer):
    """
    checkpoint_path: path to save checkpoint
    model: model that we want to load checkpoint parameters into       
    optimizer: optimizer we defined in previous training
    """
    checkpoint = torch.load(checkpoint_path)
    # initialize state_dict from checkpoint to model
    model.load_state_dict(checkpoint['state_dict'])
    # initialize optimizer from checkpoint to optimizer
    optimizer.load_state_dict(checkpoint['optimizer'])

    # return model, optimizer, epoch value
    return model, optimizer, checkpoint['epoch']


class MyFilter(object):
	"""
		Filter for logging module
	"""
	def __init__(self, level):
		self.__level = level
	def filter(self, logRecord):
		return logRecord.levelno <= self.__level


def create_logger(logname: str, logfile: str) -> logging.Logger:

	logger = logging.getLogger(logname)
	logger.setLevel(logging.INFO)

	# Create handlers
	c_handler = logging.StreamHandler()
	f_handler = logging.FileHandler(logfile)
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
	
	return logger


# get the ACGT statistics of a read matrix
def ACGT_count(submatrix: NDArray[(Any, Any), int]):
	"""
	submatrix:
		Read-SNV marix (m x n)

	Returns
		(n x 4) matrix of base counts at each SNP
	"""
	out = np.zeros((submatrix.shape[1], 4))
	for i in range(4):
		out[:, i] = (submatrix == (i + 1)).sum(axis = 0)

	return out

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


def hamming_distance(read: NDArray[(Any,), int], 
	haplo: NDArray[(Any,), int]) -> int:
	"""
	read:
		Read denoted by base at each SNP (1-D numpy array)
	haplo:
		Haplotype denoted by base at each SNP (1-D numpy array)

	Returns
		Hamming distance between read and haplotype 

	"""
	if np.shape(read) != np.shape(haplo):
		raise ValueError('Read and haplotype must be of the same dimension.')

	return sum((haplo - read)[read != 0] != 0)


def MEC(SNV_matrix: NDArray[(Any, Any), int],
        hap_matrix: NDArray[(Any, Any), int]) -> int:  # Compute MEC score
    
    """
	SNV_matrix:
		Read-SNP matrix
	hap_matrix:
		Haplotype-SNP matrix

	Returns
		MEC score for given SNV matrix and haplotypes

    """

    if np.shape(SNV_matrix)[1] != np.shape(hap_matrix)[1]:
    	raise ValueError("Different number of SNPs in reads and haplotypes.")

    res = 0
    
    for SNV_read in SNV_matrix:
        dis = [hamming_distance(SNV_read, hap) for j, hap in enumerate(hap_matrix)]
        res = res + min(dis)
        
    return res


# evaluate the correct phasing rate
def compute_cpr(recovered_haplo: NDArray[(Any, ), int], 
	true_haplo: NDArray[(Any,), int]) -> float:
	"""
	recovered_haplo:
		k x n matrix of recovered haplotypes
	true_haplo:
		True haplotypes (ground truth)

	Returns
		correct phasing rate
	"""

	if np.shape(recovered_haplo) != np.shape(true_haplo):
		raise ValueError("Input arguments should have the same shape.")

	distance_table = np.zeros((len(recovered_haplo), len(true_haplo)))
	for i, rec_hap in enumerate(recovered_haplo):
		for j, true_hap in enumerate(true_haplo):
			distance_table[i, j] = hamming_distance(rec_hap, true_hap)

	index = permutations(range(true_haplo.shape[0]))
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
	# index = (list(index))[np.argmin(np.array(distance))]  # Best one-to-one mapping
	# print(best_matching)
	cpr = 1 - min(distance) / np.size(true_haplo)

	return cpr


def read_true_hap(gt_file: str, pos_file: str) -> NDArray[(Any, Any), int]:
	"""
	gt_file:
		File containing true haplotypes (k haplotypes)
	pos_file:
		File containing SNP positions (n positions)

	Returns
		k x n numpy array of haplotypes
	"""
	if not path.isfile(gt_file):
		raise OSError("Ground truth file not found.")

	if not path.isfile(pos_file):
		raise OSError("SNP position file not found.")

	with open(pos_file, "r") as f:
		pos_str = f.readline().split()
		pos = np.array([int(ps) for ps in pos_str]) # Convert to int, assuming pos is 0-indexed

	true_hap = []
	base_int = {'A': 1, 'C': 2, 'G': 3, 'T': 4}  # mapping of base to int
	with open(gt_file, "r") as f:
		while True:
			line = f.readline()
			if not line:
				break
			if line[0] != '>':  # ignore header strings
				hap = np.array([base_int[line[p]] for p in pos])
				true_hap.append(hap.astype('int'))

	return np.array(true_hap)

		
