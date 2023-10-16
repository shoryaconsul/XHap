## Running experiments for haplotype assembly on large amounds of real data.
## read-SNP matrices are stored as sparse matrices, for efficient memory usage.
import numpy as np
import subprocess
import os
import shutil
from scipy.sparse import csr_matrix, lil_matrix, save_npz
from itertools import permutations
import argparse
from multiprocessing import Pool
import time

from helper import hamming_distance, read_sparseSNVMatrix, compute_cpr, MEC, read_hap
from read_embeddings import chunk_data
from xhap_sparse import train_xhap

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

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filehead", help="Prefix of required files",
        type=str, required=True)
    parser.add_argument("-p", "--ploidy", help="Ploidy of organism",
        default=2, type=int)
    parser.add_argument("-a", "--algo_runs", help="Number of experimental runs per dataset",
        default=1, type=int)
    parser.add_argument("-g", "--gpu", help='Number of GPUs to run XHap',
        default=1, type=int)
    parser.add_argument("--verbose", help="True for more verbose output",
        action='store_true', default=False)

    args = parser.parse_args()
    chunk_size = 100
    overlap_frac = 0.1

    # mec_expt = []
    # torch.multiprocessing.set_start_method('spawn')
    true_hap = read_hap('generate_data/NA12878/chr22/chr22_true_haplotypes.txt')

    datapath = 'generate_data/' + args.filehead + '_SNV_matrix.txt'
    savepath = 'generate_data/' + args.filehead + '_hap_matrix.npz'
    SNV_matrix = read_sparseSNVMatrix(datapath)
    # hap_matrix = np.zeros((args.ploidy, SNV_matrix.shape[1]), dtype=int)
    hap_matrix = np.load(savepath)['hap']

    def train_xhap_map(spos, SNVdata, gpu=-1, num_runs=1):
        mec_min = np.inf
        hap_matrix_best = np.zeros((args.ploidy, chunk_size), dtype=int)
        for r in range(num_runs):
            hap_matrix_run, mec_run = train_xhap(
				SNVdata, gpu=gpu, num_epoch=2000, verbose=args.verbose, num_hap=args.ploidy
				)
            if mec_run < mec_min:
                  mec_min = mec_run
                  hap_matrix_best = hap_matrix_run
        return spos, hap_matrix_best
    # train_xhap_map = partial(train_xhap, num_epoch=2000, verbose=args.verbose, num_hap=args.ploidy, num_runs=args.algo_runs)
    SNV_matrix_list = chunk_data(datapath, chunk_size=chunk_size, overlap_frac=overlap_frac)[-2:]
    gpu_list = range(7, 7 - args.gpu, -1)
    
    recon_end = 0  # End of reconstructed haplotype so far
    pool = Pool(processes=args.gpu)
    print('Created process pool')

    start_time = time.time()
    for d in range(0, len(SNV_matrix_list), args.gpu):
        chunk_starts, chunk_SNVdata = zip(*SNV_matrix_list[d:d + args.gpu])  # Chunk data for assembly
        print('Running XHap on chunks starting at: ', chunk_starts)
        res_d = pool.starmap(train_xhap_map, zip(chunk_starts, chunk_SNVdata, gpu_list))
        print('Finished running XHap on chunks starting at: ', chunk_starts)
        # Stitch haplotype chunks together
        pos_list, hap_chunk_list = zip(*res_d)
        for pos, hap_chunk in zip(pos_list, hap_chunk_list):
            np.savez('generate_data/NA12878/chr22/hap_pos_' + str(pos), hap_chunk)

            print('Chunk stats')
            print('MEC: ', MEC(SNV_matrix[:, pos:pos+hap_chunk.shape[1]].toarray(), hap_chunk))  
            print('CPR: ', compute_cpr(hap_chunk, true_hap[:, pos:pos+hap_chunk.shape[1]]))


            if pos == 0:
                hap_matrix[:, :hap_chunk.shape[1]] = hap_chunk
                recon_end = hap_chunk.shape[1]
            else:  # determine best match
                rec_hap = hap_matrix[:, pos:recon_end]
                match = best_match(rec_hap, hap_chunk[:, :recon_end - pos])
                hap_matrix[:, pos:pos + hap_chunk.shape[1]] = hap_chunk[match,:]
                recon_end = pos + hap_chunk.shape[1]
                # print('MATCH: ', match)

        print('Status so far')
        print('MEC: ', MEC(SNV_matrix[:, :recon_end].toarray(), hap_matrix[:, :recon_end]))
        print('CPR: ', compute_cpr(hap_matrix[:, :recon_end], true_hap[:, :recon_end]))

np.savez(savepath, hap = hap_matrix)   # Sve final result
print('Finished in %d seconds.' %(time.time()-start_time))
# print('MEC scores for XHap: ', mec)

# r_best = np.argmin(mec)
# print('Best MEC: %d' %(mec[r_best]))

# mec_expt.append(mec[r_best])

# print('MEC scores for real data: ', mec)
