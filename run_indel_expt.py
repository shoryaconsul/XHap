## Running experiments for indel detection

import numpy as np
import subprocess
import os
import argparse

from xhap_parallel import train_xhap
from detect_indel_shortread import parse_samfile, find_indel_shortread
from detect_indel_longread import parse_sam, find_indel_longread
from helper_indel import parse_genomes, global_align

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filehead", help="Prefix of required files",
	type=str, required=True)
parser.add_argument("-r", "--reference", help="Reference FASTA",
	type=str, required=True)
parser.add_argument("-c", "--cov", help="Coverage",
    default=20, type=float)
parser.add_argument("-i", "--indel", help="1 if insert, 2 if deleetion",
	default=0, choices=[0,1,2], type=int)
parser.add_argument("--idlen", help="Length of indel", type=int)
parser.add_argument("-n", "--num_expt", help="Number of experimental runs",
	default=1, type=int)
parser.add_argument("-g", "--gpu", help='GPU to run XHap',
    default=-1, type=int)
parser.add_argument("--long", help="True if using long reads",
    action='store_true', default=False)
args = parser.parse_args()

short_reads = not args.long  # = True for short reads, False for long reads
if args.indel == 1:
	true_ins_pos = []
	rec_ins_pos = []
	rec_ins_score = []
elif args.indel == 2:
	true_del_pos = []
	rec_del_pos = []
	rec_del_len = []

cov = np.round(args.cov/2, 1)
for r in range(args.num_expt):
	# Generate data
	os.chdir('generate_data')
	if short_reads:  # Short reads
		if args.indel == 0:  # No indel
			fhead = args.filehead + "_norm_iter" + str(r+1)
			subprocess.run(['bash', 'gendata_semiexp.bash', '-f', args.reference,
				'-o', fhead, '-n 2', '-c 20'])
		elif args.indel == 1:  # Insertion
			fhead = args.filehead + "_ins" + str(args.idlen) + "_iter" + str(r+1)
			subprocess.run(['bash', 'gendata_semiexp.bash', '-f', args.reference,
				'-o', fhead, '-n 2', '-c', str(cov), '-i', str(args.idlen)])
		else:  # Deletion
			fhead = args.filehead + "_del" + str(args.idlen) + "_iter" + str(r+1)
			subprocess.run(['bash', 'gendata_semiexp.bash', '-f', args.reference,
				'-o', fhead, '-n 2', '-c', str(cov), '-d', str(args.idlen)])
	else:  # Long reads
		if args.indel == 0:  # No indel
			fhead = args.filehead + "_norm_long_iter" + str(r+1)
			subprocess.run(['bash', 'gendata_longread_semiexp.bash', '-f', args.reference,
				'-o', fhead, '-c', str(cov), '-n 2'])
		elif args.indel == 1:  # Insertion
			fhead = args.filehead + "_ins" + str(args.idlen) + "_iter" + str(r+1)
			subprocess.run(['bash', 'gendata_longread_semiexp.bash', '-f', args.reference,
				'-o', fhead, '-n 2', '-c', str(cov), '-i', str(args.idlen)])
		else:  # Deletion
			fhead = args.filehead + "_del" + str(args.idlen) + "_iter" + str(r+1)
			subprocess.run(['bash', 'gendata_longread_semiexp.bash', '-f', args.reference,
				'-o', fhead, '-n 2', '-c', str(cov), '-d', str(args.idlen)])

	os.chdir('../')

	# Train XHap on generated data
	train_xhap(fhead, check_cpr=False, num_epoch=2000, gpu=args.gpu, num_hap=2)

	# Detect indels
	if args.indel > 0:
		ref_genome, rec_genomes = parse_genomes(fhead, args.reference)
		if short_reads:
			samfile = "generate_data/" + fhead + "/" + fhead +".sam"
			read_sam_dict, read_supp_dict, geneome_read_list = parse_samfile(
				samfile, ref_genome, rec_genomes, is_paired=True)
			ins_res, del_res = find_indel_shortread(read_sam_dict,
				read_supp_dict, geneome_read_list)
		else:
			read_idx_dict = parse_sam(fhead, ref_genome, rec_genomes)
			ins_res, del_res = find_indel_longread(fhead, read_idx_dict, rec_genomes)


	indel_file = "generate_data/" + fhead + "/indel.txt"
	with open(indel_file, "r") as f:
		line = f.readline().strip()
		indel_pos = int(line.split(':')[1])
		print('INDEL POS', indel_pos)
		if args.indel == 1:  # If insertion, look for insertion sequence
			line = f.readline().strip()
			ins_seq = line.split(':')[1].strip()

	if args.indel == 1:  # Evaluating detected insertions
		best_ins = None
		best_score = np.inf
		for ins_list in ins_res:  # Find best matching insertion
			if ins_list:
				score_list = [np.abs(p-indel_pos) + np.abs(args.idlen - ilen)
				for p, ilen, iseq in ins_list]
				best_idx = np.argmin(score_list)
				if score_list[best_idx] < best_score:
					best_ins = ins_list[best_idx]
					best_score = score_list[best_idx]

		if best_ins is not None:
			true_ins_pos.append(indel_pos)
			rec_ins_pos.append(best_ins[0])

			aln_seq1, aln_seq2, aln_score = global_align(ins_seq, best_ins[2])
			ins_score = sum([int(a == b) 
				for a, b in zip(aln_seq1, aln_seq2)])/args.idlen
			rec_ins_score.append(ins_score)

	elif args.indel == 2:  # Best matching deletion
		best_del = None
		best_score = np.inf
		for del_list in del_res:  # Find best matching insertion
			if del_list:
				score_list = [np.abs(p-indel_pos) + np.abs(args.idlen - dlen)
				for p, dlen, _ in del_list]
				best_idx = np.argmin(score_list)
				if score_list[best_idx] < best_score:
					best_del = del_list[best_idx]
					best_score = score_list[best_idx]

		if best_del is not None:
			true_del_pos.append(indel_pos)
			rec_del_pos.append(best_del[0])
			rec_del_len.append(best_del[1])

if args.indel == 1:
	dev_pos = np.array(true_ins_pos) - np.array(rec_ins_pos)
	print('INSERTIONS: ', rec_ins_pos, rec_ins_score)
	print('FNR: %.3f, Avg pos deviation: %.2f +/- %.2f' %(
		1 - len(rec_ins_pos)/args.num_expt,
		np.mean(dev_pos), np.std(dev_pos))
	)
	print('Avg alignment score: %.3f +/- %.3f' %(
		np.mean(rec_ins_score), np.std(rec_ins_score)))
elif args.indel == 2:
	dev_pos = np.array(true_del_pos) - np.array(rec_del_pos)
	print('DELETIONS: ', rec_del_pos, rec_del_len)
	print('FNR: %.3f, Avg pos deviation: %.2f +/- %.2f' %(
		1 - len(rec_del_pos)/args.num_expt,
		np.mean(dev_pos), np.std(dev_pos))
	)
	print('Avg deviation in length: %.2f +/- %.2f' %(
		np.mean(rec_del_len) - args.idlen, np.std(rec_del_len))
	)
