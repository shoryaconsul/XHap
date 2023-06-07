## Running experiments for haplotype assembly
import numpy as np
import subprocess
import os
import shutil
import argparse

from xhap_parallel import train_xhap
from detect_indel_shortread import parse_samfile, find_indel_shortread
from detect_indel_longread import parse_sam, find_indel_longread
from helper_indel import parse_genomes, global_align

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filehead", help="Prefix of required files",
	type=str, required=True)
parser.add_argument("-p", "--ploidy", help="Ploidy of organism",
    default=2, type=int)
parser.add_argument("-n", "--num_expt", help="Number of datasets",
	default=1, type=int)
parser.add_argument("-a", "--algo_runs", help="Number of experimental runs per dataset",
	default=1, type=int)
parser.add_argument("-g", "--gpu", help='GPU to run XHap',
    default=-1, type=int)
parser.add_argument("--verbose", help="True for more verbose output",
    action='store_true', default=False)

args = parser.parse_args()
mec_expt = []
for i in range(args.num_expt):

    fhead = args.filehead + str(i+1)
    mec = []
    for r in range(args.algo_runs):
        # Train XHap on generated data
        print('RUN %d for %s' %(r+1, fhead))
        mec_r, cpr_r = train_xhap(fhead, check_cpr=False, num_epoch=2000,
                                gpu=args.gpu, verbose=args.verbose, num_hap=args.ploidy)
        if len(mec) == 0 or mec_r < min(mec):
            shutil.copy('generate_data/' + fhead + '/haptest_xformer_res.npz',
                        'generate_data/' + fhead + '/haptest_xformer_res_best.npz')
        mec.append(mec_r)

    print('MEC scores for XHap: ', mec)

    r_best = np.argmin(mec)
    print('Best MEC: %d' %(mec[r_best]))

    mec_expt.append(mec[r_best])

print('MEC scores for real data: ', mec)