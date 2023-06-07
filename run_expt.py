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
parser.add_argument("-r", "--reference", help="Reference FASTA",
	type=str, required=True)
parser.add_argument("-p", "--ploidy", help="Ploidy of organism",
    default=2, type=int)
parser.add_argument("-c", "--cov", help="Coverage",
    default=10, type=float)
parser.add_argument("-n", "--num_expt", help="Number of datasets",
	default=1, type=int)
parser.add_argument("-a", "--algo_runs", help="Number of experimental runs per dataset",
	default=1, type=int)
parser.add_argument("-g", "--gpu", help='GPU to run XHap',
    default=-1, type=int)
parser.add_argument("--long", help="True if using long reads",
    action='store_true', default=False)
parser.add_argument("--verbose", help="True for more verbose output",
    action='store_true', default=False)

args = parser.parse_args()
mec_expt = []
cpr_expt = []
for i in range(args.num_expt):

    fhead = args.filehead + "_iter" + str(i+1)
    Generate data
    cov = np.round(args.cov/args.ploidy, 3)
    os.chdir('generate_data')
    if not args.long:  # Short reads
        if args.verbose:
            subprocess.run(['bash', 'gendata_semiexp.bash', '-f', args.reference,
                '-o', fhead, '-n', str(args.ploidy), '-c', str(cov), '-v'])
        else:
            subprocess.run(['bash', 'gendata_semiexp.bash', '-f', args.reference,
                '-o', fhead, '-n', str(args.ploidy), '-c', str(cov)])
    else:  # Long reads
        if args.verbose:
            subprocess.run(['bash', 'gendata_longread_semiexp.bash', '-f', args.reference,
                '-o', fhead, '-n', str(args.ploidy), '-c', str(cov), '-v'])
        else:
            subprocess.run(['bash', 'gendata_longread_semiexp.bash', '-f', args.reference,
                '-o', fhead, '-n', str(args.ploidy), '-c', str(cov)])

    os.chdir('../')

    mec = []
    cpr = []
    for r in range(args.algo_runs):
        # Train XHap on generated data
        print('RUN %d for %s' %(r+1, fhead))
        mec_r, cpr_r = train_xhap(fhead, check_cpr=True, num_epoch=2000,
                                gpu=args.gpu, verbose=args.verbose, num_hap=args.ploidy)
        if len(mec) == 0 or mec_r < min(mec):
            shutil.copy('generate_data/' + fhead + '/haptest_xformer_res.npz',
                        'generate_data/' + fhead + '/haptest_xformer_res_best.npz')
            shutil.copy('generate_data/' + fhead + '/xhap_model',
                        'generate_data/' + fhead + '/xhap_model_best')
        mec.append(mec_r)
        cpr.append(cpr_r)

    print('MEC scores for XHap: ', mec)
    print('CPR scores for XHap: ', cpr)

    r_best = np.argmin(mec)
    print('Best MEC: %.3f, Corresponding CPR: %.3f' %(mec[r_best], cpr[r_best]))

    mec_expt.append(mec[r_best])
    cpr_expt.append(cpr[r_best])

print('Average MEC: %.3f +/- %.3f, Average CPR: %.3f +/- %.3f' 
    %(np.mean(mec_expt), np.std(mec_expt),
    np.mean(cpr_expt), np.std(cpr_expt)))