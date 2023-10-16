## Running experiments for haplotype assembly
import numpy as np
import subprocess
import os
import sys
import shutil
import argparse
import time

current = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.dirname(current))
from helper import compute_cpr, MEC, compute_ver, read_hap, read_sparseSNVMatrix
from parse_hap import read_hap_blocks, compute_cpr_blocks, compute_MEC_blocks, \
                    compute_ver_blocks

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filehead", help="Prefix of required files",
	type=str, required=True)
parser.add_argument("-a", "--algo_runs", help="Number of experimental runs per dataset",
	default=1, type=int)
parser.add_argument("--mec_only", help="True when only MEC is to be computed",
    action='store_true', default=False)
parser.add_argument("--verbose", help="True for more verbose output",
    action='store_true', default=False)

hapcut_dir = 'HapCUT2'
extractHAIRS = 'HapCUT2/build/extractHAIRS'
hpop = "H-PoPG/H-PoPGv0.2.0.jar"  # Path to H-PoPG JAR file
haptree_dir = "./HapTree"  # Path to HapTree executable
ranbow = "ranbow/ranbow.py"
args = parser.parse_args()
ploidy =2  # diploid data
num_algos = 4  # XHap, CAECSeq, H-PoP, HapCUT2, HapTree (returns error)

mec_arr = np.zeros((num_algos))  # MEC
cpr_arr = np.zeros((num_algos))  # CPR
delta_mec_arr = np.zeros((num_algos))  # Delta MEC
ver_arr = np.zeros((num_algos))  # Vector error rate
mec_rate_arr = np.zeros((num_algos))  # MEC rate
ngaps_arr = np.zeros((num_algos))  # number of gaps
nblocks_arr = np.zeros((num_algos))  # number of blocks


fhead = args.filehead
mat_file = fhead + '_SNV_matrix.txt'
pos_file = fhead + '_SNV_pos.txt'
vcf_file = fhead + '_variants.vcf'
bam_file = fhead + '_reg_reads.bam'
frag_file = fhead + '_fragment_file.txt'
gt_file = fhead + '_true_haplotypes.txt'

SNV_matrix = read_sparseSNVMatrix(mat_file).toarray()

# Read best XHap result
xhap_res = np.load(fhead + '_hap_matrix.npz')
xhap_hap = xhap_res['hap']
true_hap = read_hap(gt_file)

mec_base = MEC(SNV_matrix, true_hap)
mec_xhap = MEC(SNV_matrix, xhap_hap)
mec_arr[0] = mec_xhap
delta_mec_arr[0] = mec_xhap - mec_base
cpr_arr[0] = compute_cpr(xhap_hap, true_hap)
ver_arr[0] = compute_ver(xhap_hap, true_hap)
mec_rate_arr[0] = mec_xhap/np.sum(SNV_matrix > 0)
ngaps_arr[0] = np.sum(xhap_hap == 0)
nblocks_arr[0] = np.sum(np.sum(xhap_hap != 0, axis=0) == 0) + 1


# Read best CAECSeq result
caec_hap = read_hap(fhead + '_Reconstructed_Strains.txt')
mec_caec = MEC(SNV_matrix, caec_hap)
mec_arr[1] = mec_caec
delta_mec_arr[1] = mec_caec - mec_base
cpr_arr[1] = compute_cpr(caec_hap, true_hap)
ver_arr[1] = compute_ver(caec_hap, true_hap)
mec_rate_arr[1] = mec_caec/np.sum(SNV_matrix > 0)
ngaps_arr[1] = np.sum(caec_hap == 0)
nblocks_arr[1] = np.sum(np.sum(caec_hap != 0, axis=0) == 0) + 1

os.chdir(hapcut_dir)
my_env = os.environ.copy()
if 'LD_LIBRARY_PATH' in my_env:
    my_env['LD_LIBRARY_PATH'] += os.pathsep + os.path.join(os.getcwd(), 'htslib')
else:
    my_env['LD_LIBRARY_PATH'] = os.path.join(os.getcwd(), 'htslib')

os.chdir('../')
subprocess.run([extractHAIRS, '--VCF', vcf_file, '--bam', bam_file, '--maxIS', '3000', '--out', frag_file], env=my_env)

# Run H-PoP
hpop_hapfile = fhead + '_hpop_haplotypes.txt'
for j in range(args.algo_runs):
    subprocess.run(['java', '-jar', hpop, '-p', str(ploidy), '-f', frag_file, '-o', hpop_hapfile])
    hpop_hap_blocks, hpop_gaps = read_hap_blocks(hpop_hapfile, ploidy)
    hpop_mec = compute_MEC_blocks(hpop_hap_blocks, vcf_file, SNV_matrix)
    if hpop_mec < mec_arr[2] or j == 0:
        mec_arr[2] = hpop_mec
        if not args.mec_only:
            cpr_arr[2] = compute_cpr_blocks(hpop_hap_blocks, vcf_file, true_hap)
            delta_mec_arr[2] = hpop_mec - mec_base
            ver_arr[2] = compute_ver_blocks(hpop_hap_blocks, vcf_file, true_hap)
            mec_rate_arr[2] = hpop_mec/np.sum(SNV_matrix > 0)
            ngaps_arr[2] = hpop_gaps
            nblocks_arr[2] = sum([np.sum(np.all(blk == -1, axis=0))
                                for blk in hpop_hap_blocks.values()]) + len(hpop_hap_blocks)
print('H-PoP done')

# Run HapCUT2
hapcut_hapfile = fhead + '_hapcut_haplotypes.txt'
for j in range(args.algo_runs):
    subprocess.run(['./' + hapcut_dir + '/build/HAPCUT2', '--fragments', frag_file,
                '--vcf', vcf_file, '--output', hapcut_hapfile])
    hapcut_hap_blocks, hapcut_gaps = read_hap_blocks(hapcut_hapfile, ploidy)
    hapcut_mec = compute_MEC_blocks(hapcut_hap_blocks, vcf_file, SNV_matrix)
    if hapcut_mec < mec_arr[3] or j == 0:
        mec_arr[3] = hapcut_mec
        if not args.mec_only:
            cpr_arr[3] = compute_cpr_blocks(hapcut_hap_blocks, vcf_file, true_hap)
            delta_mec_arr[3] = hapcut_mec - mec_base
            ver_arr[3] = compute_ver_blocks(hapcut_hap_blocks, vcf_file, true_hap)
            mec_rate_arr[3] = hapcut_mec/np.sum(SNV_matrix > 0)
            ngaps_arr[3] = hapcut_gaps
            nblocks_arr[3] = sum([np.sum(np.all(blk == -1, axis=0))
                            for blk in hapcut_hap_blocks.values()]) + len(hapcut_hap_blocks)
print('HapCUT2 done')

# # run HapTree
# haptree_hapdir = os.path.dirname(fhead) + '/' + 'haptree'
# haptree_hapfile = os.path.dirname(fhead) + '/' + 'haptree/HapTreeSolution'

# with open(vcf_file, 'r') as f:  # Count number of header lines in VCF file
#     vcf_head_size = 0
#     for line in f:
#         if line.startswith('#'):
#             vcf_head_size = vcf_head_size + 1
# f = open('tmp.vcf', 'w')
# subprocess.run(['tail', '-n', '+' + str(vcf_head_size + 1), vcf_file], stdout=f)  # Remove header lines
# f.close()
# for j in range(args.algo_runs):
#     print("STARTING HAPTREE")
#     t_0 = time.time()
#     subprocess.run([haptree_dir, frag_file, 'tmp.vcf', haptree_hapdir])
#     print('HapTree ran in %.2f seconds' %(time.time() - t_0))

#     haptree_hap_blocks, haptree_gaps = read_hap_blocks(haptree_hapfile, ploidy)
#     haptree_mec = compute_MEC_blocks(haptree_hap_blocks, vcf_file, SNV_matrix)
#     if haptree_mec < mec_arr[4] or j == 0:
#         mec_arr[4] = haptree_mec
#         if not args.mec_only:
#             cpr_arr[4] = compute_cpr_blocks(haptree_hap_blocks, vcf_file, true_hap)
#             delta_mec_arr[4] = haptree_mec - mec_base
#             ver_arr[4] = compute_ver_blocks(haptree_hap_blocks, vcf_file, true_hap)
#             mec_rate_arr[4] = haptree_mec/np.sum(SNV_matrix > 0)
#             ngaps_arr[4] = haptree_gaps
#             nblocks_arr[4] = sum([np.sum(np.all(blk == -1, axis=0))
#                             for blk in haptree_hap_blocks.values()]) + len(haptree_hap_blocks)
# os.remove('tmp.vcf')
# print('HapTree done')

print('Results for %s' %args.filehead)

# algo_names = ['XHap', 'CAECSeq', 'H-PoP', 'HapCUT2', 'HapTree']
algo_names = ['XHap', 'CAECSeq', 'H-PoP', 'HapCUT2']
num_algos = len(algo_names)
for row_idx in range(num_algos):
    col_idx = ~np.isinf(cpr_arr[row_idx])
    print(algo_names[row_idx] + '-----------------------')
    print('MEC: %d' %(np.mean(mec_arr[row_idx])))
    print('Delta MEC: %d' %(np.mean(delta_mec_arr[row_idx])))
    print('MEC rate: %.5f' %(np.mean(mec_rate_arr[row_idx])))
    print('CPR: %.3f' %(np.mean(cpr_arr[row_idx][col_idx])))
    print('VER: %.3f' %(np.mean(ver_arr[row_idx][col_idx])))
    print('Number of gaps: %.2f' %(np.mean(ngaps_arr[row_idx][col_idx])))
    print('Number of blocks: %.2f' %(np.mean(nblocks_arr[row_idx][col_idx])))
    print('All datasets yield valid results: ', np.all(col_idx))
