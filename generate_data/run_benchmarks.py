## Running experiments for haplotype assembly
import numpy as np
import subprocess
import os
import sys
import shutil
import argparse

current = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.dirname(current))
from helper import compute_cpr, MEC, compute_ver, read_true_hap, read_hap
from snp2vcf import gen_vcf
from parse_hap import read_hap_blocks, compute_cpr_blocks, compute_MEC_blocks, \
                    compute_ver_blocks
from parse_ranbow import *

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filehead", help="Prefix of required files",
	type=str, required=True)
parser.add_argument("-r", "--reference", help="Reference FASTA",
	type=str, required=True)
parser.add_argument("-p", "--ploidy", help="Ploidy of organism",
    default=2, type=int)
parser.add_argument("-n", "--num_expt", help="Number of datasets",
	default=1, type=int)
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
if args.ploidy > 2:
    num_algos = 4  # XHap, CAECSeq, H-PoP, Ranbow
elif args.ploidy == 2:
    num_algos = 5  # XHap, CAECSeq, H-PoP, HapCUT2, HapTree

mec_arr = np.zeros((num_algos, args.num_expt))  # MEC
cpr_arr = np.zeros((num_algos, args.num_expt))  # CPR
delta_mec_arr = np.zeros((num_algos, args.num_expt))  # Delta MEC
ver_arr = np.zeros((num_algos, args.num_expt))  # Vector error rate
mec_rate_arr = np.zeros((num_algos, args.num_expt))  # MEC rate
ngaps_arr = np.zeros((num_algos, args.num_expt))  # number of gaps
nblocks_arr = np.zeros((num_algos, args.num_expt))  # number of blocks

for i in range(args.num_expt):
    if not args.mec_only:
        fhead = args.filehead + "_iter" + str(i+1)
    else:
        fhead = args.filehead
    mat_file = fhead + '/' + fhead + '_SNV_matrix.txt'
    pos_file = fhead + '/' + fhead + '_SNV_pos.txt'
    vcf_file = fhead + '/' + fhead + '_variants.vcf'
    bam_file = fhead + '/' + fhead + '_sorted.bam'
    frag_file = fhead + '/' + fhead + '_fragment_file.txt'

    SNV_matrix = np.loadtxt(mat_file, dtype=int)

    # Read best XHap result
    xhap_res = np.load(fhead + '/haptest_xformer_res_best.npz')
    xhap_hap = xhap_res['rec_hap']
    mec_xhap = MEC(SNV_matrix, xhap_hap)
    mec_arr[0, i] = mec_xhap

    if not args.mec_only:
        true_hap = xhap_res['true_hap']

        mec_base = MEC(SNV_matrix, true_hap)
        delta_mec_arr[0, i] = mec_xhap - mec_base
        cpr_arr[0, i] = compute_cpr(xhap_hap, true_hap)
        ver_arr[0, i] = compute_ver(xhap_hap, true_hap)
        mec_rate_arr[0, i] = mec_xhap/np.sum(SNV_matrix > 0)
        ngaps_arr[0, i] = np.sum(xhap_hap == 0)
        nblocks_arr[0, i] = np.sum(np.sum(xhap_hap != 0, axis=0) == 0) + 1

    # Read best CAECSeq result
    caec_hap = read_hap(fhead + '/' + fhead + '_Reconstructed_Strains.txt')
    mec_caec = MEC(SNV_matrix, caec_hap)
    mec_arr[1, i] = mec_caec
    if not args.mec_only:
        delta_mec_arr[1, i] = mec_caec - mec_base
        cpr_arr[1, i] = compute_cpr(caec_hap, true_hap)
        ver_arr[1, i] = compute_ver(caec_hap, true_hap)
        mec_rate_arr[1, i] = mec_caec/np.sum(SNV_matrix > 0)
        ngaps_arr[1, i] = np.sum(caec_hap == 0)
        nblocks_arr[1, i] = np.sum(np.sum(caec_hap != 0, axis=0) == 0) + 1

    # Create VCF file for other benchmarks
    gen_vcf(args.reference, pos_file, mat_file, vcf_file)
    
    os.chdir(hapcut_dir)
    my_env = os.environ.copy()
    if 'LD_LIBRARY_PATH' in my_env:
        my_env['LD_LIBRARY_PATH'] += os.pathsep + os.path.join(os.getcwd(), 'htslib')
    else:
        my_env['LD_LIBRARY_PATH'] = os.path.join(os.getcwd(), 'htslib')

    os.chdir('../')

    subprocess.run([extractHAIRS, '--VCF', vcf_file, '--bam', bam_file, '--maxIS', '3000', '--out', frag_file], env=my_env)

    # Run H-PoP
    hpop_hapfile = fhead + '/' + fhead + '_hpop_haplotypes.txt'
    for j in range(args.algo_runs):
        subprocess.run(['java', '-jar', hpop, '-p', str(args.ploidy), '-f', frag_file, '-o', hpop_hapfile])
        hpop_hap_blocks, hpop_gaps = read_hap_blocks(hpop_hapfile, args.ploidy)
        hpop_mec = compute_MEC_blocks(hpop_hap_blocks, vcf_file, SNV_matrix)
        if hpop_mec < mec_arr[2,i] or j == 0:
            mec_arr[2,i] = hpop_mec
            if not args.mec_only:
                cpr_arr[2,i] = compute_cpr_blocks(hpop_hap_blocks, vcf_file, true_hap)
                delta_mec_arr[2,i] = hpop_mec - mec_base
                ver_arr[2,i] = compute_ver_blocks(hpop_hap_blocks, vcf_file, true_hap)
                mec_rate_arr[2,i] = hpop_mec/np.sum(SNV_matrix > 0)
                ngaps_arr[2, i] = hpop_gaps
                nblocks_arr[2, i] = sum([np.sum(np.all(blk == -1, axis=0))
                                    for blk in hpop_hap_blocks.values()]) + len(hpop_hap_blocks)
    
    if args.ploidy == 2:
        # Run HapCUT2
        hapcut_hapfile = fhead + '/' + fhead + '_hapcut_haplotypes.txt'
        for j in range(args.algo_runs):
            subprocess.run(['./' + hapcut_dir + '/build/HAPCUT2', '--fragments', frag_file,
                        '--vcf', vcf_file, '--output', hapcut_hapfile])
            hapcut_hap_blocks, hapcut_gaps = read_hap_blocks(hapcut_hapfile, args.ploidy)
            hapcut_mec = compute_MEC_blocks(hapcut_hap_blocks, vcf_file, SNV_matrix)
            if hapcut_mec < mec_arr[3,i] or j == 0:
                mec_arr[3,i] = hapcut_mec
                if not args.mec_only:
                    cpr_arr[3,i] = compute_cpr_blocks(hapcut_hap_blocks, vcf_file, true_hap)
                    delta_mec_arr[3,i] = hapcut_mec - mec_base
                    ver_arr[3,i] = compute_ver_blocks(hapcut_hap_blocks, vcf_file, true_hap)
                    mec_rate_arr[3,i] = hapcut_mec/np.sum(SNV_matrix > 0)
                    ngaps_arr[3, i] = hapcut_gaps
                    nblocks_arr[3, i] = sum([np.sum(np.all(blk == -1, axis=0))
                                    for blk in hapcut_hap_blocks.values()]) + len(hapcut_hap_blocks)
    
        # run HapTree
        haptree_hapdir = fhead + '/' + 'haptree'
        haptree_hapfile = fhead + '/' + 'haptree/HapTreeSolution'
        f = open('tmp.vcf', 'w')
        subprocess.run(['tail', '-n', '+2', vcf_file], stdout=f)
        f.close()
        for j in range(args.algo_runs):
            subprocess.run([haptree_dir, frag_file, 'tmp.vcf', haptree_hapdir])
            haptree_hap_blocks, haptree_gaps = read_hap_blocks(haptree_hapfile, args.ploidy)
            haptree_mec = compute_MEC_blocks(haptree_hap_blocks, vcf_file, SNV_matrix)
            if haptree_mec < mec_arr[4,i] or j == 0:
                mec_arr[4,i] = haptree_mec
                if not args.mec_only:
                    cpr_arr[4,i] = compute_cpr_blocks(haptree_hap_blocks, vcf_file, true_hap)
                    delta_mec_arr[4,i] = haptree_mec - mec_base
                    ver_arr[4,i] = compute_ver_blocks(haptree_hap_blocks, vcf_file, true_hap)
                    mec_rate_arr[4,i] = haptree_mec/np.sum(SNV_matrix > 0)
                    ngaps_arr[4, i] = haptree_gaps
                    nblocks_arr[4, i] = sum([np.sum(np.all(blk == -1, axis=0))
                                    for blk in haptree_hap_blocks.values()]) + len(haptree_hap_blocks)
        os.remove('tmp.vcf')

    if args.ploidy > 2:
        # Run Ranbow
        ranbow_par_file = fhead + '/' + 'ranbow_params'
        if 'soltub' in fhead:
            scaf_file = 'soltub_scaffold.list'
        elif 'human' in fhead:
            scaf_file = 'human_scaffold.list'

        # # Uncomment for real data
        # scaf_file = 'soltub_scaffold_1.list'
        # with open(args.reference, 'r') as f:
        #     while True:
        #         line = f.readline()
        #         if not line:
        #             break
        #         elif line[0] == '>':
        #             fheader = line[1:].strip()
                
        # with open(scaf_file, 'w') as f:
        #     f.write(fheader.split()[0])

        for j in range(args.algo_runs):
            with open(ranbow_par_file, 'w') as f:
                f.write("-ploidy %d\n" %args.ploidy)
                f.write("-noProcessor 1\n")
                f.write("-bamFile %s\n" %bam_file)
                f.write("-refFile %s\n" %args.reference)
                f.write("-vcfFile %s\n" %vcf_file)
                f.write("-selectedScf %s\n" %scaf_file)
                f.write("-outputFolderBase %s\n" %(fhead + '/ranbow'))
            
            subprocess.run(['python2', ranbow, 'hap', '-mode', 'index', '-par', ranbow_par_file,
                    '-processorIndex', '0'])
            subprocess.run(['python2', ranbow, 'hap', '-mode', 'hap', '-par', ranbow_par_file,
                    'processorIndex', '0', '-WinLen', '8'])
            subprocess.run(['python2', ranbow, 'hap', '-mode', 'collect', '-par', ranbow_par_file])
            
            # Read Ranbow results
            ranbow_hapfile = fhead + '/ranbow/ranbow.single.hap'
            ranbow_mhmfile = fhead + '/ranbow/ranbow_mhmfile.txt'

            convert_ranbow_out(ranbow_hapfile, ranbow_mhmfile)
            ranbow_hap_dict = read_a_methods_result(ranbow_mhmfile)
            ranbow_hap_list = ranbow_hap(ranbow_hap_dict, vcf_file)

            ranbow_mec = ranbow_MEC(SNV_matrix, ranbow_hap_list)
            if ranbow_mec < mec_arr[3, i] or j == 0:
                mec_arr[3,i] = ranbow_mec
                if not args.mec_only:
                    delta_mec_arr[3,i] = ranbow_mec - mec_base
                    mec_rate_arr[3,i] = ranbow_mec/np.sum(SNV_matrix > 0)
                    hap_list_len = np.array(list(map(lambda x: len(x[1]), ranbow_hap_list)))
                    if np.all(hap_list_len == hap_list_len[0]):
                        ranbow_hap_matrix = np.vstack(list(zip(*ranbow_hap_list))[1])
                        if np.shape(true_hap) == np.shape(ranbow_hap_matrix):
                            ranbow_cpr = compute_cpr(ranbow_hap_matrix, true_hap)
                            ranbow_ver = compute_ver(ranbow_hap_matrix, true_hap)
                            cpr_arr[3,i] = ranbow_cpr
                            ver_arr[3,i] = ranbow_ver
                            ngaps_arr[3, i] = np.sum(ranbow_hap_matrix == 0) + np.sum(ranbow_hap_matrix == -1)
                            nblocks_arr[3, i] = np.sum(np.sum(ranbow_hap_matrix != -1, axis=0) == 0) + 1
                    else:
                        cpr_arr[3,i] = np.inf
                        ver_arr[3,i] = np.inf
                        ngaps_arr[3, i] = np.inf
                        nblocks_arr[3, i] = np.inf
# print(cpr_arr, mec_arr, delta_mec_arr, ver_arr)
print(ver_arr)
print('Results for %s' %args.filehead)

if args.ploidy == 2:
    algo_names = ['XHap', 'CAECSeq', 'H-PoP', 'HapCUT2', 'HapTree']
else:
    algo_names = ['XHap', 'CAECSeq', 'H-PoP', 'Ranbow']
for row_idx in range(num_algos):
    col_idx = ~np.isinf(cpr_arr[row_idx])
    print(algo_names[row_idx] + '-----------------------')
    print('MEC: %d +/- %.2f' %(np.mean(mec_arr[row_idx]), np.std(mec_arr[row_idx])))
    print('Delta MEC: %d +/- %.2f' %(np.mean(delta_mec_arr[row_idx]), np.std(delta_mec_arr[row_idx])))
    print('MEC rate: %.5f +/- %.5f' %(np.mean(mec_rate_arr[row_idx]), np.std(mec_rate_arr[row_idx])))
    print('CPR: %.3f +/- %.3f' %(np.mean(cpr_arr[row_idx][col_idx]), np.std(cpr_arr[row_idx][col_idx])))
    print('VER: %.3f +/- %.3f' %(np.mean(ver_arr[row_idx][col_idx]), np.std(ver_arr[row_idx][col_idx])))
    print('Number of gaps: %.2f +/- %.2f' %(np.mean(ngaps_arr[row_idx][col_idx]), np.std(ngaps_arr[row_idx][col_idx])))
    print('Number of blocks: %.2f +/- %.2f' %(np.mean(nblocks_arr[row_idx][col_idx]), np.std(nblocks_arr[row_idx][col_idx])))
    print('All datasets yield valid results: ', np.all(col_idx))
