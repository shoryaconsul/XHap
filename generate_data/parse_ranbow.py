import numpy as np
import os
import sys

from typing import Any, Dict, Tuple, List
from nptyping import NDArray

current = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.dirname(current))
from helper import hamming_distance  #, PermDict, gen_perm, permute_distance

def convert_ranbow_out(out_file: str, hap_mhm_file: str):
    """
    This function reads the *.hap file and returns the haplotypes in a MHM file.

    Parameters
    ----------
    out_file: str
        Ranbow result file (extension: *.hap)
    hap_mhm_file: str
        Output MHM file
    """

    fi = open(out_file,'r')
    fo = open(hap_mhm_file,'w')
    lines = fi.readlines()

    haplen = int(lines[0].split()[2])
    fo.write(lines[0])
    inx = 1
    while inx < len(lines):
        if lines[inx][0]==">":
            fo.write(lines[inx])
        else:
            a = lines[inx].split()
            fo.write(a[5].split("_")[0]+"\t"+ a[1]+"\t"+ a[2]+"\t"+ a[3]+'\n')
        inx +=1


def read_a_methods_result(methods_haps_file: str) -> Dict[str, Dict[str, List[Tuple[str]]]]:
    """
    This function reads Ranbow results from the input MHM file.

    Parameters
    ----------
    methods_haps_file: str
        MHM file

    Returns
    -------
    Dict[str, Dict[str, List[Tuple[str]]]]
        [scaffold: {Block ID, [(Start SNP index, haplotype string)}]
    """

    fi = open(methods_haps_file)
    lines = fi.readlines()

    i = 0
    block = 0
    map_method = {}
    while i < len(lines):
        if lines[i][0]==">":
            _, scaf , _ = lines[i].split()
            map_method [scaf] = {}
            i += 1
        block_index, _, start, hap = lines[i].split()
        base_block_index = block_index
        map_method[scaf][block_index] = []
        while base_block_index == block_index:
            map_method[scaf][block_index] += [(start, hap)]
            i += 1
            if i >= len(lines):
                break
            block_index, _, start, hap = lines[i].split()
    return map_method


def ranbow_hap(map_method: Dict[str, Dict[str, List[Tuple[str]]]], vcf_file: str
              ) -> List[Tuple[int, NDArray[int]]]:
    """
    Convert the binary haplotype strings into haplotype arrays of integers

    Parameters
    ----------
    map_method: Dict[str, Dict[str, List[Tuple[str]]]]
        Dictionary parsed from MHM file
    vcf_file: str
        Variants VCF file
    
    Returns
    -------
    List[Tuple[int, NDArray[int]]]:
        List of tuples of (haplotype start index, haplotype array)
    """
    base_int = {'A': 1, 'C': 2, 'G': 3, 'T': 4, '-': 0}  # mapping of base to int
    variant_alleles = []

    # Read VCF file
    with open(vcf_file, "r") as f:
        line = f.readline()
        while True:
            line = f.readline()
            if not line:
                break
            else:
                line_split = line.split()
                alleles = [line_split[3]] + line_split[4].split(',')
                variant_alleles.append([base_int[b] for b in alleles])    
#     print(variant_alleles)

    blocks = [*map_method.values()][0]
    res = []
    for block_idx, hap_list in blocks.items():
        tmp = []
        for start_idx, bhap in hap_list:
            si = int(start_idx)
            hap = np.array([variant_alleles[i+si][int(b)] if b != '-' else -1 for i, b in enumerate(bhap)])
            tmp.append((si, hap))
        res = res + tmp
    return res


def ranbow_MEC(SNV_matrix: NDArray[(Any, Any), int],
               hap_list: List[Tuple[int, NDArray[int]]]) -> int:
    """
    This function computes the MEC for haplotypes reconstructed by Ranbow.

    Parameters
    ----------
    SNV_matrix: NDArray[(Any, Any), int]
        read-SNP matrix
    hap_list: List[Tuple[int, NDArray[int]]]
        List of haplotype tuples - each tuple contains starting SNP index and haplotype array
    
    Returns
    --------
    int: MEC score
    """

    def read_overlap(SNV_read: NDArray[int], si: int, haplen: int) -> bool:
        """
        True when read overlaps haplotype of length haplen starting at si.
        """
        nz_idx = np.nonzero(SNV_read)[0]
        # print(si, haplen, nz_idx)
        if len(nz_idx) == 0:
            return False
        elif si <= nz_idx[0]:
            return si + haplen >= nz_idx[0]
        elif si > nz_idx[0]:
            return si <= nz_idx[-1]

    res = 0
    for SNV_read in SNV_matrix:
        dis = [
            hamming_distance(SNV_read[si:si+len(hap)], hap)
                if read_overlap(SNV_read, si, len(hap)) else np.inf
            for si, hap in hap_list
            ]
        # if np.isinf(min(dis)):
        #     print(np.nonzero(SNV_read)[0])
        res = res + min(dis)
        # print(dis, res)
        
    return res