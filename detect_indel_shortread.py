import pysam
import re
import numpy as np
from scipy import sparse
from scipy.spatial.distance import cdist
from itertools import zip_longest
from collections import defaultdict

import subprocess
from os import remove

from matplotlib import pyplot as plt

from typing import Any, Union, List, Tuple, Dict, Optional
from nptyping import NDArray

from helper import hamming_distance
from helper_indel import *

def parse_samfile(fname: str, ref_genome: NDArray[int], rec_genomes: NDArray[(Any, Any), int],
                  is_paired: bool = True,
                 ) -> Tuple[Dict[Tuple[str, int], List],
                            Dict[Tuple[str, int], int],
                            List[List[Tuple[str, int]]]
                           ]:
    
    """
    This function parses the SAM file to create objects required by the other functions for
    indel deteciton from short reads.
    
    Parameters
    ----------
    fname: str
        Name of SAM file
    ref_genome: NDArray[int]
        Reference genome (bases converted to ints)
    rec_genomes: NDArray[(Any, Any), int]
        Recovered haplotypes (as genomes) 
    is_paired: bool
        Set to True if SAM file contains paired-end reads
    
    Returns
    -------
    Dict[Tuple[str, int], List]:
        Dictionary keyed by (read name, position) with value being list of
        CIGAR tuples, read sequence, read mapping quality and if read has
        supplementary alignments
    Dict[Tuple[str, int], int]:
        Dictionary keyed by (read name, position) with value set as the
        position of supplementary alignment
    List[List[Tuple[str, int]]]:
        List of lists of reads with each list of reads corresponding to
        one haplotype/genome
    
    """
    
    
    base_int = {'A': 1, 'C': 2, 'G': 3, 'T': 4, '-': 0}  # mapping of base to int

    samfile = pysam.AlignmentFile(fname, "r")
    sam_iter = samfile.fetch()

    read_sam_dict = {}  # Key is (read name, pos), value is list comprising CIGAR string, SEQ 
    read_idx_dict = {}  # Key is (read name, pos) and value is [row number of read, index of genome]
    read_list = []
    read_supp_dict = {}  # Key is (read name, pos) and value is position of supplementary alignment
    row_cnt = 0
    for i, x in enumerate(sam_iter):
        if i == 0:
            contig = str(x).split()[2]
        rname = x.query_name
        pos = x.reference_start  # Position
        read_sam_dict[(rname, pos)] = [x.cigartuples, x.query_sequence, x.mapping_quality,
                                       x.has_tag('SA')]  # REMOVE MAPPING QUALITY

        # if x.has_tag('SA'):
        #     print((rname, pos), ": ", x.cigarstring, "\t", x.get_tag('SA'), x.is_paired, x.is_read1, x.is_read2,
        #          x.next_reference_start, x.is_supplementary)

        q_ref_idx = x.get_aligned_pairs(matches_only=True)
        if q_ref_idx:  # If aligned pairs found
            q_idx, ref_idx = zip(*x.get_aligned_pairs(matches_only=True))
            query_bases = np.array([base_int[x.query_sequence[q_i]] for q_i in q_idx]) 

        # Check for supplementary/split alignments
        if x.has_tag('SA'):
            sa_str_split = x.get_tag('SA').split(';')[:-1]  # Dropping element due to trailing semicolon
            sa_pos_list = [int(z.split(',')[1]) -1 for z in sa_str_split]  # SAM file uses 1-based position
            read_supp_dict[(rname, pos)] = sa_pos_list[0]
            if is_paired:
                sa_pos_list.append(x.next_reference_start)  # Adding position of paired read for search
    #         print(sa_pos_list)
            split_read_found = False
            for sa_pos in sa_pos_list:
                if (rname, sa_pos) in read_idx_dict:  # Looking for primary alignment
    #                 print('Found primary alignment at %d' %sa_pos)
                    read_bases = read_list[read_idx_dict[(rname, sa_pos)]]
                    if q_ref_idx:
                        read_bases[np.array(ref_idx)] = query_bases
                    read_idx_dict[(rname, pos)] = read_idx_dict[(rname, sa_pos)]
                    split_read_found = True
                    break
            # No matching alignment
            if not split_read_found:
    #             print('Found no prior alignment')
                read_bases = np.zeros(len(ref_genome))
                if q_ref_idx:
                    read_bases[np.array(ref_idx)] = query_bases
                read_list.append(read_bases)
                read_idx_dict[(rname, pos)] = row_cnt
                row_cnt = row_cnt + 1

        elif (rname, x.next_reference_start) not in read_idx_dict:  # first paired read seen
            read_bases = np.zeros(len(ref_genome))
            if q_ref_idx:
                read_bases[np.array(ref_idx)] = query_bases
            read_list.append(read_bases)
            read_idx_dict[(rname, pos)] = row_cnt
            row_cnt = row_cnt + 1
        else:  # other paired read already seen
            read_bases = read_list[read_idx_dict[(rname, x.next_reference_start)]]
            if q_ref_idx:
                read_bases[np.array(ref_idx)] = query_bases
            read_idx_dict[(rname, pos)] = read_idx_dict[(rname, x.next_reference_start)]

    samfile.close()

    dist_genomes = cdist(np.array(read_list), rec_genomes, metric=hamming_distance)

    genome_read_list = [[] for _ in range(rec_genomes.shape[0])
                       ]  # List of lists, where each list comprises (rname, pos) associated with that genome 
    # Randomly select genome when distances to all genomes is equal
    perm = np.random.permutation(rec_genomes.shape[0])
    idx_genomes = perm[np.argmin(dist_genomes[:, perm], axis=1)]
    for (rname, pos), row_num in read_idx_dict.items():
        read_idx_dict[(rname, pos)] = [row_num, idx_genomes[row_num]]
        genome_read_list[idx_genomes[row_num]].append((rname, pos))

    return read_sam_dict, read_supp_dict, genome_read_list


def clipseq_from_cigar(pos: int, cigar: List[Tuple[int, int]], seq: str) -> Tuple[int, int, int,
                                                                                  bool, bool]:
    """
    This function processes the cigar tuple list for a given read and determines the sequence clipped
    from the read during alignment. It also returns the length of the clipping regions and whether
    the clipping occur at the start or end of the read.
    
    Parameters
    ----------
    pos: int
        Position of first non-clipped/inserted base on reference genome
    cigar: List[Tuple[int, int]]
        List of cigar tuples as returned in pysam, i.e., (operation, length)
    seq: str
        SEQ field in SAM/BAM file
    
    Returns
    -------
    pos_clip: int
        Reference position where clipped sequence starts/ends
    len_clip: int
        Length of clipped sequence
    len_match: int
        Length of matched sequence
    hard_clip: bool
        True if clipped sequence not in SEQ
    clip_at_start: bool
        True if clipped sequence occurs at start of read
    
    """
    
    c_arr, clen_arr = tuple(np.array(z) for z in zip(*cigar))
    clen_arr_clip = np.ma.array(clen_arr, mask=(c_arr == 0))  # Ignore matches
    idx_clip = np.argmax(clen_arr_clip)
    if idx_clip == 0:
        clip_atstart = True  # Clipped sequence at start of read
        len_clip = clen_arr[0]
    else:
        clip_atstart=False
        len_clip = clen_arr[idx_clip]
    
    if clip_atstart:  # Clipped sequence at start of read
        pos_clip = pos
        if c_arr[idx_clip] == 5:  # Hard clipping
            hard_clip = True
        else:  # Soft clipping
            hard_clip = False
    else:  # Clipped sequence at end of read
        pos_local = 0
        for idx, (c, clen) in enumerate(cigar[:idx_clip]):
            if c in [0, 2, 3, 7, 8]: # Cigar consumes reference bases
                pos_local = pos_local + clen
        pos_clip = pos + pos_local
        if c_arr[idx_clip] == 5:  # if Hard clipping
            hard_clip = True
        else:  # Soft clipping
            hard_clip = False
            
    len_match = int(np.sum(clen_arr) - len_clip)
    return pos_clip, len_clip, len_match, hard_clip, clip_atstart


def find_indel_readpair(pos1: int, cigar1: List[Tuple[int, int]], seq1: str,
                        pos2: int, cigar2: List[Tuple[int, int]], seq2: str
                       ) -> Tuple[bool, int, int, Optional[str]]:
    """
    This function determines whether the pair of chimeric reads indicate an insertion or deletion.
    Subsequently, it determines the start (reference) position of the indel and the length of the indel.
    In case of insertion, it also infers the sequence of the insertion.
    
    Parameters
    ----------
    pos1, pos2: int
        Reference positions of each read
    cigar1, cigar2:
        List of cigar tuples for the respective reads as returned in pysam, i.e., (operation, length)
    seq1, seq2: str
        Sequences contained in SEQ field of SAM/BAM file for each read
    
    Returns
    --------
    is_insert: bool
        True if insertion present, False if deletion occurs
    pos_insert/pos_del:
        Reference position of detected indel
    len_insert/len_del:
        Length of indel
    seq_insert: str
        Inserted sequence. None if this cannot be determined or deletion present
    
    
    """
    
    
    is_insert = False  # True if insertion detected
    seq_insert = None  # Inserted sequence 
    pos_indel1, len_clip1, len_match1, hard_clip1, clip_atstart1 = clipseq_from_cigar(pos1, cigar1, seq1)
    pos_indel2, len_clip2, len_match2, hard_clip2, clip_atstart2 = clipseq_from_cigar(pos2, cigar2, seq2)
#     print('From first read: ', pos_indel1, len_clip1, len_match1, clip_atstart1)
#     print('From second read: ', pos_indel2, len_clip2, len_match2, clip_atstart2)
    
    if np.isclose(pos_indel1, pos_indel2, atol=5):  # Indels detected at positions within 5bp of each other indicates an insert
        is_insert = True
        pos_insert = (pos_indel1 + pos_indel2) // 2
        len_insert = ((len_clip1 - len_match2) + (len_clip2 - len_match1)) // 2
    else:  # Deletion
        is_insert = False
        len_del = np.abs(pos_indel1 - pos_indel2)
        return is_insert, min(pos_indel1, pos_indel2), len_del, None
    
    # Accounting for offset in positions
    if clip_atstart1:
        len_clip1 = len_clip1 + pos_insert - pos_indel1
        len_match1 = len_match1 - (pos_insert - pos_indel1)
    else:
        len_clip1 = len_clip1 - (pos_insert - pos_indel1)
        len_match1 = len_match1 + pos_insert - pos_indel1
                                  
    if clip_atstart2:
        len_clip2 = len_clip2 + pos_insert - pos_indel2
        len_match2 = len_match2 - (pos_insert - pos_indel2)
    else:
        len_clip2 = len_clip2 - (pos_insert - pos_indel2)
        len_match2 = len_match2 + pos_insert - pos_indel2
    
    if not hard_clip1:  # Sequence not hard clipped
        if clip_atstart1:
            seq_insert = seq1[-1*len_match1 - len_insert: -1*len_match1]
        else:
            seq_insert = seq1[len_match1: len_match1 + len_insert]
    elif not hard_clip2:  # Sequence not hard clipped
        if clip_atstart2:
            seq_insert = seq2[-1*len_match2 - len_insert: -1*len_match2]
        else:
            seq_insert = seq2[len_match2: len_match2 + len_insert]
    else:  # hard clipping in both reads --> cannot find insertion sequence
        seq_insert = None

    return is_insert, pos_insert, len_insert, seq_insert


def merge_insert(ins_dict: Dict[int, List[Tuple[int, str]]], bin_size: int = 5
                ) -> List[Tuple[int, int, str]]:
    
    """
    This function aggregates nearby insertions among those found the (short) reads.In the case
    The inserted sequence for each position is recovered by the consensus of a multi-sequence
    alignemnt (MSA) if required.
    
    Parameters
    ----------
    ins_dict: Dict[int, List[Tuple[int, str]]]
        Dictionary with insert position as key and a list of insertion length
        and sequences as value
    bin_size:
        Inserts within bin_size//2 of each other are merged
    
    Returns
    -------
    List[Tuple[int, int, str]]:
        Each tuple is the insert positon, insert length and insert sequence
    
    """

    prev_pos = np.array(list(ins_dict.keys()))  # Ordered in Python 3.7 and later
    prev_dict = dict(ins_dict)
    curr_pos = []
    curr_dict = {}
    first_iter = True
    while list(prev_pos) != curr_pos:  # Keep merging indels until they are bin_size apart
        if not first_iter:
            prev_pos = np.array(curr_pos)
            prev_dict = dict(curr_dict)
        else:
            first_iter = False
        
        visited = np.zeros_like(prev_pos, dtype=bool) 
        curr_pos = []
        curr_dict = {}
        while not np.all(visited):  # Iterating until all bins visited
            idx = np.argmin(visited)
            p = prev_pos[idx]  # First unvisited position
            merge_idx = np.where(np.abs(prev_pos - p) <= bin_size//2)[0]
            merge_pos = int(np.rint(np.mean(prev_pos[merge_idx])))
            curr_pos.append(merge_pos)
            
            len_seq_tuples = []
            for idx in merge_idx:  # Collecting all lengths and sequences for merged bins
                len_seq_tuples = len_seq_tuples + prev_dict[prev_pos[idx]]
            
            curr_dict[merge_pos] = len_seq_tuples
            visited[merge_idx] = True
        
    # Find average length of insert and consensus insert sequence
    ins_len = []
    ins_seq = []
    for pos, len_seq_tuples in curr_dict.items():
        ins_len_list, ins_seq_list = tuple(map(list, zip(*len_seq_tuples)))
        ins_len.append(int(np.rint(np.mean(ins_len_list))))
        
        ########################
        tmp_ip_file = 'tmp_seq.fa'
        if len(ins_seq_list) > 1:  # Finding consensus sequence by MSA
            write_fasta(ins_seq_list, tmp_ip_file)

            # Run MAFFT for MSA
            proc = subprocess.run(['mafft', '--quiet', tmp_ip_file], capture_output=True, encoding="utf-8")
            aln_seq = read_msa(proc.stdout)
            remove(tmp_ip_file)  # Delete temporary file

            seq = []
            for seq_bases in zip_longest(*aln_seq, fillvalue='-'):
                base, base_cnt = np.unique(seq_bases, return_counts=True)
                base_mode = base[np.argmax(base_cnt)]
                if base_mode != '-': # Ignore gaps in consensu sequence
                    seq.append(base[np.argmax(base_cnt)])
            ins_seq.append(''.join(seq).upper())
        elif len(ins_seq_list) == 1:  # Only on inserted sequence
            ins_seq.extend(ins_seq_list)
    
    return list(zip(curr_pos, ins_len, ins_seq))


def merge_deletion(del_dict: Dict[int, List[int]], bin_size: int = 5
                  ) -> List[Tuple[int, int, None]]:

    """
    This function aggregates nearby deletions among those found the (short) reads.
    
    Parameters
    ----------
    ins_dict: Dict[int, List[int]]
        Dictionary with deletion position as key and a list of deletion length
        as value
    bin_size:
        Deletions within bin_size//2 of each other are merged
    
    Returns
    -------
    List[Tuple[int, int, None]]:
        Each tuple is the deletion positon and deletion length
    
    """

    prev_pos = np.array(list(del_dict.keys()))  # Ordered in Python 3.7 and later
    prev_dict = dict(del_dict)
    curr_pos = []
    curr_dict = {}
    first_iter = True
    while list(prev_pos) != curr_pos:  # Keep merging indels until they are bin_size apart
        if not first_iter:
            prev_pos = np.array(curr_pos)
            prev_dict = dict(curr_dict)
        else:
            first_iter = False
        
        visited = np.zeros_like(prev_pos, dtype=bool) 
        curr_pos = []
        curr_dict = {}
        while not np.all(visited):  # Iterating until all bins visited
            idx = np.argmin(visited)
            p = prev_pos[idx]  # First unvisited position
            merge_idx = np.where(np.abs(prev_pos - p) <= bin_size//2)[0]
            merge_pos = int(np.rint(np.mean(prev_pos[merge_idx])))
            curr_pos.append(merge_pos)
            
            len_list = []
            for idx in merge_idx:  # Collecting all lengths and sequences for merged bins
                len_list = len_list + prev_dict[prev_pos[idx]]
            
            curr_dict[merge_pos] = len_list
            visited[merge_idx] = True
        
    # Find average length of insert and consensus insert sequence
    del_len = []
    for pos, len_list in curr_dict.items():
        del_len.append(int(np.rint(np.mean(len_list))))
    
    return list(zip(curr_pos, del_len, [None]*len(curr_pos)))


def find_indel_shortread(read_sam_dict, read_supp_dict, genome_read_list
                        ) -> Tuple[List[List[Tuple[int, int, Optional[str]]]]]:
    
    """
    This function finds indels from short read data. It accomplishes this by finding supplementary
    alignments in the SAM file. These chimeric reads are inidcative of the presence of an indel.
    Based on their relative positons, we determine if it is an insertion or deletion.
    
    Parameters
    ----------
    read_sam_dict: Dict[Tuple[str, int], List]:
        Dictionary keyed by (read name, position) with value being list of
        CIGAR tuples, read sequence, read mapping quality and if read has
        supplementary alignments
    read_supp_dict: Dict[Tuple[str, int], int]:
        Dictionary keyed by (read name, position) with value set as the
        position of supplementary alignment
    genome_read_list: List[List[Tuple[str, int]]]:
        List of lists of reads with each list of reads corresponding to
        one haplotype/genome
        
    Returns
    -------
    List[List[Tuple[int, int, str]]]
        List of lists, where each list is the insertions found in that genome 
    
    List[List[Tuple[int, int, None]]]
        List of lists, where each list is the deletions found in that genome 
    
    Indels are returned as tuples of (start reference position, length of indel)
    
    
    """
    

    del_res = [[] for _ in range(len(genome_read_list))]
    ins_res = [[] for _ in range(len(genome_read_list))]
    for g_i, g_reads in enumerate(genome_read_list):    
        col_ptr = [0]
        col_idx = []  # col_idx[col_ptr[i]:col_ptr[i+1]] are the non-zero column indices for row i
        curr_ptr = 0

        read_seen = []
        insert_seq_list = []
        ins_dict = defaultdict(lambda: [])
        del_dict = defaultdict(lambda: [])

        for rname, pos in g_reads:
            if (rname, pos) in read_seen:
                continue
            else:
                is_chimeric = read_sam_dict[(rname, pos)][3]  # Read split as it overlaps with indel

            if is_chimeric:
                sa_pos = read_supp_dict[(rname, pos)]
                read_seen.extend([(rname, pos), (rname, sa_pos)])

                cigar1 = read_sam_dict[(rname, pos)][0]
                seq1 =  read_sam_dict[(rname, pos)][1]

                cigar2 = read_sam_dict[(rname, sa_pos)][0]
                seq2 = read_sam_dict[(rname, sa_pos)][1]

                is_insert, pos_indel, len_indel, seq_indel = find_indel_readpair(pos, cigar1, seq1, sa_pos, cigar2, seq2)

                if is_insert:
                    ins_dict[pos_indel].append((len_indel, seq_indel))
                else:
                    del_dict[pos_indel].append(len_indel)
        
        ins_res[g_i].extend(merge_insert(ins_dict))
        del_res[g_i].extend(merge_deletion(del_dict))
    
    return ins_res, del_res