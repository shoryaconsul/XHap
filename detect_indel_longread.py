import pysam
import re
import numpy as np
from scipy import sparse
from scipy.spatial.distance import cdist
from itertools import zip_longest, chain, product
import subprocess
from os import remove

from typing import Any, Union, List, Tuple, Dict, Optional
from nptyping import NDArray

from helper import hamming_distance
from helper_indel import *


def parse_sam(filehead: str, ref_genome: NDArray[int], rec_genomes: NDArray[(Any, Any), int]
             ) -> Dict[Tuple[str, int], List[int]]:
    
    """
    This function parses the SAM file and generates a read-index dictionary specifying which
    row of the read-SNP matrix the read is mapped to as well as the genome it is 'closest' to.
    
    Parameters
    ----------
    filehead: str
        SAM file is found as 'generate_data/filehead/filehead.sam'
    ref_genome: NDArray[int]
        Reference genome (bases converted to ints)
    rec_genomes: NDArray[(Any, Any), int]
        Recovered haplotypes (as genomes) 

    Returns
    --------
    Dict[Tuple[str, int], List[int]]:
        Read-idx dictionary
    
    """
    
    base_int = {'A': 1, 'C': 2, 'G': 3, 'T': 4, '-': 0}  # mapping of base to int
    
    samfile = pysam.AlignmentFile("generate_data/" + filehead + "/" + filehead +".sam", "r")
    sam_iter = samfile.fetch()

    read_idx_dict = {}  # Key is (read name, pos) and value is [row number of read, index of genome]
    read_list = []
    read_supp_dict = {}  # Key is (read name, pos) and value is position of supplementary alignment
    row_cnt = 0
    for i, x in enumerate(sam_iter):
        if i == 0:
            contig = str(x).split()[2]
        if not x.is_unmapped:
            rname = x.query_name
            pos = x.reference_start  # Position
            q_ref_idx = x.get_aligned_pairs(matches_only=True)
            if q_ref_idx:  # If aligned pairs found
                q_idx, ref_idx = zip(*x.get_aligned_pairs(matches_only=True))
                query_bases = np.array([base_int[x.query_sequence[q_i]] for q_i in q_idx]) 
            # q_idx, ref_idx, query_bases = zip(*map(lambda y, z: (*y, base_int[z]),
            #                                        x.get_aligned_pairs(matches_only=True), x.query_sequence))

            # Check for supplementary/split alignments
            if x.has_tag('SA'):
                sa_str_split = x.get_tag('SA').split(';')[:-1]  # Dropping element due to trailing semicolon
                sa_pos_list = [int(z.split(',')[1]) -1 for z in sa_str_split]  # SAM file uses 1-based position
                read_supp_dict[(rname, pos)] = sa_pos_list[0]
    #             sa_pos_list.append(x.next_reference_start)  # Adding position of paired read for search
                split_read_found = False
                for sa_pos in sa_pos_list:
                    if (rname, sa_pos) in read_idx_dict:  # Looking for primary alignment
                        read_bases = read_list[read_idx_dict[(rname, sa_pos)]]
                        if q_ref_idx:
                            read_bases[np.array(ref_idx)] = query_bases
                        read_idx_dict[(rname, pos)] = read_idx_dict[(rname, sa_pos)]
                        split_read_found = True
                        break
                # No matching alignment
                if not split_read_found:
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

    # Randomly select genome when distances to all genomes is equal
    dist_genomes = cdist(np.array(read_list), rec_genomes, metric=hamming_distance)
    perm = np.random.permutation(rec_genomes.shape[0])
    idx_genomes = perm[np.argmin(dist_genomes[:, perm], axis=1)]
    for (rname, pos), row_num in read_idx_dict.items():
        read_idx_dict[(rname, pos)] = [row_num, idx_genomes[row_num]]
    
    return read_idx_dict


def merge_bins(pos: NDArray[int], cnt_arr: NDArray[int], len_arr: NDArray[int], 
               seq_list: Optional[List[List[str]]], bin_size: int
              ) -> Union[Tuple[NDArray[int], NDArray[int], NDArray[int]], Optional[List[str]]]:
    
    """
    This function iteratively merges bins until bin centers are separated by at least bin_size positions.
    
    Parameters
    ----------
    pos: NDArray[int]:
        Bin centers for bins with non-zero counts
    cnt_arr: NDArray[int]
        Counts in each bin
    len_arr: NDArray[int]
        Length of each quantity (indel) in each bin
    seq_list: Optional[List[str]]
        If not None, list of lists where each list comprises inserted
        sequences starting at corresponding reference position 
    bin_size: int
        Bin size to consider when merging. All bins spaced within half
        a bin size are merged.
    
    Returns
    --------
    NDArray[int]:
        Merged values of pos, cnt_arr, len_arr
    List[str] or None:
        if not None, list of inserted sequences for each insertion position 
    
    """
    
    if len(pos) != len(cnt_arr) | len(pos) != len(len_arr):
        raise ValueError('Arguments must of same shape.')
    
    prev_cnt = list(cnt_arr)
    prev_len = list(len_arr)
    prev_pos = list(pos)

    curr_cnt = []
    curr_len = []
    curr_pos = []
    
    is_insert = (seq_list is not None)
    if is_insert:
        prev_seq = list(seq_list)
        curr_seq = []
    
    first_iter = True
    while prev_pos != curr_pos:  # Merge until bins are far enough apart
        if not first_iter:
            prev_cnt = 1*curr_cnt
            prev_len = 1*curr_len
            prev_pos = 1*curr_pos
            if is_insert:
                prev_seq = 1*curr_seq
        else:
            first_iter = False
    
        i = 0    

        curr_cnt = []
        curr_len = []
        curr_pos = []
        curr_seq = []
        while i < len(prev_pos):  
            i_next = np.searchsorted(prev_pos, prev_pos[i] + bin_size, side='right')  # Finidng bins to merge
            pos_mid = np.rint(np.mean(prev_pos[i: i_next])).astype(np.int32)
            curr_pos.append(pos_mid)
            curr_cnt.append(np.sum(prev_cnt[i: i_next]))
            curr_len.append(np.sum(prev_len[i: i_next]))
            if is_insert:  # Pooling sequences from bins that are merged
                curr_seq.append(list(chain(*prev_seq[i: i_next])))
            i = i_next
#             print(curr_seq)
        
    
    if is_insert:  # Find inserted sequences
        return_seq = []
        tmp_ip_file = 'tmp_seq.fa'
        for list_seq in curr_seq:
            if len(list_seq) > 1:  # Finding consensus sequence by MSA
                write_fasta(list_seq, tmp_ip_file)

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
                return_seq.append(''.join(seq).upper())
            elif len(list_seq) == 1:  # Only on inserted sequence
                return_seq.extend(list_seq)

        return np.array(curr_pos), np.array(curr_cnt), np.array(curr_len), return_seq
    
    else:
        return np.array(curr_pos), np.array(curr_cnt), np.array(curr_len), None


def detect_indel(bridge_cnt_arr: NDArray[int], indel_list: List[Union[Tuple[int, int], Tuple[int, int, str]]],
                 clip_list: Tuple[int, int, int], bin_size: int = 5, sig_frac: float = 0.5,
                 is_insert: bool = False
                ) -> List[Tuple[int, int, Optional[str]]]:
    
    """
    
    This functions infers the start positions and lengths of indels based on the number of reads
    bridging each position in the genome, and the start positions and lengths of indels and clippings
    deduced from the CIGAR strings.
    
    Parameters
    ----------
    bridge_cnt_arr: NDArray[int]
        The number of reads bridging each position in the genome
    indel_list: List[Union[Tuple[int, int], Tuple[int, int, str]]]
        List of (indel start position, indel length) or (indel start position, indel length, insert sequence)
    clip_list: Tuple[int, int, int]
        List of (clipping position, clipping length, direction of clipping).
        -1 indicates clipping ends at given position, +1 indicates that clipping starts at that position
    is_insert: bool
        True if looking for an insertion
    bin_size: int
        The range about each indel position to search for when merging indels. This must be odd.
        
    Returns
    --------
    List[Tuple[int, int]] or List[Tuple[int, int, str]]
        List of detected indels, namely their start positions and lengths as tuples

    """
    
    if not indel_list:  # Input indel list is empty
        return []

    indel_cnt_arr = np.zeros_like(bridge_cnt_arr, dtype=np.int32)  # Number of insertions/deletions starting at each pos
    indel_len_arr = np.zeros_like(bridge_cnt_arr, dtype=np.int32)  # Corresponding lengths of indels
    clip_cnt_arr = np.zeros_like(bridge_cnt_arr, dtype=np.int32)  # Number of reads with clippings covering each pos
    
    if is_insert:
        indel_pos_list, indel_lens, indel_seq_list = (np.array(z) for z in zip(*indel_list))
        indel_pos, indel_idx, indel_cnt = np.unique(indel_pos_list, return_inverse=True, return_counts=True)
        
        indel_seqs = []
#         print('PRINTING START')
        for p in indel_pos:
            indel_seqs.append(indel_seq_list[indel_pos_list == p])
#             print(p, indel_seq_list[indel_pos_list == p])
    else:
        indel_pos_list, indel_lens = (np.array(z) for z in zip(*indel_list))
        indel_pos, indel_idx, indel_cnt = np.unique(indel_pos_list, return_inverse=True, return_counts=True)
    
    indel_len = indel_len_arr[indel_pos]


    np.add.at(indel_len, indel_idx, indel_lens)  # Total length of insertions starting at that position

    # Aggregating nearby indels
    if is_insert:
        indel_pos, indel_cnt, indel_len, indel_seqs = merge_bins(
            indel_pos, indel_cnt, indel_len, indel_seqs, bin_size=bin_size)
    else:
        indel_pos, indel_cnt, indel_len, _ = merge_bins(indel_pos, indel_cnt, indel_len, None, bin_size=bin_size)
    indel_cnt_arr[indel_pos] = indel_cnt  # Number of insertion starting at that position
    indel_len_arr[indel_pos] = indel_len  # Number of insertion starting at that position

    clip_pos, clip_lens, clip_dirn = (np.array(z)
                                      for z in zip(*clip_list))  # Finding number of clipped bases that cover position
    for i, i_p in enumerate(indel_pos):
        clip_cnt_arr[i_p] = len(np.where((i_p >= clip_pos) & (i_p < clip_pos + clip_dirn*clip_lens))[0])
    
    indel_pos_sig = np.nonzero((indel_cnt_arr + clip_cnt_arr) > sig_frac*bridge_cnt_arr)[0]  # Finding significant indels
    indel_len = np.rint(indel_len_arr[indel_pos_sig]/indel_cnt_arr[indel_pos_sig]).astype(np.int32)  # Lenghts of selected indels
    
    if is_insert:  # Recover inserted sequences and extend to inferred length
        indel_seqs = [s for i, s in enumerate(indel_seqs) if indel_pos[i] in indel_pos_sig]
        return_seqs = []
        for s, l in zip(indel_seqs, indel_len):
            if l > len(s):  # Inferred sequence is smaller than that computed
                return_seqs.append(s + '-'*(l - l))
            else:  # Inferred sequence is longer than that computed
                return_seqs.append(s[:l])
        return list(zip(indel_pos_sig, indel_len, return_seqs))
    else:
        return list(zip(indel_pos_sig, indel_len, [None]*len(indel_len)))


def find_indel_longread(filehead: str, read_idx_dict: Dict[Tuple[str, int], List[int]],
                        rec_genomes: NDArray[(Any, Any), int],
                        ins_min : int = 10, del_min : int = 10, bin_size: int = 5,
                        sig_frac : float = 0.5
                       ) -> Tuple[List[List[Tuple[int, int, Optional[str]]]]]:
    
    """
    This function finds indels from long read data. It parses the SAM file and finds positions
    indicated as having indels by the CIGAR strings. It subsequently identifies 'significant'
    indels as those having a minimum length (can be set to different values for insertions and
    deletions) and having at least sig_frac number of reads covering that position exhibiting
    either a clipping or insertion/deletion at that position.
    
    Parameters
    ----------
    filehead: str
        Output folder where all XHap-related files are stored
    read_idx_dict:  Dict[Tuple[str, int], List[int]]
        Keyed by (QNAME, POS) with value = [row number of read in read-SNP matrix,
        index of genome for read]. (QNAME, POS) serves as a unique identifier for
        each read
    rec_genomes: NDArray[(Any, Any), int]
        Recovered haplotypes (as genomes)        
    ins_min: int
        Minimum length of insertions to be considered
    del_min: int
        Minimum length of deletions to be considered
    sig_frac: float
        Minimum fraction of reads at a position that have to exhibit indel or cllipping
        for indel to be considered significant
    
    Returns
    -------    
    List[List[Tuple[int, int, str]]]
        List of lists, where each list is the insertions found in that genome 
    
    List[List[Tuple[int, int, None]]]
        List of lists, where each list is the deletions found in that genome 
    
    Indels are returned as tuples of (start reference position, length of indel)

    
    """
                
    samfile = pysam.AlignmentFile("generate_data/" + filehead + "/" + filehead +".sam", "r")
    sam_iter = samfile.fetch()

    len_ins_min = ins_min  # Insertion must be at least 10bp
    len_del_min = del_min  # Deletion must be at least 10bp
    indel_bin_size = bin_size  # Size of bins when determining indels (must be odd)

    bridge_cnt_arr = np.zeros_like(rec_genomes, dtype=np.int32)
    del_list = [[] for _ in range(np.shape(rec_genomes)[0])]
    ins_list = [[] for _ in range(np.shape(rec_genomes)[0])]
    clip_list = [[] for _ in range(np.shape(rec_genomes)[0])]

    for i, x in enumerate(sam_iter):
        if not x.is_unmapped:
            rname = x.query_name
            pos = x.reference_start  # Position
            genome_idx = read_idx_dict[(rname, pos)][1]  # Genome the read is attributed to

            c_arr, clen_arr = tuple(np.array(z) for z in zip(*x.cigartuples))
            c_mask = (c_arr == 1)  # Masking alignements that do not consume reference base
            for c in [4, 5, 6]:
                c_mask[c_arr == c] = True
            clen_arr_ref = 1*clen_arr    
            clen_arr_ref[c_mask] = 0  # These alignments do not consume reference base
            local_refpos = np.cumsum(np.insert(clen_arr_ref, 0, 0))[:-1]
            
            c_mask = (c_arr == 2)  # Masking elements that do not consume query base
            for c in [3, 5, 6]:
                c_mask[c_arr == c] = True
            clen_arr_query = 1*clen_arr
            clen_arr_query[c_mask] = 0  # These alignments do not consume query base
            local_querypos = np.cumsum(np.insert(clen_arr_query, 0, 0))[:-1]

            pos_del = pos + local_refpos[c_arr == 2]  # Positions of deletions along reference
            len_del = clen_arr[c_arr == 2]  # Length of deletions
            sel_del = (len_del > len_del_min)
            del_list[genome_idx].extend(zip(pos_del[sel_del], len_del[sel_del]))
            
            pos_ins = pos + local_refpos[c_arr == 1]  # Positions of insertions along reference
            len_ins = clen_arr[c_arr == 1]
            sel_ins = (len_ins > len_ins_min)
            pos_query = local_querypos[c_arr == 1]  # Positions of insertions along query
            
            def seq_ins_gen(pos_query, len_ins):
                for p, l in zip(pos_query, len_ins):
                    yield x.query_sequence[p: p + l]

            ins_list[genome_idx].extend(zip(pos_ins[sel_ins], len_ins[sel_ins],
                                            seq_ins_gen(pos_query[sel_ins], len_ins[sel_ins])))

            mask_clip = (c_arr == 4) | (c_arr == 5)
            pos_clip = pos + local_refpos[mask_clip]
            len_clip = clen_arr[mask_clip]
            dirn_clip = (pos_clip == pos)*(-1) + (pos_clip != pos)*1  # Storing if clipped segment starts or ends at pos_clip
            clip_list[genome_idx].extend(zip(pos_clip, len_clip, dirn_clip))

            bridge_start = x.reference_start + 1
            bridge_end = x.reference_end - 1
            bridge_cnt_arr[genome_idx, bridge_start: bridge_end] = (
                bridge_cnt_arr[genome_idx, bridge_start: bridge_end] + 1)
    
    del_res = [[] for _ in range(np.shape(rec_genomes)[0])]
    ins_res = [[] for _ in range(np.shape(rec_genomes)[0])] 
    
    for g_i in range(np.shape(rec_genomes)[0]):  # For each haplotype
        ins_res[g_i].extend(detect_indel(bridge_cnt_arr[g_i], ins_list[g_i], clip_list[g_i], bin_size, sig_frac,
                                         is_insert=True))
        del_res[g_i].extend(detect_indel(bridge_cnt_arr[g_i], del_list[g_i], clip_list[g_i], bin_size, sig_frac))
    
    return ins_res, del_res