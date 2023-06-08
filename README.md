# XHap: Haplotype assembly sing long-distance read correlations learned by transformers

## About
XHap is a framework based on transformers (from natural language processing) for diploid and polyploid haplotype assembly. The framework is capable of haplotype reconstruction from either short reads (e.g. Illumina, Roche) or long reads (e.g. PacBio, ONT) or a combination thereof.

The current implementation of XHap uses Python3, PyTorch and PyTorch ROCm (for AMD GPUs). Both CPU and GPU implementations are available in _xhap.py_ and _xhap_parallel.py_ respectively.

## Dependencies
- PyTorch >= 1.10
- PyTorch ROCm >= 2.0.1 (to use AMD GPUs)
- Numpy
- Scipy
- C++
- Samtools
- MAFFT

Where possible, additional dependencies have been included in the GitHub repository.

## Assumed directory structure
All the scripts included in this repository assume that the XHap source code is stored in the current working directory and the data files are stored in a subdirectory [data] of the directory _generate_data_. Consequently, the resulting data files can be found in _generate_data/[data]_.

> **Note:** This structure can be easily changed in the provided scripts by changing _generate_data_ to the desired folder in the respective files. 

## Input
The provided pipeline for XHap takes a tab-separated file containing read-SNP matrix as the input. Details on how to obtain this matrix can be found in the **generate_data** folder.

## Output
Each round of training XHap yields results stored in the following files saved in the corresponding data directory:
- **xhap_model:** Stores the state_dict for both the convolutional (embedAE) and transformer encoder (corr_xformer) layers in XHap.
- **haptest_xformer_res.npz:** NPZ file storing the reconstructed haplotypes (rec_hap), the read attributions (rec_hap_origin) and if applicable, the ground truth haplotyes (true_hap).

## Usage
The function _train_xhap_ in xhap (or xhap_parallel) can be invoked to run XHap on the data in the folder specified by _outhead_. This function also takes in the following parameters:
- **d_model**: Embeding size for each read
- **num_hap**: Number of haplotypes (ploidy of organism)
- **num_epoch**: Number of training epochs for XHap
- **check_cpr**: Set to true if the ground truth is present

## Included scripts
There are several Python scripts included to run XHap end-to-end -- from data generation through haplotype assembly. These scripts can also be used to replicate the experiments described in the associated manuscript.

1. **run_expt.py:** Running experiments on semi-experimental data (includes data generation)
 1. **run_expt_real.py:** Running experiments on experimental data (assumes data processing has been done)
1. **run_indel_expt.py:** Running experiments to validate indel detection included in the XHap pipeline

## Citation