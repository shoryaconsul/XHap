# Data simulation and benchmarking
This folder contains all the code for the algorithms used to benchmark the performance of XHap, as well as the tools needed to simulate semi-experimental and experimental data generation.

## Included bencmakrs
1. [H-PoP](https://github.com/MinzhuXie/H-PoPG)
2. [HapCUT2](https://github.com/vibansal/HapCUT2)
3. [HapTree](http://cb.csail.mit.edu/cb/haptree/)
4. [Ranbow](https://github.com/moeinzadeh/Ranbow)

All the abovementioned packages but HapCUT2 have been included in this repository. HapCUT2 requires the htslib library to be installed. The Bash script _install_hapcut2.sh_ has been included to enable easy installation of HapCUT2 in this directory.

## Data simulation

### Semi-experimental data
Semi-experimental data generation comprises the following steps:
1. Generating haplotypes from selected reference genome using _haplotypegenerator.py_ from [HaploSim](https://github.com/EhsanMotazedi/Haplosim).
2. Generating short or long reads using [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) or P[BSIM2](https://github.com/yukiteruono/pbsim2) respectively.
3. Aligning reads to reference genome using BWA-MEM.
4. Variant calling and generating read-SNP matrix using ExtractMatrix (or ExtractMatrix_longread for long reads).

The provided Bash script **gendata_semiexp.bash** can be used to perform the aforemnetioned steps. The script has to be supplied with the following arguments (those in bold are required): 

| Option | Description |
|--------|-------------|
| **-f** | Path to reference genome |
| **-o** | Output directory |
| **-n** | Ploidy |
| **-c** | Coverage per haplotype, i.e., coverage/n|
| -i | Lenght of insertion to be placed in one of the generated haplotype sequences |
| -d | Lenght of deletion to be placed in one of the generated haplotype sequences |
| -v | Added if more verbose output is desired (useful for debugging)|

The simulated files can be found in _generate_data/[out] where [out] is the value provided to the option -o. 
> Most of the files found in this subdirectory are prefixed with [out]. For example, if [out] is specifed to be 'test', then the data files in _generate_data/test_ will mostly be of the form _test*.*_.

### Experimental data
Similarly, the Bash script **gendata_exp.bash** can be used to process experimental data (stored in a subdirectory _SRR6173308_). It only takes the first three of the above arguments.

### Miscellaneous scripts
The script **run_caecseq_expt.py** is a convenient way to run CAECSeq on multiple datasets. This may be advisable as methods based on deep learning usually take longer than other algorithms.

**run_benchmarks.py** can be used to run the aforementioned algorithms on the specified data. The location of the various benchmarks can be altered by modifying the path specifications at the start of this Python script. The arguments for this script are listed below. 

| Option | Long option | Description |
|--------|-------------|-------------|
| -f | --filehead  |  Prefix of data files ([out])|
| -r | --reference |  Path to reference genome |
| -p | --ploidy| Ploidy of organism |
| -n | --num_expt | Number of datasets to run benchmarks over |
|-a | --algo_runs | Number of times to run each benchmark per dataset |
||--mec-only| Specify if only MEC is to be computed (usually when ground truth is absent) |
||--verbose| Specify for more verbose output (for debugging)


> When running over multiple datasets, it is assumed that the datasets are named as *[prefix]_iter[number]*. For example, if the prefix is 'test', then the data subdirectories should be named _test_iter1, test_iter2, etc._

### Included reference genomes
- **solanum_tuberosum.fa:** 10 kbp sample of Solanum Tuberosum Chromosome 5 genome (used for semi-experimental data generation)
- **human_sample.fa:**: 100 kbp sample of GrCh38 genome (used for semi-experimental data generation)
- **solanum_tuberosum_[*].fa:** 10 kbp samples of Solanum Tuberosum Chromosome 5 genome (used for experimental data processings)