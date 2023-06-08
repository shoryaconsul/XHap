# Data simulation and benchmarking
This folder contains all the code for the algorithms used to benchmark the performance of XHap, as well as the tools needed to simulate semi-experimental and experimental data generation.

## Included bencmakrs
1. [H-PoP](https://github.com/MinzhuXie/H-PoPG)
2. [HapCUT2](https://github.com/vibansal/HapCUT2)
3. [HapTree](http://cb.csail.mit.edu/cb/haptree/)
4. [Ranbow](https://github.com/moeinzadeh/Ranbow)

All the abovementioned packages but HapCUT2 have been included in this repository. HapCUT2 requires the htslib library to be installed. The Bash script _install_hapcut2.sh_ has been included to enable easy installation of HapCUT2 in this directory.

## Data simulation
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