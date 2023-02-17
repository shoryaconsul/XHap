#!/bin/bash

RANBOW="ranbow/ranbow.py"  # Path to Ranbow Python file
PAR="ranbow_params"  # Path to Ranbow params file
SCAF="ranbow_scaffold.list"  # Path to Ranbow scaffolds file

NUMPROC=1  # Number of processors

while getopts n:b:r:v:o:w: FLAG
do
  case "${FLAG}" in
    n) PLOIDY=${OPTARG};;  # Ploidy
    b) BAMFILE=${OPTARG};;  # Sorted BAM file
    r) REFFILE=${OPTARG};;  # Reference FASTA file
    v) VCFFILE=${OPTARG};;  # Variants VCF file
    o) OUTBASE=${OPTARG};;  # Output folder base + prefix
    w) WINLEN=${OPTARG};;  # Window length for Ranbow
    *) echo "Invalid command line option: -$FLAG" ;;
  esac
done

awk '{print $1}' "$REFFILE"".fai" > $SCAF

# Writing required parameters into Ranbow parameter file
echo "-ploidy $PLOIDY" > $PAR
echo "-noProcessor $NUMPROC" >> $PAR
echo "-bamFile $BAMFILE" >> $PAR
echo "-refFile $REFFILE" >> $PAR
echo "-vcfFile $VCFFILE" >> $PAR
echo "-selectedScf $SCAF" >> $PAR
echo "-outputFolderBase $OUTBASE" >> $PAR

python2 $RANBOW hap -mode index -par $PAR -processorIndex 0
python2 $RANBOW hap -mode hap -par $PAR -processorIndex 0 -WinLen $WINLEN
echo "Finished running RANBOW --------------------------------"
python2 $RANBOW hap -mode collect -par $PAR