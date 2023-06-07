#!/bin/bash

BWA="bwa-0.7.17/bwa"  # Path to BEA aligner
EXTRACT_MATRIX="./ExtractMatrix"  # Path to ExtractMatrix executable
CAECSEQ="./CAECseq_haplo.py"  # Path to CAECSeq python file
EXTRACTHAIRS="./HapCUT2/build/extractHAIRS"  # Path to extractHAIRs executable from HapCUT2
VERBOSE=false

while getopts f:o:n:v FLAG
do
  case "${FLAG}" in
    f) REF=${OPTARG};;  # Path to reference genome
    o) OUTHEAD=${OPTARG};;  # Base name of generated haplotype files
    n) HAPNUM=${OPTARG};;  # Number of haplotypes to generate
    v) VERBOSE=true;;
    *) echo "Invalid command line option: -$FLAG" ;;
  esac
done


if [[ ! -f $REF ]];then
  echo "Reference genome $REF not found"
  exit 1
fi

if [[ -f "$REF.fai" ]];then
  echo "Reference genome index found"
else
  samtools faidx $REF  # Index reference genome if not already done
fi

HAPFOLDER="./$OUTHEAD"
if [[ -d $HAPFOLDER ]];then
  rm -r $HAPFOLDER
fi

# Move generate haplotype files into new folder
mkdir -p ./tmp
mv "$HAPFOLDER"* ./tmp
mv ./tmp $HAPFOLDER
rm "$HAPFOLDER/"*".fai"  # Remove extraneous faidx files

# Align reads
$BWA index $REF
$BWA mem -t 5 $REF "SRR6173308/SRR6173308_"*".fastq" | samtools view -h -F 4 > "$HAPFOLDER/$OUTHEAD.sam" # t is number of threads
# samtools view -h -F 4 "$HAPFOLDER/$OUTHEAD.sam" > "$HAPFOLDER/$OUTHEAD.sam"  # Remove unmapped reads
samtools sort "$HAPFOLDER/$OUTHEAD.sam" -o "$HAPFOLDER/$OUTHEAD""_sorted.bam"


# Generate read-SNP matrix
THRESH=0.1
HAPLEN=$(awk '{print length }' $REF | tail -1)  # Find length of reference genome
$EXTRACT_MATRIX -f $REF -s "$HAPFOLDER/$OUTHEAD.sam" -t $THRESH -z "$HAPFOLDER/$OUTHEAD" -k $HAPNUM \
-b 0 -e $HAPLEN -q 0 -l 30 -i 560 

# Create VCF file from read-SNP matrix
python snp2vcf.py -r $REF -p "$HAPFOLDER/$OUTHEAD""_SNV_pos.txt" -m "$HAPFOLDER/$OUTHEAD""_SNV_matrix.txt" -o "$HAPFOLDER/$OUTHEAD""_variants.vcf"

: << 'COMMENT'
# CAECSeq
python $CAECSEQ "$HAPFOLDER/$OUTHEAD" -k $HAPNUM -n 5

COMMENT

