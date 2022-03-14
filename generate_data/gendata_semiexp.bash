#!/bin/bash

HAPGEN="HaploSim/haplogenerator.py" # Path to haplotypegenerator.py
ART="art_bin_ChocolateCherryCake/art_illumina"  # Path to ART
BWA="bwa-0.7.17/bwa"  # Path to BEA aligner
FREEBAYES="./freebayes-1.3.6"  # Path to FreeBayes
EXTRACT_MATRIX="./ExtractMatrix"  # Path to ExtractMatrix executable
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

: << 'COMMENT'
REF=$2  # Path to reference genome
OUTHEAD=$3  # Base name of generated haplotype files
HAPNUM=$4  # Number of haplotypes to generate
COMMENT

MODEL="lognormal"  # Model used by haplotype generator (poisson or lognormal)


if [[ ! -f $REF ]];then
  echo "Reference genome $REF not found"
  exit 1
fi

if [[ -f "$REF.fai" ]];then
  echo "Reference genome index found"
else
  samtools faidx $REF  # Index reference genome if not already done
fi

# Generate haplotypes
if [[ $MODEL == "poisson" ]];then
  python2 $HAPGEN -f $REF -o $OUTHEAD --model $MODEL -s "[0.04,0.01,0.01]" -m "{'A':'GCT','G':'ACT','C':'TAG','T':'ACG'}" -p $HAPNUM
elif [[ $MODEL == "lognormal" ]];then
  python2 $HAPGEN -f $REF -o $OUTHEAD --model $MODEL -s "[2.70,0,0]" --sdlog "[0.693,0,0]" -m "{'A':'GCT','G':'ACT','C':'TAG','T':'ACG'}" -p $HAPNUM
else
  echo "Invalid model specified"
  exit 1
fi
echo "Haplotypes generation complete"

HAPFOLDER="./$OUTHEAD"
if [[ -d $HAPFOLDER ]];then
  rm -r $HAPFOLDER
fi

# Move generate haplotype files into new folder
mkdir -p ./tmp
mv "$HAPFOLDER"* ./tmp
mv ./tmp $HAPFOLDER
rm "$HAPFOLDER/"*".fai"  # Remove extraneous faidx files

cat "$HAPFOLDER/"*".fa" > "$HAPFOLDER/combined.fa"  # Combined FASTA files


# Generate paired-end reads reads using ART
# -l [read length] -f [coverage] -m [mean insert length] -s [sd of insert length]
if [[ $VERBOSE = true ]];then
  $ART -p -na -i "$HAPFOLDER/combined.fa" -l 250 -f 60 -m 550 -s 5 -o "$HAPFOLDER/haptest"
else   
  $ART -p -na -q -i "$HAPFOLDER/combined.fa" -l 250 -f 60 -m 550 -s 5 -o "$HAPFOLDER/haptest"
fi
rm "$HAPFOLDER/"*".aln"
echo "Generated reads"

# Align reads
RG_HEADER="@RG\tID:rg1\tSM:sample1"
$BWA index $REF
$BWA mem -t 5 -R $RG_HEADER $REF "$HAPFOLDER/"*".fq" > "$HAPFOLDER/$OUTHEAD.sam"  # t is number of threads
samtools sort "$HAPFOLDER/$OUTHEAD.sam" -o "$HAPFOLDER/$OUTHEAD""_sorted.bam"

# Generate read-SNP matrix

# Variant calling
$FREEBAYES -f $REF -p 2 "$HAPFOLDER/$OUTHEAD""_sorted.bam" > "$HAPFOLDER/$OUTHEAD.vcf"

: << 'COMMENT'
COMMENT
