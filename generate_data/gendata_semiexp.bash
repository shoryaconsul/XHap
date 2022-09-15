#!/bin/bash

HAPGEN="HaploSim/haplogenerator.py" # Path to haplotypegenerator.py
ART="art_bin_ChocolateCherryCake/art_illumina"  # Path to ART
BWA="bwa-0.7.17/bwa"  # Path to BEA aligner
FREEBAYES="./freebayes-1.3.6"  # Path to FreeBayes
EXTRACT_MATRIX="./ExtractMatrix"  # Path to ExtractMatrix executable
CAECSEQ="./CAECseq_haplo.py"  # Path to CAECSeq python file
HAPCUT="./HapCUT2"  # Folder containing HapCUT2 and htslib source files
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
  python2 $HAPGEN -f $REF -o $OUTHEAD --model $MODEL -s "[4.63,0,0]" --sdlog "[0.693,0,0]" -m "{'A':'GCT','G':'ACT','C':'TAG','T':'ACG'}" -p $HAPNUM
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
  $ART -p -na -i "$HAPFOLDER/combined.fa" -l 150 --seqSys HS25 -f 5 -m 550 -s 5 -o "$HAPFOLDER/hap"
else   
  $ART -p -na -q -i "$HAPFOLDER/combined.fa" -l 150 --seqSys HS25 -f 5 -m 550 -s 5 -o "$HAPFOLDER/hap"
fi
rm "$HAPFOLDER/"*".aln"
echo "Generated reads"

# Align reads
RG_HEADER="@RG\tID:rg1\tSM:sample1"
$BWA index $REF
$BWA mem -t 5 -R $RG_HEADER $REF "$HAPFOLDER/hap"*".fq" > "$HAPFOLDER/$OUTHEAD.sam"  # t is number of threads
samtools sort "$HAPFOLDER/$OUTHEAD.sam" -o "$HAPFOLDER/$OUTHEAD""_sorted.bam"

# Generate read-SNP matrix
THRESH=0.2
HAPLEN=$(awk '{print length }' $REF | tail -1)  # Find length of reference genome
$EXTRACT_MATRIX -f $REF -s "$HAPFOLDER/$OUTHEAD.sam" -t $THRESH -z "$HAPFOLDER/$OUTHEAD" -k $HAPNUM \
-b 0 -e $HAPLEN -q 0 -l 100 -i 560 

# CAECSeq
python $CAECSEQ "$HAPFOLDER/$OUTHEAD" -k $HAPNUM -n 3

# Variant calling
$FREEBAYES -f $REF -p $HAPNUM "$HAPFOLDER/$OUTHEAD""_sorted.bam" > "$HAPFOLDER/$OUTHEAD""_variants.vcf"

# Generate SNP fragment matrtix for SDHap/HapCUT
cd $HAPCUT
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$(pwd)/htslib"  # Add to path
cd ..
$EXTRACTHAIRS --VCF "$HAPFOLDER/$OUTHEAD""_variants.vcf" --bam "$HAPFOLDER/$OUTHEAD""_sorted.bam" --maxIS 3000 --out "$HAPFOLDER/$OUTHEAD""_fragment_file.txt"

: << 'COMMENT'

B=$(echo | awk -v c=$C -v d=$D '{print c / d}')

COMMENT
