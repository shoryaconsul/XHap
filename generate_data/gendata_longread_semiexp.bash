#!/bin/bash

HAPGEN="HaploSim/haplogenerator.py" # Path to haplotypegenerator.py
PBSIM="./pbsim2/src/pbsim"  # Path of PBSIM2 executable
BWA="bwa-0.7.17/bwa"  # Path to BEA aligner
FREEBAYES="./freebayes-1.3.6"  # Path to FreeBayes
EXTRACT_MATRIX="./ExtractMatrix_longread"  # Path to ExtractMatrix executable
CAECSEQ="./CAECseq_haplo.py"  # Path to CAECSeq python file
HAPCUT="./HapCUT2"  # Folder containing HapCUT2 and htslib source files
EXTRACTHAIRS="./HapCUT2/build/extractHAIRS"  # Path to extractHAIRs executable from HapCUT2
VERBOSE=false
INDEL=100  # Length of insertion/deletion

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
  python2 $HAPGEN -f $REF -o $OUTHEAD --model $MODEL -s "[3.03,0,0]" --sdlog "[1.293,0,0]" -m "{'A':'GCT','G':'ACT','C':'TAG','T':'ACG'}" -p $HAPNUM
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

# Add indels to one of the genomes
python3 add_indel.py -f "$HAPFOLDER/combined.fa" -i $INDEL

# Generate long reads using PBSIM2
# length_* paramters are related to read length
$PBSIM --depth 13.3 --prefix "$HAPFOLDER/$OUTHEAD" --hmm_model pbsim2/data//P6C4.model "$HAPFOLDER/combined.fa"
for i in `ls "$HAPFOLDER/$OUTHEAD"_*.fastq`; do cat $i >> "$HAPFOLDER/$OUTHEAD.fq"; done
# rm "$HAPFOLDER/$OUTHEAD"_*.fastq "$HAPFOLDER/$OUTHEAD"_*.ref "$HAPFOLDER/$OUTHEAD"_*.maf


# Align reads
# RG_HEADER="@RG\tID:rg1\tSM:sample1"
$BWA index $REF
$BWA mem -t 5 $REF "$HAPFOLDER/$OUTHEAD.fq" > "$HAPFOLDER/$OUTHEAD.sam"  # t is number of threads
samtools sort "$HAPFOLDER/$OUTHEAD.sam" -o "$HAPFOLDER/$OUTHEAD""_sorted.bam"
echo "Generated reads"

# Generate read-SNP matrix
THRESH=0.2
HAPLEN=$(awk '{print length }' $REF | tail -1)  # Find length of reference genome
$EXTRACT_MATRIX -f $REF -s "$HAPFOLDER/$OUTHEAD.sam" -t $THRESH -z "$HAPFOLDER/$OUTHEAD" -k $HAPNUM \
-b 0 -e $HAPLEN -q 0 -l 100 -i 560 

# CAECSeq
python $CAECSEQ "$HAPFOLDER/$OUTHEAD" -k $HAPNUM -n 5

: << 'COMMENT'
# Variant calling
$FREEBAYES -f $REF -p 2 "$HAPFOLDER/$OUTHEAD""_sorted.bam" > "$HAPFOLDER/$OUTHEAD""_variants.vcf"

# Generate SNP fragment matrtix for SDHap/HapCUT
cd $HAPCUT
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$(pwd)/htslib"  # Add to path
cd ..
$EXTRACTHAIRS --VCF "$HAPFOLDER/$OUTHEAD""_variants.vcf" --bam "$HAPFOLDER/$OUTHEAD""_sorted.bam" --maxIS 3000 --out "$HAPFOLDER/$OUTHEAD""_fragment_file.txt"


B=$(echo | awk -v c=$C -v d=$D '{print c / d}')

COMMENT
