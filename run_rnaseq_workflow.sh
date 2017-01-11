#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -N Workflow
#$ -e logs/error_$TASK_ID.txt
#$ -o logs/stdout_$TASK_ID.txt
#$ -pe threaded 8
#$ -l mem_free=16G

# Source conda environment
source activate rnaseq_workflow


# - - - - - - - - - - -
# Variables to set
# - - - - - - - - - - -
qualityCutoff=20
trimLength=50


# - - - - - - - - - - -
# Array job variables
# - - - - - - - - - - -
files=(input/*.fastq.gz)
number=`expr $SGE_TASK_ID - 1`

# Path to current working file
sequence=${files[number]}

# Name of current working file
BN=`basename $sequence .fastq.gz`


# - - - - - - - - - - - -
# Verify workflow can run
# - - - - - - - - - - - -

# No Trim Galore! in PATH
if [[ -z $(type STAR) ]]; then
    echo -e "No STAR installation found! Exiting now."
    exit 1

# No Trim Galore! in PATH
elif [[ -z $(type trim_galore) ]]; then
    echo -e "No TrimeGaore! installation found! Exiting now."
    exit 1

# No cutadapt in PATH
elif [[ -z $(type cutadapt) ]]; then
    echo -e "No cutadapt installation found! Exiting now."
    exit 1

# No sortmerna in PATH
elif [[ -z $(type sortmerna) ]]; then
    echo -e "No SortMeRNA installation found! Exiting now."
    exit 1

# No featureCounts in PATH
elif [[ -z $(type featureCounts) ]]; then
    printf "No featureCounts installation found! Exiting now."
    exit 1

# No FASTQ in PATH
elif [[ -z $(type fastqc) ]]; then
    echo -e "No FASTQC installation found! Exiting now."
    exit 1

# No Index created
elif [[ -z "$(ls -A index/)" ]]; then
    echo -e "Error genome index is empty. You must first generate an index before alignment."
    exit 1
fi


# - - - - - - - - - - - - - - - - -
# Print cluster info
# - - - - - - - - - - - - - - - - -
echo -e "- - - - Diagnostics - - - - -"
echo -e "Number of slots: $NSLOTS"
echo -e "Number of hosts: $NHOSTS"
echo -e "Number in Queue: $QUEUE"
echo -e "OS type: $SGE_ARCH \n"
echo -e "Run parameters:"
echo -e "Host name: "$(hostname -f)" "
echo -e "Date started: $(date)"
echo -e "Currently in:" $(pwd)
echo -e "Python version: "$(python -V 2>&1)" "
echo -e "Working on SGE#: "${SGE_TASK_ID}" "
echo -e "Working on file: "${files[number]}""
echo -e "- - - - - - - - - - - - - - - \n"


# - - - - - - - - - - - -
# Set default variables
# - - - - - - - - - - - -

# Input/Output
outputFolder=$(pwd)/"output"
mkdir -p $outputFolder

# Genome/Annotation
genomeIndexDir=$(pwd)/"index/"

# QC data
outputQcFolder="${outputFolder}/1_initial_qc"
outputTrimFolder="${outputFolder}/2_trimmed_output"

# SortMeRNA
sortMeRnaAligned="${outputFolder}/3_rRNA/aligned/"
sortMeRnaFiltered="${outputFolder}/3_rRNA/filtered/"
sortMeRnaLogs="${outputFolder}/3_rRNA/logs/"
sortmernaDB="tools/sortmerna-2.1"

# Alignment
alignedSequences="${outputFolder}/4_aligned_sequences"
alignedBAM="${outputFolder}/4_aligned_sequences/aligned_bam/"
alignedLog="${outputFolder}/4_aligned_sequences/aligned_logs/"
#alignedStat="${outputFolder}/4_aligned_sequences/aligned_stats/"
alignedCounts="${outputFolder}/4_aligned_sequences/aligned_counts/"


# - - - - - - - - - - -
# Run FASTQC
# - - - - - - - - - - -
echo -e "Quality analysis of raw reads... \n"
mkdir -p $outputQcFolder

# Run Quality Analysis on Raw Reads
fastqc \
-o $outputQcFolder \
--noextract \
-t $NSLOTS \
$sequence


# - - - - - - - - - - -
# Run Trim Galore!
# - - - - - - - - - - -
echo -e "Trimming raw reads... \n"
mkdir -p $outputTrimFolder

# Run Trim Galore to remove adapters and low base quality scores
trim_galore \
--quality $qualityCutoff \
--fastqc \
--length $trimLength \
--output_dir $outputTrimFolder \
$sequence

# unzip all sequences (if needed)
gunzip $outputTrimFolder/${BN}_trimmed.fq.gz


# - - - - - - - - - - -
# Run SortMeRNA
# - - - - - - - - - - -
echo -e "Removing rRNA sequences... \n"
mkdir -p $sortMeRnaAligned
mkdir -p $sortMeRnaFiltered
mkdir -p $sortMeRnaLogs

# Initialize Database
sortmernaREF=${sortmernaDB}/rRNA_databases/silva-arc-16s-id95.fasta,${sortmernaDB}/index/silva-arc-16s-id95:\
${sortmernaDB}/rRNA_databases/silva-arc-23s-id98.fasta,${sortmernaDB}/index/silva-arc-23s-id98:\
${sortmernaDB}/rRNA_databases/silva-bac-16s-id90.fasta,${sortmernaDB}/index/silva-bac-16s-id95:\
${sortmernaDB}/rRNA_databases/silva-bac-23s-id98.fasta,${sortmernaDB}/index/silva-bac-23s-id98:\
${sortmernaDB}/rRNA_databases/silva-euk-18s-id95.fasta,${sortmernaDB}/index/silva-euk-18s-id95:\
${sortmernaDB}/rRNA_databases/silva-euk-28s-id98.fasta,${sortmernaDB}/index/silva-euk-28s-id98

# sortmerna
sortmerna \
--ref $sortmernaREF \
--reads $outputTrimFolder/${BN}_trimmed.fq \
--aligned ${sortMeRnaAligned}/${BN}_aligned \
--other ${sortMeRnaFiltered}/${BN}_filtered \
--fastx \
--log \
-a $NSLOTS \
-v

# Move Log Files into correct order
mv -v ${sortMeRnaAligned}/${BN}_aligned.log $sortMeRnaLogs

# Removed cutadapt Trimmed Files to save space
rm -rf $outputTrimFolder/${BN}_trimmed.fq

# Gzip the sequences that aligned to save space
gzip ${sortMeRnaAligned}/${BN}_aligned.fq


# - - - - - - - - - - -
# Run STAR-aligner
# - - - - - - - - - - -
echo -e "Aligning sequences... \n"
mkdir -p $alignedSequences
mkdir -p $alignedBAM
mkdir -p $alignedLog
mkdir -p $alignedStat
mkdir -p $alignedCounts

# STAR
STAR \
--genomeDir $genomeIndexDir \
--readFilesIn ${sortMeRnaFiltered}/${BN}_filtered.fq \
--runThreadN $NSLOTS \
--outFileNamePrefix $alignedSequences/${BN} \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts


# - - - - - - - - - - - - - -
# Move and clean up
# - - - - - - - - - - - - - -

# Print progress message
echo -e "Moving files around and cleaning up... \n"

# Gzip the sequences that were filtered to save space
gzip ${sortMeRnaFiltered}/${BN}_filtered.fq

# Move Sorted BAM
mv -v $alignedSequences/${BN}Aligned.sortedByCoord.out.bam $alignedBAM/

# Rename Sorted BAM
mv -v $alignedBAM/${BN}Aligned.sortedByCoord.out.bam $alignedBAM/${BN}.bam

# Move Final Log.out file
mv -v $alignedSequences/${BN}Log.final.out $alignedLog/

# Move Log.out file
mv -v $alignedSequences/${BN}*Log.out $alignedLog/

# Move (.Tab) Files
mv -v $alignedSequences/${BN}*out.tab $alignedCounts/

# Remove leftover TEMP folders
rm -rf $alignedSequences/${BN}*_STARtmp

# Remove *Log.progress.out file
rm -rf $alignedSequences/${BN}*Log.progress.out


# End of script
echo -e "RNAseq completed! \n"
echo -e "Date finished: $(date)"
deactivate
