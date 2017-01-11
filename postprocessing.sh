#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -N Output_Processing
#$ -e Post_error.txt
#$ -o Post_stdout.txt
#$ -pe threaded 8
#$ -l mem_free=12G

# - - - - - - - - - - - - - - - - - - - - - - - -
# @ Thomas W. Battaglia
# @ tb1280@nyu.edu

# Description:
# This is the workflow script for processing and
# analyzing RNAseq data from raw FASTQ sequences.
# There must be an index made for both the STAR
# aligner as well as the SortMeRNA removal tool.
# - - - - - - - - - - - - - - - - - - - - - - - -

# Source virtual environment that contains all
# the required tools.
source activate rnaseq_workflow


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
echo -e "- - - - - - - - - - - - - - - \n"


# - - - - - - - - - - - -
# Verify workflow can run
# - - - - - - - - - - - -

# No featureCounts in PATH
if [[ -z $(type featureCounts) ]]; then
    echo "No featureCounts installation found! Exiting now."
    exit 1

# No MultiQC in PATH
elif [[ -z $(type multiqc) ]]; then
    echo "No MultiQC installation found! Exiting now."
    pip install multiqc
fi


# - - - - - - - - - - - -
# Set default variables
# - - - - - - - - - - - -

# Input/Output
outputFolder=$(pwd)/"output"
mkdir -p $outputFolder

# Annotation file
annotationFile=$(pwd)/$(ls -d -1 annotation/*)

# Alignment
alignedSequences="${outputFolder}/4_aligned_sequences"
alignedBAM="${outputFolder}/4_aligned_sequences/aligned_bam/"

# featureCounts
finalCounts="${outputFolder}/5_final_counts"

# MultiQC
multiQc="${outputFolder}/6_multiQC"


# - - - - - - - - - - - - - -
# Run Subread (featureCounts)
# - - - - - - - - - - - - - -
echo -e "Summarizing gene counts... \n"
mkdir -p $finalCounts

# Change directory
cd output/4_aligned_sequences/aligned_bam

# Store list of files as a variable
dirlist=$(ls -t ./*.bam | tr '\n' ' ')

# Run featureCounts
featureCounts \
-a $annotationFile \
-o $finalCounts/final_counts.txt \
-g 'gene_name' \
-T $NSLOTS \
$dirlist

# Change back to main directory
cd ../../..


# - - - - - - - - - - - - - -
# Run MultiQC
# - - - - - - - - - - - - - -
echo "Running multiQC..."
mkdir -p $multiQc

# Run multiqc
multiqc $outputFolder \
--outdir $multiQc

# End of script
echo "Processing completed!"
echo "Date finished: $(date)"
deactivate
