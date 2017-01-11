#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -N STAR_Index
#$ -e star_index_error.txt
#$ -o star_index_stdout.txt
#$ -pe threaded 12
#$ -l mem_free=24G

# - - - - - - - - - - - - - - - - - - - - - - - - - - -
# @ Thomas W. Battaglia
# tb1280@nyu.edu

# This script will generate index files for STAR-aligner
# before use with RNAseq data alignment. This file can
# can be submitted as a job if on a proper cluster
# enviroment.
# - - - - - - - - - - - - - - - - - - - - - - - - - -

# Source virtualevn
source activate rnaseq_workflow


# - - - - - - - - - - - -
# Set variables
# - - - - - - - - - - - -
export genomeFile=$(ls -d -1 genome/*)
export annotationFile=$(ls -d -1 annotation/*)


# - - - - - - - - - - - -
# Verify script can run
# - - - - - - - - - - - -

# No STAR installation
if [[ -z $(STAR --version) ]]; then
    echo "No STAR installation found! Exiting now."
    exit 1
fi


# Build an index for alignment. Only needs to be run once.
STAR \
--runMode genomeGenerate \
--genomeDir $(pwd)/genome_data/index \
--genomeFastaFiles $(pwd)/$genomeFile \
--sjdbGTFfile $(pwd)/$annotationFile \
--runThreadN $NSLOTS

# Clean up temporary folder
rm -rf _STARtmp/

# Message completed
echo 'STAR Index generated!'
deactivate
