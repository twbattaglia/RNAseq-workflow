Differential Gene Expression using RNA-Seq (Workflow)
================

##### Thomas W. Battaglia (02/15/17)

### Introduction

RNAseq is becoming the one of the most prominent methods for measuring celluar responses. Not only does RNAseq have the ability to analyze differences in gene expression between samples, but can discover new isoforms and analyze SNP variations. This tutorial will cover the basic workflow for processing and analyzing differential gene expression data and is meant to give a general method for setting up an environment and running alignment tools. Be aware that is not meant to be used for all types of analyses and data-types, and the alignment tools are not for every analysis. Additionally, this tutorial is focused on giving a general sense of the flow when performing these analysis. For larger scale studies, it is highly reccomended to use a HPC environment for increased RAM and computational power.

### Getting Setup

#### 1. Installating Miniconda (if needed)

Miniconda is a comprehensive and easy to use package manager for Python (among other things). Miniconda is meant to replace your current Python installation with one that has more features and is modular, so you can delete it without any damage to your system. Not only does it allow you to install Python packages, you can create virtual environments and have access to large bioinformatics repositories (Bioconda).

``` bash
# Download the Miniconda3 installer to your home directory
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/minoconda.sh

# Run the miniconda installation
bash minoconda.sh -b -f -p ~/miniconda

# Add miniconda to the system path
echo 'PATH="$HOME/miniconda/bin:$PATH' >> ~/.bash_profile

# Source system file to activate miniconda
~/.bash_profile

# Add bioinformatic channels for downloading required packages
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda
```

#### 2. Setting Up the Folder Structure

Organizing is key to proper reproducible researcher. During the processing and analysis steps, many files are created. To best organize the analysis and increase the reprducibility of your analysis, it is best to use a simple folder structure. The struture allows other researchers and collaborators to find certain files or follow the steps used. The structure within this repository is just one way of organizing the data, but you can choose whichever way is the most comfortable.

``` bash
# Install git (if needed)
conda install --yes -c anaconda git wget

# Clone this repository with folder structure into the current working folder
git clone https://github.com/twbattaglia/RNAseq-workflow new_workflow
```

#### 3. Download Host Genome

To find either differentially expressed genes or isoform transcripts, you first need a reference genome to compare to. For any alignment, we need the host genome in `.fasta` format, but we also need an annotation file in `.GTF/.GFF`, which relates the coordinates in the genome to an annotated gene identifier. Both of these files are required to perform an alignment and gene abundance count. Be aware that the different resources (Ensembl, UCSC, RefSeq, Gencode) have different versions of the same species genome and annotation files cannot be mixed between versions. In this workflow, we will focus on the Gencode's genome. (<https://www.gencodegenes.org/>)

See here for a listing of genomes/annotation beyond mouse and human: <http://useast.ensembl.org/info/data/ftp/index.html>

##### Mouse (Gencode)

``` bash
# Download genome fasta file
wget -P genome/ ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/GRCm38.p5.genome.fa.gz

# Download annotation file
wget -P annotation/ ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.annotation.gtf.gz

# Decompress files for use with programs
gunzip genome/GRCm38.p4.genome.fa.gz
gunzip annotation/gencode.vM12.annotation.gtf.gz
```

##### Human (Gencode)

``` bash
# Download genome fasta file
wget -p genome/ ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh38.p7.genome.fa.gz

# Download annotation
wget -P annotation/ ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz

# Decompress files for use with programs
gunzip genome/GRCh38.p7.genome.fa.gz
gunzip annotation/gencode.v25.annotation.gtf.gz
```

------------------------------------------------------------------------

### Workflow

![RNAseq Workflow](README_files/rnaseq_workflow.jpg)

------------------------------------------------------------------------
