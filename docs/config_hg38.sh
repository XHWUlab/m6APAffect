#!/usr/bin/env bash

# ============================================================================
# Project configuration
# ============================================================================

PROJECT_DIR="/path/to/project"

# One accession per line.
ACCESSION_LIST="${PROJECT_DIR}/SRR_Acc_List.txt"

# Tab-separated sample metadata file.
SAMPLE_METADATA="${PROJECT_DIR}/sample_metadata.tsv"

# Human reference genome and gene annotation.
GENOME_FASTA="/path/to/Homo_sapiens.GRCh38.dna.toplevel.fa"
GTF_FILE="/path/to/Homo_sapiens.GRCh38.113.gtf"

# rRNA reference sequence.
RRNA_FASTA="/path/to/rRNA_hg38.fasta"

# QAPA annotation files.
POLYASITE_BED="/path/to/clusters.hg38.bed"
GENCODE_BASIC="/path/to/gencode.basic.txt"
ENSEMBL_IDENTIFIERS="/path/to/ensembl_identifiers.txt"

# Sequencing mode:
#   PE = paired-end sequencing
#   SE = single-end sequencing
SEQUENCING_MODE="PE"

# Thread configuration.
THREADS_DOWNLOAD=8
THREADS_FASTQ=8
THREADS_FASTP=8
THREADS_FASTQC=12
THREADS_CUTADAPT=6
THREADS_BOWTIE2=8
THREADS_HISAT2=8
THREADS_SALMON=8

# HISAT2 multimapping parameter.
HISAT2_K=1

