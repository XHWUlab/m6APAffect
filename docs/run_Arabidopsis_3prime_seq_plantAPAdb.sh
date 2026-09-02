#!/bin/bash -l

# ============================================================================
# Script name:
#   run_Arabidopsis_3prime_seq_plantAPAdb.sh
#
# Purpose:
#   Process Arabidopsis 3'-end sequencing data and generate PAC datasets
#   for downstream plantAPAdb and m6APAreg analyses.
#
# Supported data types:
#   1. A-seq data
#   2. PAS-seq data
#
# Main steps:
#   1. Remove adapters
#   2. Reverse-complement the A-seq R2 reads
#   3. Merge processed A-seq R1 and R2 reads
#   4. Build a STAR genome index
#   5. Align reads to the Arabidopsis genome
#   6. Convert SAM files into PAT, PA, and PAC files using plantAPAdb
#
# Required programs:
#   - cutadapt
#   - fastx_clipper
#   - fastx_reverse_complement
#   - STAR
#   - perl
#   - bedtools
#   - bedmap
#   - sort-bed
#   - awk
#   - sed
#   - plantAPAdb scripts:
#       MAP_parseSAM2PAT.pl
#       PAT_setIP_big.pl
#       PAT2PA2PAC.sh
#
# Usage:
#   bash run_Arabidopsis_3prime_seq_plantAPAdb.sh \
#       <project_dir> \
#       <genome_fasta> \
#       <annotation_gtf> \
#       <plantAPAdb_dir>
#
# Example:
#   bash run_Arabidopsis_3prime_seq_plantAPAdb.sh \
#       /path/to/Ara \
#       /path/to/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
#       /path/to/Arabidopsis_thaliana.TAIR10.42.gtf \
#       /path/to/plantAPAdb
#
# Important:
#   - Input FASTQ files must be placed in <project_dir>/raw_fastq.
#   - STAR output files are written to <project_dir>/star.
#   - plantAPAdb output files are written to:
#       <project_dir>/plantAPAdb/Aseq
#       <project_dir>/plantAPAdb/PASseq
#   - The chromosome names in the FASTA, GTF, and SAM files must be
#     mutually consistent.
# ============================================================================

set -euo pipefail

# ============================================================================
# 1. Read command-line arguments
# ============================================================================

if [[ "$#" -ne 4 ]]; then
    echo "Usage:"
    echo "  bash $0 <project_dir> <genome_fasta> <annotation_gtf> <plantAPAdb_dir>"
    exit 1
fi

project_dir="$1"
genome_fasta="$2"
annotation_gtf="$3"
plantAPAdb_dir="$4"

# ============================================================================
# 2. Configuration
# ============================================================================

# Number of threads used by adapter trimming.
THREADS_CUTADAPT=6

# Number of threads used to build the STAR index.
THREADS_STAR_INDEX=16

# Number of threads used for STAR alignment.
THREADS_STAR_ALIGN=20

# Distance used to merge nearby poly(A) sites into PACs.
PAC_DISTANCE=5

# Adapter sequences.
# A-seq reads contain the expected 5'-terminal sequence NNNNTTT.
ASEQ_ADAPTER="^NNNNTTT"

# PAS-seq reads contain the Illumina adapter sequence at the 3' end.
PASSEQ_ADAPTER="AGATCGGAAGAG"

# Input and output directories.
raw_fastq_dir="${project_dir}/raw_fastq"
trimmed_fastq_dir="${project_dir}/trimmed_fastq"
merged_fastq_dir="${project_dir}/merged_fastq"
star_dir="${project_dir}/star"
star_index_dir="${project_dir}/star_index"

plantAPAdb_aseq_dir="${project_dir}/plantAPAdb/Aseq"
plantAPAdb_passeq_dir="${project_dir}/plantAPAdb/PASseq"

# Create output directories.
mkdir -p \
    "${trimmed_fastq_dir}" \
    "${merged_fastq_dir}" \
    "${star_dir}" \
    "${star_index_dir}" \
    "${plantAPAdb_aseq_dir}" \
    "${plantAPAdb_passeq_dir}"

# ============================================================================
# 3. Check required files and programs
# ============================================================================

echo "Checking input files and required programs..."

required_commands=(
    cutadapt
    fastx_clipper
    fastx_reverse_complement
    STAR
    perl
    bedtools
    bedmap
    sort-bed
    awk
    sed
)

for command_name in "${required_commands[@]}"; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: ${command_name}"
        exit 1
    fi
done

if [[ ! -f "${genome_fasta}" ]]; then
    echo "ERROR: Genome FASTA file not found:"
    echo "  ${genome_fasta}"
    exit 1
fi

if [[ ! -f "${annotation_gtf}" ]]; then
    echo "ERROR: Annotation GTF file not found:"
    echo "  ${annotation_gtf}"
    exit 1
fi

if [[ ! -d "${plantAPAdb_dir}" ]]; then
    echo "ERROR: plantAPAdb directory not found:"
    echo "  ${plantAPAdb_dir}"
    exit 1
fi

if [[ ! -f "${plantAPAdb_dir}/MAP_parseSAM2PAT.pl" ]]; then
    echo "ERROR: MAP_parseSAM2PAT.pl not found in:"
    echo "  ${plantAPAdb_dir}"
    exit 1
fi

if [[ ! -f "${plantAPAdb_dir}/PAT_setIP_big.pl" ]]; then
    echo "ERROR: PAT_setIP_big.pl not found in:"
    echo "  ${plantAPAdb_dir}"
    exit 1
fi

if [[ ! -f "${plantAPAdb_dir}/PAT2PA2PAC.sh" ]]; then
    echo "ERROR: PAT2PA2PAC.sh not found in:"
    echo "  ${plantAPAdb_dir}"
    exit 1
fi

if [[ ! -d "${raw_fastq_dir}" ]]; then
    echo "ERROR: Raw FASTQ directory not found:"
    echo "  ${raw_fastq_dir}"
    exit 1
fi

echo "Input files and programs are available."

# ============================================================================
# 4. Build the STAR genome index
# ============================================================================

# The index is built only when the STAR index directory is empty.
# Remove this condition if the index needs to be rebuilt.

if [[ -z "$(find "${star_index_dir}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    echo "========================================================================"
    echo "Building STAR genome index"
    echo "========================================================================"

    STAR \
        --runMode genomeGenerate \
        --runThreadN "${THREADS_STAR_INDEX}" \
        --genomeFastaFiles "${genome_fasta}" \
        --sjdbGTFfile "${annotation_gtf}" \
        --sjdbGTFtagExonParentTranscript Parent \
        --genomeDir "${star_index_dir}"
else
    echo "========================================================================"
    echo "Existing STAR index detected"
    echo "Index directory: ${star_index_dir}"
    echo "========================================================================"
fi

# ============================================================================
# 5. Process A-seq data
# ============================================================================

echo "========================================================================"
echo "Processing A-seq data"
echo "========================================================================"

# ---------------------------------------------------------------------------
# A-seq sample definitions
#
# Sample names:
#   APA-seq-Colo-1
#   APA-seq-Colo-2
#   APA-seq-cpsf30-L-1
#   APA-seq-cpsf30-L-2
#
# The original files are:
#   CRR151243_f1.fq / CRR151243_r2.fq
#   CRR151244_f1.fq / CRR151244_r2.fq
#   CRR151245_f1.fq / CRR151245_r2.fq
#   CRR151246_f1.fq / CRR151246_r2.fq
# ---------------------------------------------------------------------------

aseq_ids=(
    CRR151243
    CRR151244
    CRR151245
    CRR151246
)

aseq_labels=(
    APA-seq-Colo-1
    APA-seq-Colo-2
    APA-seq-cpsf30-L-1
    APA-seq-cpsf30-L-2
)

for index in "${!aseq_ids[@]}"; do

    sample_id="${aseq_ids[$index]}"
    sample_label="${aseq_labels[$index]}"

    f1="${raw_fastq_dir}/${sample_id}_f1.fq"
    r2="${raw_fastq_dir}/${sample_id}_r2.fq"

    trimmed_f1="${trimmed_fastq_dir}/${sample_id}_f1.trimmed.fq"
    untrimmed_f1="${trimmed_fastq_dir}/${sample_id}_f1.untrimmed.fq"

    trimmed_r2="${trimmed_fastq_dir}/${sample_id}_r2.trimmed.fq"
    untrimmed_r2="${trimmed_fastq_dir}/${sample_id}_r2.untrimmed.fq"

    reverse_r2="${trimmed_fastq_dir}/${sample_id}_r2.reverse_complement.fq"

    merged_fastq="${merged_fastq_dir}/${sample_label}.fq"

    if [[ ! -f "${f1}" ]]; then
        echo "ERROR: A-seq forward-read file not found:"
        echo "  ${f1}"
        exit 1
    fi

    if [[ ! -f "${r2}" ]]; then
        echo "ERROR: A-seq reverse-read file not found:"
        echo "  ${r2}"
        exit 1
    fi

    echo "Processing A-seq sample: ${sample_label}"

    # Remove the 5'-terminal NNNNTTT sequence from the forward reads.
    # --no-trim preserves the sequence while requiring the adapter to be
    # present at the 5' end. Reads without the expected sequence are saved
    # separately.
    cutadapt \
        -g "${ASEQ_ADAPTER}" \
        -O 7 \
        -m 15 \
        -j "${THREADS_CUTADAPT}" \
        --no-trim \
        --untrimmed-output "${untrimmed_f1}" \
        -o "${trimmed_f1}" \
        "${f1}"

    # Reverse-complement the A-seq R2 reads before adapter processing.
    fastx_reverse_complement \
        -i "${r2}" \
        -o "${reverse_r2}"

    # Remove the 5'-terminal NNNNTTT sequence from the reverse-complemented
    # R2 reads.
    cutadapt \
        -g "${ASEQ_ADAPTER}" \
        -O 7 \
        -m 15 \
        -j "${THREADS_CUTADAPT}" \
        --no-trim \
        --untrimmed-output "${untrimmed_r2}" \
        -o "${trimmed_r2}" \
        "${reverse_r2}"

    # Combine the processed forward and reverse-complemented reads.
    # plantAPAdb will treat the resulting file as one PAT-generating input.
    cat \
        "${trimmed_f1}" \
        "${trimmed_r2}" \
        > "${merged_fastq}"

    echo "Generated merged A-seq FASTQ:"
    echo "  ${merged_fastq}"
done

# ============================================================================
# 6. Align A-seq reads with STAR
# ============================================================================

echo "========================================================================"
echo "Aligning A-seq reads with STAR"
echo "========================================================================"

aseq_sam_files=()

for index in "${!aseq_labels[@]}"; do

    sample_label="${aseq_labels[$index]}"
    merged_fastq="${merged_fastq_dir}/${sample_label}.fq"
    star_prefix="${star_dir}/Aseq_${sample_label}_"

    echo "Aligning A-seq sample: ${sample_label}"

    STAR \
        --runThreadN "${THREADS_STAR_ALIGN}" \
        --genomeDir "${star_index_dir}" \
        --readFilesIn "${merged_fastq}" \
        --outFileNamePrefix "${star_prefix}" \
        --outMultimapperOrder Random \
        --outFilterMultimapNmax 1 \
        --outSAMtype SAM

    sam_file="${star_prefix}Aligned.out.sam"

    if [[ ! -f "${sam_file}" ]]; then
        echo "ERROR: STAR did not generate the expected SAM file:"
        echo "  ${sam_file}"
        exit 1
    fi

    aseq_sam_files+=("${sam_file}")
done

# ============================================================================
# 7. Convert A-seq alignments into plantAPAdb PAC files
# ============================================================================

echo "========================================================================"
echo "Running plantAPAdb for A-seq data"
echo "========================================================================"

# PAT2PA2PAC.sh expects its SAM input list as one comma-separated argument.
aseq_sam_list="$(IFS=,; echo "${aseq_sam_files[*]}")"

(
    cd "${plantAPAdb_dir}"

    bash PAT2PA2PAC.sh \
        "${genome_fasta}" \
        "${star_dir}" \
        "${plantAPAdb_aseq_dir}" \
        "${PAC_DISTANCE}" \
        "${aseq_sam_list}"
)

# Rename the generated header labels to match the A-seq sample names.
if [[ -f "${plantAPAdb_aseq_dir}/all.PAC.header" ]]; then

    {
        printf "chr\tUPA_start\tUPA_end\tstrand\tPAnum\ttot_tagnum\tcoord\trefPAnum"

        for sample_label in "${aseq_labels[@]}"; do
            printf "\t%s" "${sample_label}"
        done

        printf "\n"
    } > "${plantAPAdb_aseq_dir}/all.PAC.header"
fi

# ============================================================================
# 8. Process PAS-seq data
# ============================================================================

echo "========================================================================"
echo "Processing PAS-seq data"
echo "========================================================================"

# ---------------------------------------------------------------------------
# PAS-seq sample definitions
#
# Control samples:
#   SRR11068176
#   SRR11068177
#   SRR11068178
#
# cpsf30-l samples:
#   SRR11068179
#   SRR11068180
#   SRR11068181
# ---------------------------------------------------------------------------

passeq_ids=(
    SRR11068176
    SRR11068177
    SRR11068178
    SRR11068179
    SRR11068180
    SRR11068181
)

passeq_labels=(
    col-rep1
    col-rep2
    col-rep3
    cpsf30-l-rep1
    cpsf30-l-rep2
    cpsf30-l-rep3
)

for index in "${!passeq_ids[@]}"; do

    sample_id="${passeq_ids[$index]}"
    sample_label="${passeq_labels[$index]}"

    input_fastq="${raw_fastq_dir}/${sample_id}.fastq"
    output_fastq="${trimmed_fastq_dir}/PASseq_${sample_label}.fastq"

    if [[ ! -f "${input_fastq}" ]]; then
        echo "ERROR: PAS-seq FASTQ file not found:"
        echo "  ${input_fastq}"
        exit 1
    fi

    echo "Processing PAS-seq sample: ${sample_label}"

    # Remove the 3' Illumina adapter and discard reads shorter than 20 nt.
    fastx_clipper \
        -a "${PASSEQ_ADAPTER}" \
        -l 20 \
        -i "${input_fastq}" \
        -o "${output_fastq}"
done

# ============================================================================
# 9. Align PAS-seq reads with STAR
# ============================================================================

echo "========================================================================"
echo "Aligning PAS-seq reads with STAR"
echo "========================================================================"

passeq_sam_files=()

for index in "${!passeq_ids[@]}"; do

    sample_id="${passeq_ids[$index]}"
    sample_label="${passeq_labels[$index]}"

    input_fastq="${trimmed_fastq_dir}/PASseq_${sample_label}.fastq"
    star_prefix="${star_dir}/PASseq_${sample_label}_"

    echo "Aligning PAS-seq sample: ${sample_label}"

    STAR \
        --runThreadN "${THREADS_STAR_ALIGN}" \
        --genomeDir "${star_index_dir}" \
        --readFilesIn "${input_fastq}" \
        --outFileNamePrefix "${star_prefix}" \
        --outMultimapperOrder Random \
        --outFilterMultimapNmax 1 \
        --outSAMtype SAM

    sam_file="${star_prefix}Aligned.out.sam"

    if [[ ! -f "${sam_file}" ]]; then
        echo "ERROR: STAR did not generate the expected SAM file:"
        echo "  ${sam_file}"
        exit 1
    fi

    passeq_sam_files+=("${sam_file}")
done

# ============================================================================
# 10. Convert PAS-seq alignments into plantAPAdb PAC files
# ============================================================================

echo "========================================================================"
echo "Running plantAPAdb for PAS-seq data"
echo "========================================================================"

passeq_sam_list="$(IFS=,; echo "${passeq_sam_files[*]}")"

(
    cd "${plantAPAdb_dir}"

    bash PAT2PA2PAC.sh \
        "${genome_fasta}" \
        "${star_dir}" \
        "${plantAPAdb_passeq_dir}" \
        "${PAC_DISTANCE}" \
        "${passeq_sam_list}"
)

# Rename the generated header labels to match the PAS-seq sample names.
if [[ -f "${plantAPAdb_passeq_dir}/all.PAC.header" ]]; then

    {
        printf "chr\tUPA_start\tUPA_end\tstrand\tPAnum\ttot_tagnum\tcoord\trefPAnum"

        for sample_label in "${passeq_labels[@]}"; do
            printf "\t%s" "${sample_label}"
        done

        printf "\n"
    } > "${plantAPAdb_passeq_dir}/all.PAC.header"
fi

# ============================================================================
# 11. Final output summary
# ============================================================================

echo
echo "========================================================================"
echo "Pipeline completed successfully"
echo "========================================================================"

echo "A-seq outputs:"
echo "  ${plantAPAdb_aseq_dir}/all.PAC.info"
echo "  ${plantAPAdb_aseq_dir}/all.PAC.header"
echo "  ${plantAPAdb_aseq_dir}/all.PAC.PAcount"
echo "  ${plantAPAdb_aseq_dir}/all.PAC.PATcount"

echo
echo "PAS-seq outputs:"
echo "  ${plantAPAdb_passeq_dir}/all.PAC.info"
echo "  ${plantAPAdb_passeq_dir}/all.PAC.header"
echo "  ${plantAPAdb_passeq_dir}/all.PAC.PAcount"
echo "  ${plantAPAdb_passeq_dir}/all.PAC.PATcount"

echo
echo "STAR alignment files:"
echo "  ${star_dir}"
