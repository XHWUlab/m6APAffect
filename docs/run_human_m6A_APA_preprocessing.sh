#!/usr/bin/env bash

# ============================================================================
# Script name:
#   run_human_m6A_APA_preprocessing.sh
#
# Purpose:
#   Process human cancer and matched normal RNA-seq-related data for:
#
#   1. m6A peak identification
#   2. APA quantification using QAPA and Salmon
#   3. Generation of sorted BAM files for m6Aexpress
#   4. Generation of QAPAresults.txt files for downstream m6APAreg analysis
#
# The script is designed for one project or one cancer type at a time.
#
# Required software:
#   - SRA Toolkit: prefetch, fasterq-dump
#   - fastp
#   - FastQC
#   - cutadapt
#   - Bowtie2
#   - HISAT2
#   - samtools
#   - QAPA
#   - Salmon
#   - m6Aexpress in R
#
# Required input files:
#   1. SRR_Acc_List.txt
#   2. sample_metadata.tsv
#   3. Human genome FASTA
#   4. Human GTF annotation
#   5. rRNA reference FASTA
#   6. PolyASite BED file
#   7. GENCODE gene annotation table
#   8. Ensembl identifier table
#
# Usage:
#   bash run_human_m6A_APA_preprocessing.sh config.sh
#
# Example:
#   bash run_human_m6A_APA_preprocessing.sh config_hg38.sh
#
# Notes:
#   - The default workflow assumes paired-end sequencing.
#   - For single-end data, set SEQUENCING_MODE="SE" in config.sh.
#   - Each sample must have a unique sample ID.
#   - Sample grouping for m6Aexpress is controlled by sample_metadata.tsv.
# ============================================================================

set -Eeuo pipefail

# ============================================================================
# 1. Read configuration
# ============================================================================

if [[ "$#" -ne 1 ]]; then
    echo "Usage:"
    echo "  bash $0 <configuration_file>"
    exit 1
fi

config_file="$1"

if [[ ! -f "${config_file}" ]]; then
    echo "ERROR: Configuration file not found:"
    echo "  ${config_file}"
    exit 1
fi

# shellcheck source=/dev/null
source "${config_file}"

# ============================================================================
# 2. Check mandatory configuration variables
# ============================================================================

required_variables=(
    PROJECT_DIR
    ACCESSION_LIST
    SAMPLE_METADATA
    GENOME_FASTA
    GTF_FILE
    RRNA_FASTA
    POLYASITE_BED
    GENCODE_BASIC
    ENSEMBL_IDENTIFIERS
    SEQUENCING_MODE
)

for variable_name in "${required_variables[@]}"; do
    if [[ -z "${!variable_name:-}" ]]; then
        echo "ERROR: Required configuration variable is empty:"
        echo "  ${variable_name}"
        exit 1
    fi
done

# ============================================================================
# 3. Software parameters
# ============================================================================

THREADS_DOWNLOAD="${THREADS_DOWNLOAD:-8}"
THREADS_FASTQ="${THREADS_FASTQ:-8}"
THREADS_FASTP="${THREADS_FASTP:-8}"
THREADS_FASTQC="${THREADS_FASTQC:-12}"
THREADS_BOWTIE2="${THREADS_BOWTIE2:-8}"
THREADS_HISAT2="${THREADS_HISAT2:-8}"
THREADS_SALMON="${THREADS_SALMON:-8}"

# Maximum number of multimapping locations allowed by HISAT2.
HISAT2_K="${HISAT2_K:-1}"

# ============================================================================
# 4. Directory layout
# ============================================================================

RAW_SRA_DIR="${PROJECT_DIR}/raw_sra"
RAW_FASTQ_DIR="${PROJECT_DIR}/raw_fastq"
FASTP_DIR="${PROJECT_DIR}/fastp"
FASTQC_DIR="${PROJECT_DIR}/fastqc"
CUTADAPT_DIR="${PROJECT_DIR}/cutadapt"
RRNA_DIR="${PROJECT_DIR}/rm_rRNA_fastq"
ALIGNMENT_DIR="${PROJECT_DIR}/alignment"
BAM_DIR="${PROJECT_DIR}/bam"
M6A_DIR="${PROJECT_DIR}/m6A_results"

QAPA_DIR="${PROJECT_DIR}/QAPA"
QAPA_LIBRARY_DIR="${QAPA_DIR}/library"
SALMON_DIR="${QAPA_DIR}/salmon"
QAPA_PROJECT_DIR="${QAPA_DIR}/project"

INDEX_DIR="${PROJECT_DIR}/index"
RRNA_INDEX_DIR="${INDEX_DIR}/rRNA_hg38"
HISAT2_INDEX_DIR="${INDEX_DIR}/hg38"
SALMON_INDEX_DIR="${QAPA_LIBRARY_DIR}/utr_library"

mkdir -p \
    "${RAW_SRA_DIR}" \
    "${RAW_FASTQ_DIR}" \
    "${FASTP_DIR}" \
    "${FASTQC_DIR}" \
    "${CUTADAPT_DIR}" \
    "${RRNA_DIR}" \
    "${ALIGNMENT_DIR}" \
    "${BAM_DIR}" \
    "${M6A_DIR}" \
    "${QAPA_LIBRARY_DIR}" \
    "${SALMON_DIR}" \
    "${QAPA_PROJECT_DIR}" \
    "${INDEX_DIR}" \
    "${RRNA_INDEX_DIR}" \
    "${HISAT2_INDEX_DIR}"

LOG_FILE="${PROJECT_DIR}/human_m6A_APA_preprocessing.log"

exec > >(tee -a "${LOG_FILE}") 2>&1

echo "========================================================================"
echo "Human m6A and APA preprocessing pipeline"
echo "========================================================================"
echo "Project directory: ${PROJECT_DIR}"
echo "Sequencing mode: ${SEQUENCING_MODE}"
echo "Started at: $(date)"
echo

# ============================================================================
# 5. Check required programs
# ============================================================================

echo "Checking required programs..."

required_commands=(
    prefetch
    fasterq-dump
    fastp
    fastqc
    cutadapt
    bowtie2
    bowtie2-build
    hisat2
    hisat2-build
    samtools
    qapa
    salmon
    Rscript
)

for command_name in "${required_commands[@]}"; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
        echo "ERROR: Required program was not found:"
        echo "  ${command_name}"
        exit 1
    fi
done

echo "All required programs are available."

# ============================================================================
# 6. Validate input files
# ============================================================================

echo "Checking input files..."

input_files=(
    "${ACCESSION_LIST}"
    "${SAMPLE_METADATA}"
    "${GENOME_FASTA}"
    "${GTF_FILE}"
    "${RRNA_FASTA}"
    "${POLYASITE_BED}"
    "${GENCODE_BASIC}"
    "${ENSEMBL_IDENTIFIERS}"
)

for input_file in "${input_files[@]}"; do
    if [[ ! -f "${input_file}" ]]; then
        echo "ERROR: Input file not found:"
        echo "  ${input_file}"
        exit 1
    fi
done

echo "All input files are available."

# ============================================================================
# 7. Download SRA files
# ============================================================================

echo "========================================================================"
echo "Downloading SRA files"
echo "========================================================================"

# Each line of ACCESSION_LIST should contain one accession, for example:
#
# SRR19688220
# SRR19688221
# SRR19688222
#
# Empty lines and lines beginning with '#' are ignored.

grep -vE '^[[:space:]]*($|#)' "${ACCESSION_LIST}" \
    | while read -r accession; do

        echo "Downloading: ${accession}"

        prefetch \
            --max-size 100G \
            --output-directory "${RAW_SRA_DIR}" \
            "${accession}"

    done

# ============================================================================
# 8. Convert SRA files to FASTQ
# ============================================================================

echo "========================================================================"
echo "Converting SRA files to FASTQ"
echo "========================================================================"

if [[ "${SEQUENCING_MODE}" == "PE" ]]; then

    for sra_file in "${RAW_SRA_DIR}"/*.sra; do

        if [[ ! -f "${sra_file}" ]]; then
            continue
        fi

        accession="$(basename "${sra_file}" .sra)"

        echo "Converting paired-end SRA: ${accession}"

        fasterq-dump \
            --split-files \
            --threads "${THREADS_FASTQ}" \
            --outdir "${RAW_FASTQ_DIR}" \
            "${sra_file}"

    done

elif [[ "${SEQUENCING_MODE}" == "SE" ]]; then

    for sra_file in "${RAW_SRA_DIR}"/*.sra; do

        if [[ ! -f "${sra_file}" ]]; then
            continue
        fi

        accession="$(basename "${sra_file}" .sra)"

        echo "Converting single-end SRA: ${accession}"

        fasterq-dump \
            --threads "${THREADS_FASTQ}" \
            --outdir "${RAW_FASTQ_DIR}" \
            "${sra_file}"

    done

else
    echo "ERROR: SEQUENCING_MODE must be PE or SE."
    exit 1
fi

# ============================================================================
# 9. FASTP quality control and adapter removal
# ============================================================================

echo "========================================================================"
echo "Running fastp"
echo "========================================================================"

if [[ "${SEQUENCING_MODE}" == "PE" ]]; then

    for r1_file in "${RAW_FASTQ_DIR}"/*_1.fastq; do

        if [[ ! -f "${r1_file}" ]]; then
            continue
        fi

        accession="$(basename "${r1_file}" _1.fastq)"
        r2_file="${RAW_FASTQ_DIR}/${accession}_2.fastq"

        if [[ ! -f "${r2_file}" ]]; then
            echo "ERROR: Matching R2 file not found for ${accession}:"
            echo "  ${r2_file}"
            exit 1
        fi

        fastp \
            --in1 "${r1_file}" \
            --in2 "${r2_file}" \
            --out1 "${FASTP_DIR}/${accession}_1.fastq" \
            --out2 "${FASTP_DIR}/${accession}_2.fastq" \
            --thread "${THREADS_FASTP}" \
            --html "${FASTP_DIR}/${accession}.fastp.html" \
            --json "${FASTP_DIR}/${accession}.fastp.json"
    done

else

    for input_file in "${RAW_FASTQ_DIR}"/*.fastq; do

        if [[ ! -f "${input_file}" ]]; then
            continue
        fi

        accession="$(basename "${input_file}" .fastq)"

        fastp \
            --in1 "${input_file}" \
            --out1 "${FASTP_DIR}/${accession}.fastq" \
            --thread "${THREADS_FASTP}" \
            --html "${FASTP_DIR}/${accession}.fastp.html" \
            --json "${FASTP_DIR}/${accession}.fastp.json"
    done

fi

# ============================================================================
# 10. Run FastQC
# ============================================================================

echo "========================================================================"
echo "Running FastQC"
echo "========================================================================"

fastqc \
    --threads "${THREADS_FASTQC}" \
    --outdir "${FASTQC_DIR}" \
    "${FASTP_DIR}"/*.fastq

# ============================================================================
# 11. Optional additional poly(A) trimming with cutadapt
# ============================================================================

echo "========================================================================"
echo "Running additional poly(A) trimming"
echo "========================================================================"

# fastp is used as the primary quality-control and adapter-removal step.
# cutadapt is applied here to remove remaining long poly(A) sequences.
#
# The adapter sequence A{10} is treated as a 3' adapter.
# Reads are retained in their original orientation.

if [[ "${SEQUENCING_MODE}" == "PE" ]]; then

    for r1_file in "${FASTP_DIR}"/*_1.fastq; do

        if [[ ! -f "${r1_file}" ]]; then
            continue
        fi

        accession="$(basename "${r1_file}" _1.fastq)"
        r2_file="${FASTP_DIR}/${accession}_2.fastq"

        cutadapt \
            -a "A{10}" \
            -j "${THREADS_CUTADAPT:-6}" \
            -o "${CUTADAPT_DIR}/${accession}_1.fastq" \
            "${r1_file}"

        cutadapt \
            -a "A{10}" \
            -j "${THREADS_CUTADAPT:-6}" \
            -o "${CUTADAPT_DIR}/${accession}_2.fastq" \
            "${r2_file}"
    done

else

    for input_file in "${FASTP_DIR}"/*.fastq; do

        if [[ ! -f "${input_file}" ]]; then
            continue
        fi

        accession="$(basename "${input_file}" .fastq)"

        cutadapt \
            -a "A{10}" \
            -j "${THREADS_CUTADAPT:-6}" \
            -o "${CUTADAPT_DIR}/${accession}.fastq" \
            "${input_file}"
    done

fi

# ============================================================================
# 12. Build Bowtie2 rRNA index
# ============================================================================

echo "========================================================================"
echo "Building Bowtie2 rRNA index"
echo "========================================================================"

if [[ ! -f "${RRNA_INDEX_DIR}.1.bt2" &&
      ! -f "${RRNA_INDEX_DIR}.1.bt2l" ]]; then

    bowtie2-build \
        --threads "${THREADS_BOWTIE2}" \
        "${RRNA_FASTA}" \
        "${RRNA_INDEX_DIR}"
else
    echo "Existing Bowtie2 rRNA index detected."
fi

# ============================================================================
# 13. Remove rRNA reads
# ============================================================================

echo "========================================================================"
echo "Removing rRNA reads"
echo "========================================================================"

if [[ "${SEQUENCING_MODE}" == "PE" ]]; then

    for r1_file in "${CUTADAPT_DIR}"/*_1.fastq; do

        if [[ ! -f "${r1_file}" ]]; then
            continue
        fi

        accession="$(basename "${r1_file}" _1.fastq)"
        r2_file="${CUTADAPT_DIR}/${accession}_2.fastq"

        bowtie2 \
            -x "${RRNA_INDEX_DIR}" \
            --un-conc-gz "${RRNA_DIR}/${accession}_rmrRNA.fastq.gz" \
            -1 "${r1_file}" \
            -2 "${r2_file}" \
            -p "${THREADS_BOWTIE2}" \
            -S "${RRNA_DIR}/${accession}_rRNA.sam"

        rm -f "${RRNA_DIR}/${accession}_rRNA.sam"
    done

else

    for input_file in "${CUTADAPT_DIR}"/*.fastq; do

        if [[ ! -f "${input_file}" ]]; then
            continue
        fi

        accession="$(basename "${input_file}" .fastq)"

        bowtie2 \
            -x "${RRNA_INDEX_DIR}" \
            --un-gz "${RRNA_DIR}/${accession}_rmrRNA.fastq.gz" \
            -U "${input_file}" \
            -p "${THREADS_BOWTIE2}" \
            -S "${RRNA_DIR}/${accession}_rRNA.sam"

        rm -f "${RRNA_DIR}/${accession}_rRNA.sam"
    done

fi

# ============================================================================
# 14. Build HISAT2 genome index
# ============================================================================

echo "========================================================================"
echo "Building HISAT2 genome index"
echo "========================================================================"

if [[ ! -f "${HISAT2_INDEX_DIR}.1.ht2" &&
      ! -f "${HISAT2_INDEX_DIR}.1.ht2l" ]]; then

    hisat2-build \
        -p "${THREADS_HISAT2}" \
        "${GENOME_FASTA}" \
        "${HISAT2_INDEX_DIR}"
else
    echo "Existing HISAT2 genome index detected."
fi

# ============================================================================
# 15. Align rRNA-depleted reads to the human genome
# ============================================================================

echo "========================================================================"
echo "Aligning rRNA-depleted reads to the human genome"
echo "========================================================================"

if [[ "${SEQUENCING_MODE}" == "PE" ]]; then

    for r1_file in "${RRNA_DIR}"/*_rmrRNA.fastq.1.gz; do

        if [[ ! -f "${r1_file}" ]]; then
            continue
        fi

        basename_r1="$(basename "${r1_file}")"
        accession="${basename_r1%_rmrRNA.fastq.1.gz}"

        r2_file="${RRNA_DIR}/${accession}_rmrRNA.fastq.2.gz"
        sam_file="${ALIGNMENT_DIR}/${accession}.sam"
        bam_file="${BAM_DIR}/${accession}.bam"
        sorted_bam="${BAM_DIR}/${accession}.sorted.bam"

        hisat2 \
            -p "${THREADS_HISAT2}" \
            -x "${HISAT2_INDEX_DIR}" \
            -1 "${r1_file}" \
            -2 "${r2_file}" \
            -k "${HISAT2_K}" \
            -S "${sam_file}"

        samtools view \
            -@ "${THREADS_HISAT2}" \
            -b \
            "${sam_file}" \
            -o "${bam_file}"

        samtools sort \
            -@ "${THREADS_HISAT2}" \
            "${bam_file}" \
            -o "${sorted_bam}"

        samtools index \
            -@ "${THREADS_HISAT2}" \
            "${sorted_bam}"

        rm -f "${sam_file}" "${bam_file}"
    done

else

    for input_file in "${RRNA_DIR}"/*_rmrRNA.fastq.gz; do

        if [[ ! -f "${input_file}" ]]; then
            continue
        fi

        basename_input="$(basename "${input_file}")"
        accession="${basename_input%_rmrRNA.fastq.gz}"

        sam_file="${ALIGNMENT_DIR}/${accession}.sam"
        bam_file="${BAM_DIR}/${accession}.bam"
        sorted_bam="${BAM_DIR}/${accession}.sorted.bam"

        hisat2 \
            -p "${THREADS_HISAT2}" \
            -x "${HISAT2_INDEX_DIR}" \
            -U "${input_file}" \
            -k "${HISAT2_K}" \
            -S "${sam_file}"

        samtools view \
            -@ "${THREADS_HISAT2}" \
            -b \
            "${sam_file}" \
            -o "${bam_file}"

        samtools sort \
            -@ "${THREADS_HISAT2}" \
            "${bam_file}" \
            -o "${sorted_bam}"

        samtools index \
            -@ "${THREADS_HISAT2}" \
            "${sorted_bam}"

        rm -f "${sam_file}" "${bam_file}"
    done

fi

# ============================================================================
# 16. Generate an m6Aexpress sample manifest
# ============================================================================

echo "========================================================================"
echo "Generating m6Aexpress sample manifest"
echo "========================================================================"

# SAMPLE_METADATA must contain the following tab-separated columns:
#
# sample_id    cancer_type    condition    assay    replicate
#
# Example:
#
# sample_id       cancer_type   condition   assay   replicate
# SRR19688220    glioma        normal      IP      1
# SRR19688221    glioma        normal      INPUT   1
# SRR19688222    glioma        cancer      IP      1
# SRR19688223    glioma        cancer      INPUT   1
#
# The manifest is used to construct the R input vectors for m6Aexpress.

awk -F '\t' '
BEGIN {
    OFS = "\t"
}
!/^#/ && NR > 1 && NF >= 5 {
    print $1, $2, $3, $4, $5
}
' "${SAMPLE_METADATA}" \
    > "${M6A_DIR}/m6Aexpress_sample_manifest.tsv"

# ============================================================================
# 17. Create the m6Aexpress R script
# ============================================================================

echo "========================================================================"
echo "Creating m6Aexpress analysis script"
echo "========================================================================"

m6a_r_script="${M6A_DIR}/run_m6Aexpress.R"

cat > "${m6a_r_script}" <<EOF
#!/usr/bin/env Rscript

# ============================================================================
# m6Aexpress peak-calling script
# ============================================================================

library(m6Aexpress)

project_dir <- "${PROJECT_DIR}"
gtf_file <- "${GTF_FILE}"
metadata_file <- "${M6A_DIR}/m6Aexpress_sample_manifest.tsv"
output_dir <- "${M6A_DIR}"

metadata <- read.delim(
    metadata_file,
    header = FALSE,
    sep = "\\t",
    stringsAsFactors = FALSE,
    col.names = c(
        "sample_id",
        "cancer_type",
        "condition",
        "assay",
        "replicate"
    )
)

# Extract the cancer type represented in this project.
cancer_types <- unique(metadata\$cancer_type)

if (length(cancer_types) != 1L) {
    stop(
        "The m6Aexpress input must contain exactly one cancer type. ",
        "Detected: ",
        paste(cancer_types, collapse = ", ")
    )
}

cancer_type <- cancer_types[[1]]

# Construct sorted BAM paths.
metadata\$bam_file <- file.path(
    "${BAM_DIR}",
    paste0(metadata\$sample_id, ".sorted.bam")
)

if (any(!file.exists(metadata\$bam_file))) {
    missing_files <- metadata\$bam_file[!file.exists(metadata\$bam_file)]
    stop(
        "The following sorted BAM files are missing:\\n",
        paste(missing_files, collapse = "\\n")
    )
}

# Define BAM files for the four required groups.
normal_ip <- metadata\$bam_file[
    metadata\$condition == "normal" &
        metadata\$assay == "IP"
]

normal_input <- metadata\$bam_file[
    metadata\$condition == "normal" &
        metadata\$assay == "INPUT"
]

cancer_ip <- metadata\$bam_file[
    metadata\$condition == "cancer" &
        metadata\$assay == "IP"
]

cancer_input <- metadata\$bam_file[
    metadata\$condition == "cancer" &
        metadata\$assay == "INPUT"
]

if (length(normal_ip) == 0L ||
    length(normal_input) == 0L ||
    length(cancer_ip) == 0L ||
    length(cancer_input) == 0L) {
    stop(
        "Each project must contain normal IP, normal INPUT, ",
        "cancer IP, and cancer INPUT BAM files."
    )
}

# Identify differential methylation peaks between the normal and cancer
# conditions.
peak_information <- Get_peakinfor(
    normal_ip,
    normal_input,
    cancer_ip,
    cancer_input,
    GENE_ANNO_GTF = gtf_file,
    UCSC_TABLE_NAME = NA
)

saveRDS(
    peak_information,
    file = file.path(
        output_dir,
        paste0(cancer_type, "_peak_information.rds")
    )
)

# Call enriched peaks in cancer samples.
cancer_peak_results <- exomepeak(
    GENE_ANNO_GTF = gtf_file,
    IP_BAM = cancer_ip,
    INPUT_BAM = cancer_input,
    OUTPUT_DIR = file.path(
        output_dir,
        paste0(cancer_type, "_cancer_peaks")
    )
)

# Call enriched peaks in normal samples.
normal_peak_results <- exomepeak(
    GENE_ANNO_GTF = gtf_file,
    IP_BAM = normal_ip,
    INPUT_BAM = normal_input,
    OUTPUT_DIR = file.path(
        output_dir,
        paste0(cancer_type, "_normal_peaks")
    )
)

saveRDS(
    cancer_peak_results,
    file = file.path(
        output_dir,
        paste0(cancer_type, "_cancer_peak_results.rds")
    )
)

saveRDS(
    normal_peak_results,
    file = file.path(
        output_dir,
        paste0(cancer_type, "_normal_peak_results.rds")
    )
)
EOF

# ============================================================================
# 18. Build the QAPA 3' UTR library
# ============================================================================

echo "========================================================================"
echo "Building QAPA 3' UTR library"
echo "========================================================================"

QAPA_UTR_BED="${QAPA_LIBRARY_DIR}/output_utrs.bed"
QAPA_UTR_BED_MODIFIED="${QAPA_LIBRARY_DIR}/output_utrs_modified.bed"
QAPA_UTR_FASTA="${QAPA_LIBRARY_DIR}/output_sequences.fa"

if [[ ! -f "${QAPA_UTR_BED_MODIFIED}" ]]; then

    qapa build \
        --db "${ENSEMBL_IDENTIFIERS}" \
        -g "${POLYASITE_BED}" \
        -p "${GENCODE_BASIC}" \
        > "${QAPA_UTR_BED}"

    # QAPA output can contain UCSC-style chromosome names beginning with
    # "chr", whereas the Ensembl FASTA may use names without this prefix.
    # Modify this command only if the FASTA chromosome names do not contain
    # the "chr" prefix.
    sed 's/^chr//' \
        "${QAPA_UTR_BED}" \
        > "${QAPA_UTR_BED_MODIFIED}"
else
    echo "Existing QAPA UTR BED file detected."
fi

if [[ ! -f "${QAPA_UTR_FASTA}" ]]; then

    qapa fasta \
        -f "${GENOME_FASTA}" \
        "${QAPA_UTR_BED_MODIFIED}" \
        "${QAPA_UTR_FASTA}"
else
    echo "Existing QAPA UTR FASTA detected."
fi

if [[ ! -f "${SALMON_INDEX_DIR}/versionInfo.json" ]]; then

    salmon index \
        -t "${QAPA_UTR_FASTA}" \
        -i "${SALMON_INDEX_DIR}"
else
    echo "Existing Salmon index detected."
fi

# ============================================================================
# 19. Quantify APA usage with Salmon
# ============================================================================

echo "========================================================================"
echo "Quantifying APA usage with Salmon"
echo "========================================================================"

if [[ "${SEQUENCING_MODE}" == "PE" ]]; then

    for r1_file in "${CUTADAPT_DIR}"/*_1.fastq; do

        if [[ ! -f "${r1_file}" ]]; then
            continue
        fi

        accession="$(basename "${r1_file}" _1.fastq)"
        r2_file="${CUTADAPT_DIR}/${accession}_2.fastq"

        salmon quant \
            -l A \
            -i "${SALMON_INDEX_DIR}" \
            --validateMappings \
            -p "${THREADS_SALMON}" \
            -1 "${r1_file}" \
            -2 "${r2_file}" \
            -o "${SALMON_DIR}/salmon_${accession}"
    done

else

    for input_file in "${CUTADAPT_DIR}"/*.fastq; do

        if [[ ! -f "${input_file}" ]]; then
            continue
        fi

        accession="$(basename "${input_file}" .fastq)"

        salmon quant \
            -l ISR \
            -i "${SALMON_INDEX_DIR}" \
            --validateMappings \
            -p "${THREADS_SALMON}" \
            -r "${input_file}" \
            -o "${SALMON_DIR}/salmon_${accession}"
    done

fi

# ============================================================================
# 20. Create QAPA sample directories
# ============================================================================

echo "========================================================================"
echo "Preparing QAPA sample directories"
echo "========================================================================"

# sample_metadata.tsv must contain a column called "sample_index".
#
# Example:
#
# sample_index    sample_id       cancer_type   condition   assay   replicate
# sample1         SRR19688220    glioma        normal      INPUT   1
# sample2         SRR19688221    glioma        normal      IP      1
# sample3         SRR19688222    glioma        cancer      INPUT   1
# sample4         SRR19688223    glioma        cancer      IP      1
#
# Each sample_index corresponds to one QAPA sample directory.

awk -F '\t' '
BEGIN {
    OFS = "\t"
}
!/^#/ && NR > 1 && NF >= 1 {
    print $1
}
' "${SAMPLE_METADATA}" \
    | while read -r sample_index; do

        sample_id="$(
            awk -F '\t' -v sample_index="${sample_index}" '
            NR > 1 && $1 == sample_index {
                print $2
                exit
            }
            ' "${SAMPLE_METADATA}"
        )"

        if [[ -z "${sample_id}" ]]; then
            echo "ERROR: Cannot find sample ID for ${sample_index}"
            exit 1
        fi

        mkdir -p "${QAPA_PROJECT_DIR}/${sample_index}"

        quant_file="${SALMON_DIR}/salmon_${sample_id}/quant.sf"

        if [[ ! -f "${quant_file}" ]]; then
            echo "ERROR: Salmon quantification file not found:"
            echo "  ${quant_file}"
            exit 1
        fi

        cp "${quant_file}" \
            "${QAPA_PROJECT_DIR}/${sample_index}/quant.sf"
    done

# ============================================================================
# 21. Generate QAPAresults.txt
# ============================================================================

echo "========================================================================"
echo "Generating QAPAresults.txt"
echo "========================================================================"

(
    cd "${QAPA_PROJECT_DIR}"

    qapa quant \
        --db "${ENSEMBL_IDENTIFIERS}" \
        ./sample*/quant.sf \
        > ./QAPAresults.txt
)

# ============================================================================
# 22. Final output summary
# ============================================================================

echo
echo "========================================================================"
echo "Pipeline completed successfully"
echo "========================================================================"

echo "Main output directories:"
echo "  Raw FASTQ:       ${RAW_FASTQ_DIR}"
echo "  Clean FASTQ:     ${CUTADAPT_DIR}"
echo "  rRNA-depleted:   ${RRNA_DIR}"
echo "  Sorted BAM:      ${BAM_DIR}"
echo "  m6A results:     ${M6A_DIR}"
echo "  QAPA results:    ${QAPA_PROJECT_DIR}"

echo
echo "Important files:"
echo "  m6Aexpress script:"
echo "    ${m6a_r_script}"
echo
echo "  QAPA result table:"
echo "    ${QAPA_PROJECT_DIR}/QAPAresults.txt"
echo
echo "Finished at: $(date)"
