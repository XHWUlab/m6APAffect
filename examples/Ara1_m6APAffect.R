library(m6APAreg)
library(movAPA)

# Set the directory containing the Arabidopsis input files.
# Replace this path with the corresponding directory on your computer.
project_dir <- "path/to/Ara1"

setwd(project_dir)

# Input files:
#   peak.rda
#   all.PAC.PATcount
#   all.PAC.header
gtf_file <- "path/to/Arabidopsis_thaliana.TAIR10.62.gtf.gz"

# -------------------------------------------------------------------------
# 1. Process m6A-seq peak data
# -------------------------------------------------------------------------
# Only control m6A-seq data are available for this Arabidopsis dataset.
# Therefore, the treated m6A samples are not provided and treats is NULL.
load("peak.rda")

peak_info <- peak
peak_df <- peak_info[[1]]
lib_sizes <- peak_info[[2]]

ctrls <- c("IP1", "IP2", "Input1", "Input2")
treats <- NULL

# Check whether all required m6A-seq sample columns are available.
missing_m6a_samples <- setdiff(ctrls, colnames(peak_df))
if (length(missing_m6a_samples) > 0) {
  stop(
    "The following m6A-seq sample columns are missing: ",
    paste(missing_m6a_samples, collapse = ", ")
  )
}

# Reorder library sizes according to the sample order required by the model.
if (length(lib_sizes) != length(ctrls)) {
  stop("The length of lib_sizes must equal the number of control samples.")
}

names(lib_sizes) <- ctrls
lib_sizes <- lib_sizes[ctrls]

# Convert the peak table into a m6A PACdataset.
# DE testing and DE-based filtering are disabled because only control
# m6A-seq data are available.
m6ads <- m6APAreg::m6AExpress2PACds(
  peakDf = peak_df,
  ctrls = ctrls,
  treats = treats,
  libSizes = lib_sizes,
  doDE = FALSE,
  filterDE = FALSE
)

# Make chromosome names consistent with the APA annotation.
# Modify this step if the GTF file uses the "chr" prefix.
m6ads@anno$chr <- sub("^chr", "", m6ads@anno$chr)

# Rename the two control IP samples for the M3 control group.
m6ads_m3 <- m6APAreg::setPACdsSmpInfo(
  m6ads,
  smpInfo = cbind(
    old = c("IP1", "IP2"),
    new = c("ctrlM31", "ctrlM32"),
    group = c("ctrlM3", "ctrlM3")
  )
)

# Annotate m6A sites with the reference GTF file.
m6ads_m3 <- annotatePAC(m6ads_m3, gtf_file)

# Arabidopsis gene names in this dataset contain additional transcript
# information. Keep the first nine characters as the gene identifier.
m6ads_m3@anno$gene <- substring(m6ads_m3@anno$gene, 1, 9)

# geneID is not needed for the downstream gene-level m6A calculation.
m6ads_m3@anno$geneID <- NULL

# -------------------------------------------------------------------------
# 2. Process plantAPAdb poly(A) site data
# -------------------------------------------------------------------------
# Read the poly(A) site count matrix and its sample header.
apa_counts <- read.table(
  file.path(project_dir, "all.PAC.PATcount"),
  header = FALSE,
  stringsAsFactors = FALSE
)

apa_header <- read.table(
  file.path(project_dir, "all.PAC.header"),
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(apa_counts) <- as.character(apa_header[[1]])

# The last four columns contain the four APA sample count libraries.
colnames(apa_counts)[9:12] <- paste0(
  "sample", 1:4, ".CNT"
)

# Convert the APA table into a PACdataset.
apa_pac <- readPACds(apa_counts)

# Annotate APA sites using the same GTF file.
apa_pac <- annotatePAC(apa_pac, gtf_file)

# Keep the same Arabidopsis gene identifier format used for m6A sites.
apa_pac@anno$gene <- substring(apa_pac@anno$gene, 1, 9)

# Combine the required annotation columns and APA count columns.
apa <- cbind(
  apa_pac@anno[, c(1:6, 10)],
  apa_pac@counts[, 3:6]
)

# Filter APA sites by total expression.
# Sites above the 95th percentile are retained, following the original
# Arabidopsis analysis workflow.
apa$count <- rowSums(apa[, 8:11])
apa <- apa[
  apa$count > quantile(apa$count, probs = 0.95),
  ,
  drop = FALSE
]

apa$count <- NULL
apa$gene_QAPA <- apa$gene
rownames(apa) <- seq_len(nrow(apa))

# -------------------------------------------------------------------------
# 3. Convert APA counts into an APA PACdataset
# -------------------------------------------------------------------------
apa_input <- apa
count_suffix <- "CNT"

count_columns <- grep(
  paste0(count_suffix, "$"),
  colnames(apa_input)
)

if (length(count_columns) != 4) {
  stop("Exactly four APA count columns ending with 'CNT' are required.")
}

# Remove the ".CNT" suffix before constructing the sample metadata.
old_sample_names <- colnames(apa_input)[count_columns]
new_sample_names <- sub(
  paste0("\\.", count_suffix, "$"),
  "",
  old_sample_names
)

colnames(apa_input)[count_columns] <- new_sample_names
apa_input[, count_columns] <- lapply(
  apa_input[, count_columns, drop = FALSE],
  round
)

# Remove APA sites with zero counts in all samples.
zero_count_rows <- which(
  rowSums(apa_input[, count_columns, drop = FALSE]) == 0
)

if (length(zero_count_rows) > 0) {
  apa_input <- apa_input[-zero_count_rows, , drop = FALSE]
  count_columns <- match(
    new_sample_names,
    colnames(apa_input)
  )

  cat(
    "Removed APA sites with zero counts in all samples: ",
    length(zero_count_rows),
    "\n",
    sep = ""
  )
}

apa_col_data <- data.frame(
  group = colnames(apa_input)[count_columns],
  row.names = colnames(apa_input)[count_columns],
  stringsAsFactors = FALSE
)

apa_pac <- readPACds(
  apa_input,
  apa_col_data,
  PAname = "coord"
)

# The APA dataset contains two control and two treated sample columns.
# The treated APA samples are available even though treated m6A-seq data
# are unavailable.
apads_m3 <- m6APAreg::setPACdsSmpInfo(
  apa_pac,
  smpInfo = cbind(
    old = paste0("sample", 1:4),
    new = c("ctrlM31", "ctrlM32", "treatM31", "treatM32"),
    group = c("ctrlM3", "ctrlM3", "treatM3", "treatM3")
  )
)

# Standardize the coordinate column names before annotation.
colnames(apads_m3@anno)[4:5] <- c("start", "end")

apads_m3 <- annotatePAC(apads_m3, gtf_file)

# Use the original QAPA gene identifier for downstream APA analysis.
apads_m3@anno$gene <- apads_m3@anno$gene_QAPA

# Retain poly(A) sites located in the annotated 3' UTR.
apads_m3 <- ext3UTRPACds(apads_m3, 300)

# Extract 3' UTR APA sites and define proximal-distal APA pairs.
apads_utr_m3 <- movAPA::get3UTRAPAds(apads_m3)

apads_pd_m3 <- get3UTRAPApd(
  apads_utr_m3,
  minDist = 50,
  maxDist = 5000,
  minRatio = 0,
  fixDistal = FALSE,
  addCols = "pd"
)

# Identify differential APA sites between the control and treated groups.
apads_utr_m3_de <- get3UTRAPApdDE(
  apads_pd_m3,
  pthd = 0.05,
  filterDE = TRUE
)

# Calculate relative usage difference at the gene level.
p_m3 <- getRUDperGene(apads_utr_m3_de)

# -------------------------------------------------------------------------
# 4. Filter m6A sites and calculate gene-level m6A values
# -------------------------------------------------------------------------
# Apply the original 75th-percentile expression filter to m6A sites.
m6a_sample_columns <- c(
  "IP1",
  "IP2",
  "Input1",
  "Input2"
)

missing_m6a_count_columns <- setdiff(
  m6a_sample_columns,
  colnames(m6ads_m3@anno)
)

if (length(missing_m6a_count_columns) > 0) {
  stop(
    "The following m6A count columns are missing from m6ads_m3@anno: ",
    paste(missing_m6a_count_columns, collapse = ", ")
  )
}

m6ads_m3@anno$count <- rowSums(
  m6ads_m3@anno[, m6a_sample_columns, drop = FALSE]
)

m6a_keep <- m6ads_m3@anno$count >
  quantile(m6ads_m3@anno$count, probs = 0.75)

m6ads_m3@anno <- m6ads_m3@anno[m6a_keep, , drop = FALSE]
m6ads_m3@counts <- m6ads_m3@counts[m6a_keep, , drop = FALSE]

m6ads_m3@anno$count <- NULL

# Calculate gene-level m6A values for genes with differential APA.
m_m3 <- getM6AperGene(
  m6ads_m3,
  apads_utr_m3_de
)

# -------------------------------------------------------------------------
# 5. Merge m6A and APA measurements
# -------------------------------------------------------------------------
mp_all_m3 <- merge(
  m_m3,
  p_m3,
  by.x = "gene",
  by.y = "gene"
)

# No treated m6A-seq data are available for this dataset.
# These zero values are placeholders required to match the model input
# structure; they do not represent measured treated m6A signals.
mp_all_m3$treatM31_m6A <- 0
mp_all_m3$treatM32_m6A <- 0

# Ensure model variables are numeric.
numeric_columns <- setdiff(
  colnames(mp_all_m3),
  "gene"
)

mp_all_m3[numeric_columns] <- lapply(
  mp_all_m3[numeric_columns],
  as.numeric
)

# -------------------------------------------------------------------------
# 6. Fit the m6APAreg model
# -------------------------------------------------------------------------
models_fit_m3 <- fitQuasiGLM(mp_all_m3)

models_m3 <- addPost2FitModels(
  models_fit_m3,
  robust = FALSE
)

res_m3 <- testBetas(models_m3)

res_mp_m3 <- merge(
  res_m3,
  mp_all_m3,
  by.x = "gene",
  by.y = "gene"
)

# Calculate p-values and q-values for the m6APAreg test.
res_mp_m3 <- m6APAreg::.statPQ(
  res_mp_m3,
  pthd = 0.05,
  qthd = 0.05
)

# Keep the columns required for reporting.
res_mp_m3 <- m6APAreg::.subsetResMp(res_mp_m3)

# Select significant m6APA-regulated genes.
res_mp_m3_qv <- m6APAreg::.subsetResMp(
  res_mp_m3,
  qval1 = 0.05
)

# Save results inside the project directory.
write.csv(
  res_mp_m3_qv,
  file = file.path(project_dir, "Ara1.m6APAreg.csv"),
  row.names = FALSE
)

# -------------------------------------------------------------------------
# 7. Summarize the results
# -------------------------------------------------------------------------
# The following calculation independently verifies the number of genes
# with q-value < 0.05.
qvalue_count <- sum(
  qvalue::qvalue(res_mp_m3$pval1)$qvalues < 0.05,
  na.rm = TRUE
)

# Preserve the complete result table and apply the significance threshold.
res_mp_m3_limited <- res_mp_m3[
  !is.na(res_mp_m3$qval1) &
    res_mp_m3$qval1 < 0.05,
  ,
  drop = FALSE
]

# Count positive and negative regulatory effects.
positive_regulators <- sum(
  res_mp_m3_limited$b1 > 0,
  na.rm = TRUE
)

negative_regulators <- sum(
  res_mp_m3_limited$b1 < 0,
  na.rm = TRUE
)

# Differential APA genes.
DP <- unique(apads_utr_m3_de@anno$gene_QAPA)

# Compare differential APA genes with significant m6APAffect genes.
DP_m6APAreg <- intersect(
  DP,
  res_mp_m3_limited$gene
)

nonm6APAreg <- setdiff(
  DP,
  res_mp_m3_limited$gene
)

# Save gene lists for downstream comparison.
write.table(
  DP_m6APAreg,
  file = file.path(project_dir, "DP_m6APAreg_genes.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  nonm6APAreg,
  file = file.path(project_dir, "nonm6APAreg_genes.txt"),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Print a compact summary.
cat("Genes with q-value < 0.05:", qvalue_count, "\n")
cat("Positive regulatory effects:", positive_regulators, "\n")
cat("Negative regulatory effects:", negative_regulators, "\n")
cat("Differential APA genes:", length(DP), "\n")
cat("DP genes also detected as m6APAffect genes:", length(DP_m6APAreg), "\n")
cat("DP genes not detected as m6APAffect genes:", length(nonm6APAreg), "\n")

