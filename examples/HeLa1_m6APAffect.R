library(m6APAreg)
library(movAPA)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)

# =========================================================================
# Configuration
# =========================================================================

# Replace these paths with the directories containing the HeLa1 input files.
project_dir <- "path/to/HeLa1"
qapa_root <- "path/to/HeLa1"
m6a_dir <- project_dir

setwd(project_dir)

# Use the hg38 annotation system for both APA and m6A datasets.
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# =========================================================================
# 1. Add Salmon counts to the QAPA results
# =========================================================================

# M3 KDM3 knockdown
m6APAreg::addQAPARawCounts(
  qfile = file.path(qapa_root, "M3_KD_LF", "QAPAresults.txt"),
  sfs = setNames(
    file.path(
      qapa_root,
      "M3_KD_LF",
      paste0("sample", c(5, 6, 9, 10)),
      "quant.sf"
    ),
    paste0("sample", c(5, 6, 9, 10))
  ),
  ofile = file.path(project_dir, "M3KD_qapa_CNT.txt")
)

# M14 KDM14 knockdown
m6APAreg::addQAPARawCounts(
  qfile = file.path(qapa_root, "M14_KD_LF", "QAPAresults.txt"),
  sfs = setNames(
    file.path(
      qapa_root,
      "M14_KD_LF",
      paste0("sample", c(1, 2, 7, 8)),
      "quant.sf"
    ),
    paste0("sample", c(1, 2, 7, 8))
  ),
  ofile = file.path(project_dir, "M14KD_qapa_CNT.txt")
)

# Combined M3/M14 knockdown
m6APAreg::addQAPARawCounts(
  qfile = file.path(qapa_root, "M3M14_KD_LF", "QAPAresults.txt"),
  sfs = setNames(
    file.path(
      qapa_root,
      "M3M14_KD_LF",
      paste0(
        "sample",
        c(1, 2, 5, 6, 7, 8, 9, 10)
      ),
      "quant.sf"
    ),
    paste0(
      "sample",
      c(1, 2, 5, 6, 7, 8, 9, 10)
    )
  ),
  ofile = file.path(project_dir, "M3M14KD_qapa_CNT.txt")
)

# WTAP knockdown
m6APAreg::addQAPARawCounts(
  qfile = file.path(qapa_root, "WTAP_KD_LF", "QAPAresults.txt"),
  sfs = setNames(
    file.path(
      qapa_root,
      "WTAP_KD_LF",
      paste0("sample", c(3, 4, 7, 8)),
      "quant.sf"
    ),
    paste0("sample", c(3, 4, 7, 8))
  ),
  ofile = file.path(project_dir, "WTAPKD_qapa_CNT.txt")
)

# =========================================================================
# 2. Convert QAPA count tables into APA PACdatasets
# =========================================================================

# M3 KDM3 knockdown
qapa_m3 <- read.table(
  file.path(project_dir, "M3KD_qapa_CNT.txt"),
  header = TRUE,
  sep = "\t",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("M3 QAPA rows: ", nrow(qapa_m3), "\n", sep = "")
cat("M3 QAPA columns: ", ncol(qapa_m3), "\n", sep = "")

apads_m3 <- m6APAreg::QAPA2PACds(
  qapa_m3,
  vcol = "CNT"
)

apads_m3 <- m6APAreg::setPACdsSmpInfo(
  apads_m3,
  smpInfo = cbind(
    old = paste0("sample", c(5, 6, 9, 10)),
    new = c(
      "KDM3_1",
      "KDM3_2",
      "ctrlM3_1",
      "ctrlM3_2"
    ),
    group = c(
      "KDM3",
      "KDM3",
      "ctrlM3",
      "ctrlM3"
    )
  )
)

apads_m3@anno$gene_QAPA <- apads_m3@anno$gene

# M14 KDM14 knockdown
qapa_m14 <- read.table(
  file.path(project_dir, "M14KD_qapa_CNT.txt"),
  header = TRUE,
  sep = "\t",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("M14 QAPA rows: ", nrow(qapa_m14), "\n", sep = "")
cat("M14 QAPA columns: ", ncol(qapa_m14), "\n", sep = "")

apads_m14 <- m6APAreg::QAPA2PACds(
  qapa_m14,
  vcol = "CNT"
)

apads_m14 <- m6APAreg::setPACdsSmpInfo(
  apads_m14,
  smpInfo = cbind(
    old = paste0("sample", c(1, 2, 7, 8)),
    new = c(
      "KDM14_1",
      "KDM14_2",
      "ctrlM14_1",
      "ctrlM14_2"
    ),
    group = c(
      "KDM14",
      "KDM14",
      "ctrlM14",
      "ctrlM14"
    )
  )
)

apads_m14@anno$gene_QAPA <- apads_m14@anno$gene

# Combined M3/M14 knockdown
qapa_m3m14 <- read.table(
  file.path(project_dir, "M3M14KD_qapa_CNT.txt"),
  header = TRUE,
  sep = "\t",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Combined QAPA rows: ", nrow(qapa_m3m14), "\n", sep = "")
cat("Combined QAPA columns: ", ncol(qapa_m3m14), "\n", sep = "")

apads_m3m14 <- m6APAreg::QAPA2PACds(
  qapa_m3m14,
  vcol = "CNT"
)

apads_m3m14 <- m6APAreg::setPACdsSmpInfo(
  apads_m3m14,
  smpInfo = cbind(
    old = paste0(
      "sample",
      c(1, 2, 5, 6, 7, 8, 9, 10)
    ),
    new = c(
      "KDM14_1",
      "KDM14_2",
      "KDM3_1",
      "KDM3_2",
      "ctrlM14_1",
      "ctrlM14_2",
      "ctrlM3_1",
      "ctrlM3_2"
    ),
    group = c(
      "KDM14",
      "KDM14",
      "KDM3",
      "KDM3",
      "ctrlM14",
      "ctrlM14",
      "ctrlM3",
      "ctrlM3"
    )
  )
)

apads_m3m14@anno$gene_QAPA <- apads_m3m14@anno$gene

# WTAP knockdown
qapa_wtap <- read.table(
  file.path(project_dir, "WTAPKD_qapa_CNT.txt"),
  header = TRUE,
  sep = "\t",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("WTAP QAPA rows: ", nrow(qapa_wtap), "\n", sep = "")
cat("WTAP QAPA columns: ", ncol(qapa_wtap), "\n", sep = "")

apads_wtap <- m6APAreg::QAPA2PACds(
  qapa_wtap,
  vcol = "CNT"
)

apads_wtap <- m6APAreg::setPACdsSmpInfo(
  apads_wtap,
  smpInfo = cbind(
    old = paste0("sample", c(3, 4, 7, 8)),
    new = c(
      "KDWTAP_1",
      "KDWTAP_2",
      "ctrlWTAP_1",
      "ctrlWTAP_2"
    ),
    group = c(
      "KDWTAP",
      "KDWTAP",
      "ctrlWTAP",
      "ctrlWTAP"
    )
  )
)

apads_wtap@anno$gene_QAPA <- apads_wtap@anno$gene

# =========================================================================
# 3. Convert m6A data into m6A PACdatasets
# =========================================================================

# M3 KDM3 knockdown
m6a_m3_data <- m6APAreg::.getRdaData(
  file.path(m6a_dir, "M3KD-Hela1.rda")
)

m6ads_m3 <- m6APAreg::m6AExpress2PACds(
  peakDf = m6a_m3_data[[1]],
  ctrls = m6a_m3_data[[2]],
  treats = m6a_m3_data[[3]],
  libSizes = m6a_m3_data[[4]],
  doDE = TRUE,
  filterDE = FALSE
)

m6ads_m3 <- m6APAreg::setPACdsSmpInfo(
  m6ads_m3,
  smpInfo = cbind(
    old = c(
      "IP1",
      "IP2",
      "Treated_IP1",
      "Treated_IP2"
    ),
    new = c(
      "ctrlM3_1",
      "ctrlM3_2",
      "KDM3_1",
      "KDM3_2"
    ),
    group = c(
      "ctrlM3",
      "ctrlM3",
      "KDM3",
      "KDM3"
    )
  )
)

# M14 KDM14 knockdown
m6a_m14_data <- m6APAreg::.getRdaData(
  file.path(m6a_dir, "M14KD-Hela1.rda")
)

m6ads_m14 <- m6APAreg::m6AExpress2PACds(
  peakDf = m6a_m14_data[[1]],
  ctrls = m6a_m14_data[[2]],
  treats = m6a_m14_data[[3]],
  libSizes = m6a_m14_data[[4]],
  doDE = TRUE,
  filterDE = FALSE
)

m6ads_m14 <- m6APAreg::setPACdsSmpInfo(
  m6ads_m14,
  smpInfo = cbind(
    old = c(
      "IP1",
      "IP2",
      "Treated_IP1",
      "Treated_IP2"
    ),
    new = c(
      "ctrlM14_1",
      "ctrlM14_2",
      "KDM14_1",
      "KDM14_2"
    ),
    group = c(
      "ctrlM14",
      "ctrlM14",
      "KDM14",
      "KDM14"
    )
  )
)

# WTAP knockdown
m6a_wtap_data <- m6APAreg::.getRdaData(
  file.path(m6a_dir, "WTAPKD-Hela1.rda")
)

m6ads_wtap <- m6APAreg::m6AExpress2PACds(
  peakDf = m6a_wtap_data[[1]],
  ctrls = m6a_wtap_data[[2]],
  treats = m6a_wtap_data[[3]],
  libSizes = m6a_wtap_data[[4]],
  doDE = TRUE,
  filterDE = FALSE
)

m6ads_wtap <- m6APAreg::setPACdsSmpInfo(
  m6ads_wtap,
  smpInfo = cbind(
    old = c(
      "IP1",
      "IP2",
      "Treated_IP1",
      "Treated_IP2"
    ),
    new = c(
      "ctrlWTAP_1",
      "ctrlWTAP_2",
      "KDWTAP_1",
      "KDWTAP_2"
    ),
    group = c(
      "ctrlWTAP",
      "ctrlWTAP",
      "KDWTAP",
      "KDWTAP"
    )
  )
)

# =========================================================================
# 4. Filter DE-m6A sites
# =========================================================================

# M3 DE-m6A genes
de_genes_m3 <- unique(
  m6ads_m3@anno$gene[
    !is.na(m6ads_m3@anno$DE_pvalue) &
      m6ads_m3@anno$DE_pvalue < 0.05
  ]
)

m6ads_m3_de <- m6APAreg::subsetPACds(
  m6ads_m3,
  genes = de_genes_m3
)

# M14 DE-m6A genes
de_genes_m14 <- unique(
  m6ads_m14@anno$gene[
    !is.na(m6ads_m14@anno$DE_pvalue) &
      m6ads_m14@anno$DE_pvalue < 0.05
  ]
)

m6ads_m14_de <- m6APAreg::subsetPACds(
  m6ads_m14,
  genes = de_genes_m14
)

# WTAP DE-m6A genes
de_genes_wtap <- unique(
  m6ads_wtap@anno$gene[
    !is.na(m6ads_wtap@anno$DE_pvalue) &
      m6ads_wtap@anno$DE_pvalue < 0.05
  ]
)

m6ads_wtap_de <- m6APAreg::subsetPACds(
  m6ads_wtap,
  genes = de_genes_wtap
)

# =========================================================================
# 5. Annotate APA and m6A datasets using hg38
# =========================================================================

# APA datasets
apads_m3 <- annotatePAC(apads_m3, txdb)
apads_m14 <- annotatePAC(apads_m14, txdb)
apads_wtap <- annotatePAC(apads_wtap, txdb)

# m6A datasets
m6ads_m3 <- annotatePAC(m6ads_m3, txdb)
m6ads_m14 <- annotatePAC(m6ads_m14, txdb)
m6ads_wtap <- annotatePAC(m6ads_wtap, txdb)

m6ads_m3_de <- annotatePAC(m6ads_m3_de, txdb)
m6ads_m14_de <- annotatePAC(m6ads_m14_de, txdb)
m6ads_wtap_de <- annotatePAC(m6ads_wtap_de, txdb)

# Extend the 3' UTR annotation by 2000 nt.
apads_m3 <- ext3UTRPACds(apads_m3, 2000)
apads_m14 <- ext3UTRPACds(apads_m14, 2000)
apads_wtap <- ext3UTRPACds(apads_wtap, 2000)

m6ads_m3 <- ext3UTRPACds(m6ads_m3, 2000)
m6ads_m14 <- ext3UTRPACds(m6ads_m14, 2000)
m6ads_wtap <- ext3UTRPACds(m6ads_wtap, 2000)

m6ads_m3_de <- ext3UTRPACds(m6ads_m3_de, 2000)
m6ads_m14_de <- ext3UTRPACds(m6ads_m14_de, 2000)
m6ads_wtap_de <- ext3UTRPACds(m6ads_wtap_de, 2000)

# Remove m6A sites that cannot be assigned to a gene.
m6ads_m3 <- m6ads_m3[
  !is.na(m6ads_m3@anno$gene)
]

m6ads_m14 <- m6ads_m14[
  !is.na(m6ads_m14@anno$gene)
]

m6ads_wtap <- m6ads_wtap[
  !is.na(m6ads_wtap@anno$gene)
]

m6ads_m3_de <- m6ads_m3_de[
  !is.na(m6ads_m3_de@anno$gene)
]

m6ads_m14_de <- m6ads_m14_de[
  !is.na(m6ads_m14_de@anno$gene)
]

m6ads_wtap_de <- m6ads_wtap_de[
  !is.na(m6ads_wtap_de@anno$gene)
]

# =========================================================================
# 6. Identify differential APA events
# =========================================================================

# M3 differential APA
apads_m3_utr <- movAPA::get3UTRAPAds(apads_m3)

apads_m3_pd <- get3UTRAPApd(
  apads_m3_utr,
  minDist = 50,
  maxDist = 5000,
  minRatio = 0.05,
  fixDistal = FALSE,
  addCols = "pd"
)

apads_m3_de <- get3UTRAPApdDE(
  apads_m3_pd,
  pthd = 0.05,
  filterDE = TRUE
)

# M14 differential APA
apads_m14_utr <- movAPA::get3UTRAPAds(apads_m14)

apads_m14_pd <- get3UTRAPApd(
  apads_m14_utr,
  minDist = 50,
  maxDist = 5000,
  minRatio = 0.05,
  fixDistal = FALSE,
  addCols = "pd"
)

apads_m14_de <- get3UTRAPApdDE(
  apads_m14_pd,
  pthd = 0.05,
  filterDE = TRUE
)

# WTAP differential APA
apads_wtap_utr <- movAPA::get3UTRAPAds(apads_wtap)

apads_wtap_pd <- get3UTRAPApd(
  apads_wtap_utr,
  minDist = 50,
  maxDist = 5000,
  minRatio = 0.05,
  fixDistal = FALSE,
  addCols = "pd"
)

apads_wtap_de <- get3UTRAPApdDE(
  apads_wtap_pd,
  pthd = 0.05,
  filterDE = TRUE
)

# Save intermediate objects for reproducibility.
save(
  m6ads_m3,
  m6ads_m14,
  m6ads_wtap,
  apads_m3,
  apads_m14,
  apads_wtap,
  file = file.path(project_dir, "m6A_APA_HeLa1.rda")
)

save(
  m6ads_m3_de,
  m6ads_m14_de,
  m6ads_wtap_de,
  apads_m3_de,
  apads_m14_de,
  apads_wtap_de,
  file = file.path(project_dir, "DE_m6A_APA_HeLa1.rda")
)

# =========================================================================
# 7. Calculate gene-level m6A scores and APA usage
# =========================================================================

# M3 gene-level m6A and APA data
p_m3 <- getRUDperGene(apads_m3_de)

m_m3 <- getM6AperGene(
  m6ads_m3_de,
  apads_m3_de
)

mp_all_m3 <- merge(
  m_m3,
  p_m3,
  by = "gene"
)

# M14 gene-level m6A and APA data
p_m14 <- getRUDperGene(apads_m14_de)

m_m14 <- getM6AperGene(
  m6ads_m14_de,
  apads_m14_de
)

mp_all_m14 <- merge(
  m_m14,
  p_m14,
  by = "gene"
)

# WTAP gene-level m6A and APA data
p_wtap <- getRUDperGene(apads_wtap_de)

m_wtap <- getM6AperGene(
  m6ads_wtap_de,
  apads_wtap_de
)

mp_all_wtap <- merge(
  m_wtap,
  p_wtap,
  by = "gene"
)

# =========================================================================
# 8. Prepare Entrez-to-gene-symbol mapping
# =========================================================================

gene_annotation <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = AnnotationDbi::keys(
    org.Hs.eg.db,
    keytype = "ENTREZID"
  ),
  columns = c(
    "ENTREZID",
    "ENSEMBL",
    "SYMBOL"
  ),
  keytype = "ENTREZID"
)

colnames(gene_annotation) <- c(
  "gene_entrezid",
  "gene_ensembl",
  "gene_symbol"
)

gene_annotation <- gene_annotation[
  !duplicated(gene_annotation$gene_entrezid),
  ,
  drop = FALSE
]

# =========================================================================
# 9. Fit the m6APAreg model for each knockdown condition
# =========================================================================

# M3 model
mp_all_m3[
  setdiff(colnames(mp_all_m3), "gene")
] <- lapply(
  mp_all_m3[
    setdiff(colnames(mp_all_m3), "gene")
  ],
  as.numeric
)

models_fit_m3 <- fitQuasiGLM(mp_all_m3)

models_m3 <- addPost2FitModels(
  models_fit_m3,
  robust = FALSE
)

result_m3 <- testBetas(models_m3)

result_m3 <- merge(
  result_m3,
  mp_all_m3,
  by = "gene"
)

result_m3 <- m6APAreg::.statPQ(
  result_m3,
  pthd = 0.05,
  qthd = 0.2
)

result_m3 <- m6APAreg::.subsetResMp(result_m3)

result_m3 <- merge(
  result_m3,
  gene_annotation,
  by.x = "gene",
  by.y = "gene_entrezid",
  all.x = TRUE,
  all.y = FALSE
)

result_m3_significant <- m6APAreg::.subsetResMp(
  result_m3,
  pval1 = 0.05
)

write.csv(
  result_m3,
  file = file.path(project_dir, "M3KD_HeLa1.resMp.csv"),
  row.names = FALSE
)

write.csv(
  result_m3_significant,
  file = file.path(project_dir, "M3KD_HeLa1.m6APAreg.csv"),
  row.names = FALSE
)

# M14 model
mp_all_m14[
  setdiff(colnames(mp_all_m14), "gene")
] <- lapply(
  mp_all_m14[
    setdiff(colnames(mp_all_m14), "gene")
  ],
  as.numeric
)

models_fit_m14 <- fitQuasiGLM(mp_all_m14)

models_m14 <- addPost2FitModels(
  models_fit_m14,
  robust = FALSE
)

result_m14 <- testBetas(models_m14)

result_m14 <- merge(
  result_m14,
  mp_all_m14,
  by = "gene"
)

result_m14 <- m6APAreg::.statPQ(
  result_m14,
  pthd = 0.05,
  qthd = 0.2
)

result_m14 <- m6APAreg::.subsetResMp(result_m14)

result_m14 <- merge(
  result_m14,
  gene_annotation,
  by.x = "gene",
  by.y = "gene_entrezid",
  all.x = TRUE,
  all.y = FALSE
)

result_m14_significant <- m6APAreg::.subsetResMp(
  result_m14,
  pval1 = 0.05
)

write.csv(
  result_m14,
  file = file.path(project_dir, "M14KD_HeLa1.resMp.csv"),
  row.names = FALSE
)

write.csv(
  result_m14_significant,
  file = file.path(project_dir, "M14KD_HeLa1.m6APAreg.csv"),
  row.names = FALSE
)

# WTAP model
mp_all_wtap[
  setdiff(colnames(mp_all_wtap), "gene")
] <- lapply(
  mp_all_wtap[
    setdiff(colnames(mp_all_wtap), "gene")
  ],
  as.numeric
)

models_fit_wtap <- fitQuasiGLM(mp_all_wtap)

models_wtap <- addPost2FitModels(
  models_fit_wtap,
  robust = FALSE
)

result_wtap <- testBetas(models_wtap)

result_wtap <- merge(
  result_wtap,
  mp_all_wtap,
  by = "gene"
)

result_wtap <- m6APAreg::.statPQ(
  result_wtap,
  pthd = 0.05,
  qthd = 0.2
)

result_wtap <- m6APAreg::.subsetResMp(result_wtap)

result_wtap <- merge(
  result_wtap,
  gene_annotation,
  by.x = "gene",
  by.y = "gene_entrezid",
  all.x = TRUE,
  all.y = FALSE
)

result_wtap_significant <- m6APAreg::.subsetResMp(
  result_wtap,
  pval1 = 0.05
)

write.csv(
  result_wtap,
  file = file.path(project_dir, "WTAPKD_HeLa1.resMp.csv"),
  row.names = FALSE
)

write.csv(
  result_wtap_significant,
  file = file.path(project_dir, "WTAPKD_HeLa1.m6APAreg.csv"),
  row.names = FALSE
)

# =========================================================================
# 10. Summarize the results
# =========================================================================

# M3 summary
cat("\nM3 KDM3 knockdown\n")

cat(
  "Genes with pval1 < 0.05: ",
  nrow(result_m3_significant),
  "\n",
  sep = ""
)

cat(
  "Positive regulatory effects: ",
  sum(result_m3_significant$b1 > 0, na.rm = TRUE),
  "\n",
  sep = ""
)

cat(
  "Negative regulatory effects: ",
  sum(result_m3_significant$b1 < 0, na.rm = TRUE),
  "\n",
  sep = ""
)

# M14 summary
cat("\nM14 KDM14 knockdown\n")

cat(
  "Genes with pval1 < 0.05: ",
  nrow(result_m14_significant),
  "\n",
  sep = ""
)

cat(
  "Positive regulatory effects: ",
  sum(result_m14_significant$b1 > 0, na.rm = TRUE),
  "\n",
  sep = ""
)

cat(
  "Negative regulatory effects: ",
  sum(result_m14_significant$b1 < 0, na.rm = TRUE),
  "\n",
  sep = ""
)

# WTAP summary
cat("\nWTAP knockdown\n")

cat(
  "Genes with pval1 < 0.05: ",
  nrow(result_wtap_significant),
  "\n",
  sep = ""
)

cat(
  "Positive regulatory effects: ",
  sum(result_wtap_significant$b1 > 0, na.rm = TRUE),
  "\n",
  sep = ""
)

cat(
  "Negative regulatory effects: ",
  sum(result_wtap_significant$b1 < 0, na.rm = TRUE),
  "\n",
  sep = ""
)
