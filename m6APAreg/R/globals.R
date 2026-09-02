utils::globalVariables(c(
  "Get_peak_infor",
  "UPA_end",
  "UPA_start",
  "chr",
  "coord",
  "subjectHits",
  "tottag",
  "strand"
))

.isChrConsistent <- function(pacds, obj, col = NULL, allin = FALSE) {
  if (is.vector(pacds)) {
    pacchr <- pacds
  } else if (
    inherits(pacds, "data.frame") ||
    inherits(obj, "matrix")
  ) {
    pacchr <- unique(as.character(pacds$chr))
  } else if (inherits(pacds, "PACdataset")) {
    pacchr <- unique(as.character(pacds@anno$chr))
  } else {
    stop("Unsupported pacds object.")
  }

  if (inherits(obj, "BSgenome")) {
    objchr <- as.character(GenomeInfoDb::seqnames(obj))
  } else if (inherits(obj, "FaFile")) {
    objchr <- names(GenomeInfoDb::seqlengths(obj))
  } else if (
    inherits(obj, "data.frame") ||
    inherits(obj, "matrix")
  ) {
    if (is.null(col)) {
      if ("seqnames" %in% colnames(obj)) {
        col <- "seqnames"
      } else if ("chr" %in% colnames(obj)) {
        col <- "chr"
      }
    }

    if (is.null(col) || !(col %in% colnames(obj))) {
      stop("Chromosome column was not found in obj.")
    }

    objchr <- unique(as.character(obj[, col]))
  } else {
    objchr <- as.character(obj)
  }

  if (isTRUE(allin)) {
    return(all(as.character(pacchr) %in% as.character(objchr)))
  }

  length(base::intersect(
    as.character(pacchr),
    as.character(objchr)
  )) > 0
}
