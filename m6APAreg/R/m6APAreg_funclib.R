#' @importFrom magrittr %>%
#' @importFrom movAPA readPACds subsetPACds faFromPACds
#' @importFrom movAPA annotatePAC ext3UTRPACds get3UTRAPAds
#' @importFrom utils read.table head write.table combn read.csv
#' @importFrom methods new setClass
#' @importFrom GenomicRanges GRanges reduce resize mcols "mcols<-" findOverlaps
#' @importFrom IRanges IRanges subsetByOverlaps width
#' @importFrom GenomeInfoDb seqnames
#' @importFrom Biostrings vmatchPattern
#' @importFrom dplyr group_by summarise across everything n
#' @importFrom limma squeezeVar
#' @importFrom ChIPseeker annotatePeak
#' @importFrom qvalue qvalue
#' @importFrom stats sd glm fisher.test median quantile pnorm dnorm p.adjust pt
#' @importFrom S4Vectors elementNROWS
#' @importFrom stats start
#' @importFrom locfdr locfdr
#' @importFrom BiocIO import
NULL


# ----- utils ----

# validate.arg(NULL, c('a','b','c', 'no')) --> a
# validate.arg('x', c('a','b','c', 'no')) --> error
# validate.arg(NULL, c('a','b','c', 'no'), null2default=FALSE)  --> error
# validate.arg(NULL, c('a','b','c', 'no', NULL), null2default=FALSE)  --> error
validate.arg <- function(arg, choices, several.ok = FALSE, null2default = TRUE, lc = TRUE) {
  if (is.null(arg)) {
    if (null2default) arg=choices[1]
  }
  if (lc) arg=tolower(arg)
  arg = match.arg(arg, choices, several.ok)
  return(arg)
}


AinB<-function(A, B, all=T){
  if (is.null(A)) return(TRUE)
  if (length(A)==0) return(TRUE)
  #cat('all', all ,class(A), class(B), '\n')
  if (is.factor(A))   A=as.character(A)
  if (is.factor(B))   B=as.character(B)
  if (!(is.vector(A) & is.vector(B))) {
    return(F)
  }
  #print(A); print(B)
  x=sum(A %in% B)
  if ((all & x==length(A)) | (!all & x>0)) {
    return(T)
  } else {
    return(F)
  }
}


## format sample names into three columns df: old new group
## if group is not provided, then group=new
#.formatSmpInfo(c('smp1','smp2'))
#.formatSmpInfo(data.frame(cbind(1:4, 5:8, rep('a',4))))
#.formatSmpInfo(data.frame(cbind(1:4, 5:8)))
#.formatSmpInfo(data.frame(cbind(1:4, 'a')))
#.formatSmpInfo(data.frame(cbind(old=1:4, new='a'))) ##error
.formatSmpInfo<-function (smpInfo) {
  old=NA; new=NA; group=NA
  if (is.vector(smpInfo)) smpInfo=matrix(smpInfo, ncol=1)
  if (ncol(smpInfo)==1) {
    old=smpInfo[, 1]
    new=old
    group=new
  } else if (ncol(smpInfo)>=3) {
    old=smpInfo[, 1]
    new=smpInfo[, 2]
    group=smpInfo[, 3]
  } else if (ncol(smpInfo)==2) {
    if (AinB(c('old','new'), colnames(smpInfo))) {
      old=smpInfo[, 'old']
      new=smpInfo[, 'new']
      group=rep(1, length(old))
    } else if (AinB(c('old','group'), colnames(smpInfo))) {
      old=smpInfo[, 'old']
      new=old
      group=smpInfo[, 'group']
    } else {
      if ( sum(duplicated(smpInfo[, 2]))>0) {
        group=smpInfo[, 2]
        old=smpInfo[, 1]
        new=old
      } else {
        old=smpInfo[, 1]
        new=smpInfo[, 2]
        group=new
      }
    }
  }
  if (anyNA(c(old, new, group))) stop("Wrong smpInfo format, should be 2/3 columns, with old/new/group headers!\n")
  if ( sum(duplicated(old) | duplicated(new))>0) stop("Error smpInfo: duplicated old/new sample names!\n")
  return(data.frame(old, new, group))
}


.validSmpInfo<-function(d, smpInfo) {
  smpInfo=.formatSmpInfo(smpInfo)
  
  if (inherits(d, "data.frame")) {
    if (!AinB(smpInfo$old, colnames(d))) return(F)
  } else if (inherits(d, "PACdataset")) {
    if (!AinB(smpInfo$old, colnames(d@counts))) return(F)
  } else {
    stop("d should be data.frame or PACdataset!\n")
  }
  return(TRUE)
}

## ----------- any2PACds -------------

#' Set sample information for a PACdataset
#'
#' @param d A PACdataset object.
#' @param smpInfo A data.frame containing old, new, and group columns.
#' @return A PACdataset object.
#' @export
setPACdsSmpInfo<-function(d, smpInfo) {
  smpInfo=.formatSmpInfo(smpInfo)
  if (!AinB(colnames(d@counts), smpInfo$old)) {
    stop("setPACdsSmpInfo error: smpInfo$old not all in PACdataset d!\n")
  }
  rownames(d@colData)=smpInfo$new
  d@colData$group=smpInfo$group
  d@counts=d@counts[, smpInfo$old] #order
  colnames(d@counts)=smpInfo$new
  return(d)
}

# d= file or data.frame, .<vcol>$ is the sample columns
# vcol: TPM, CNT.. col in the QAPA file
# return: PACds with zero rows removed
# example:
# d=QAPA2PACds('sim3_qapa_hg38.txt')
# library('BSgenome.Hsapiens.UCSC.hg38')
# bsgenome=BSgenome.Hsapiens.UCSC.hg38
# fafiles=faFromPACds(d, bsgenome, what='updn', fapre='pac', byGrp = 'strand')
# plotATCGforFAfile(fafiles, mergePlots = TRUE, filepre='atcg')


#' convert QAPA results to PACdataset
#'
#' QAPA2PACds convert QAPA results to PACdataset.
#'
#' @param d A file name of QAPA output or a data.frame object.
#' The file or data.frame contains columns: Transcript Gene Gene_Name Chr Strand UTR3.Start UTR3.End, and sample columns like "sample1.PAU sample2.PAU..." or "sample1.TPM sample2.TPM ...".
#' @param vcol The string to extract sample columns. For example, vcol='TPM' means to get all columns containing 'TPM'.
#' @return PACdataset with all-zero rows removed.
#' @export
#' @examples
#' \dontrun{
#' apads=QAPA2PACds('sim3_qapa_hg38.txt', vcol='TPM')
#' }
#' @family APAdata functions
QAPA2PACds <- function(d, vcol = 'TPM') {
  
  if (is.character(d)) {
    d = read.table(d, header = TRUE)
    cat('read QAPA:', nrow(d), 'rows\n')
  }

  d = d[, c(
    'Chr', 'Strand', 'UTR3.End', 'UTR3.Start',
    'APA_ID', 'Transcript', 'Gene', 'Gene_Name',
    colnames(d)[grep(paste0(vcol, '$'), colnames(d))]
  )]
  
  d$coord = d$UTR3.End
  d$coord[d$Strand == '-'] = d$UTR3.Start[d$Strand == '-']
  
  cid = grep(paste0(vcol, '$'), colnames(d))
  
  colnames(d)[grep(paste0(vcol, '$'), colnames(d))] =
    gsub(
      paste0('.', vcol),
      '',
      colnames(d)[grep(paste0(vcol, '$'), colnames(d))],
      fixed = TRUE
    )
  
  colnames(d) = c(
    'chr', 'strand', 'UTR3_end', 'UTR3_start',
    'PAID', 'transcript', 'gene', 'gene_symbol',
    colnames(d)[cid], 'coord'
  )
  
  ## Round TPM values to integers
  d[, cid] = round(d[, cid])

  ## Remove rows in which all expression values are zero
  rid = which(rowSums(d[, cid]) == 0)
  if (length(rid) > 0) {
    d = d[-rid, ]
    cat('remove counts=0 rows:', length(rid), '\n')
  }
  
  rid = grep('chr', d$chr)
  if (length(rid) == 0) {
    cat('add chr to chromosome names: 1 --> chr1\n')
    d$chr = paste0('chr', d$chr)
  }
  
  colDataFile = as.data.frame(
    matrix(
      colnames(d)[cid],
      ncol = 1,
      dimnames = list(colnames(d)[cid], 'group')
    )
  )
  
  d = readPACds(d, colDataFile, PAname = 'coord')
  d@anno$gene_QAPA = d@anno$gene
  
  return(d)
}






#' convert APAtrap results to PACdataset
#'
#' APAtrap2PACds convert APAtrap results to PACdataset.
#'
#' @param d A file name of APAtrap output.
#' @return PACdataset.
#' @export
#' @examples
#' \dontrun{
#' apads=APAtrap2PACds('APAtrap.txt')
#' }
#' @family APAdata functions
APAtrap2PACds <- function(d) {
  
  if (is.character(d)) {
    d = read.table(d, header = TRUE)
    cat('read APAtrap:', nrow(d), 'rows\n')
  }
  
  d$chr = NULL
  d$strand = NULL
  d$gene = NULL
  d$gene_id = NULL
  d$GENE = NULL
  d$Mean_Squared_Error = NULL
  
  # Extract chromosome and strand information
  chrs = as.data.frame(
    matrix(
      unlist(strsplit(d$Gene, '|', fixed = TRUE)),
      ncol = 4,
      byrow = TRUE
    )
  )
  
  colnames(chrs) = c('transcript_id', 'gene', 'chr', 'strand')
  d = cbind(chrs, d)
  
  # Some genes have duplicated rows, such as RNVU1. These may be
  # non-coding RNAs; retain only the first occurrence.
  if (sum(duplicated(d$gene)) > 0) {
    cat(
      "Remove",
      sum(duplicated(d$gene)),
      'duplicated gene rows\n'
    )
  }
  
  d = d[!duplicated(d$gene), ]
  
  # Use the latter part of Loci as UTR_start and UTR_end.
  # Regardless of strand, UTR_end is the larger genomic coordinate.
  loci = gsub('.*:', '', d$Loci)
  
  loci = as.data.frame(
    matrix(
      unlist(strsplit(loci, '-', fixed = TRUE)),
      ncol = 2,
      byrow = TRUE
    )
  )
  
  head(loci)
  colnames(loci) = c('UTR_start', 'UTR_end')
  loci$UTR_start = as.numeric(loci$UTR_start)
  loci$UTR_end = as.numeric(loci$UTR_end)
  d = cbind(loci, d)
  
  dold = d
  
  # Use Predicted_APA as the PA coordinates. If a gene has multiple
  # PA sites, duplicate the other information for each PA site.
  # Add UTR_end for the plus strand or UTR_start for the minus strand
  # as the final PA site.
  d$Predicted_APA[d$strand == '+'] = paste0(
    d$Predicted_APA[d$strand == '+'],
    ',',
    d$UTR_end[d$strand == '+']
  )
  
  d$Predicted_APA[d$strand == '-'] = paste0(
    d$Predicted_APA[d$strand == '-'],
    ',',
    d$UTR_start[d$strand == '-']
  )
  
  PA = strsplit(d$Predicted_APA, ',', fixed = TRUE)
  nPA = unlist(lapply(PA, length))
  PA = unlist(PA)  # Individual PA sites
  
  d = d[rep(seq_len(nrow(d)), times = nPA), ]
  d$coord = as.numeric(PA)
  
  # Expression columns such as Group_1_1_Separate_Exp already
  # include expression values for distal PA sites.
  ec = colnames(d)[
    grep('Separate_Exp', colnames(d), fixed = TRUE)
  ]
  
  cat('Add', length(ec), 'expression columns\n')
  
  for (i in ec) {
    exp = dold[, i]
    
    # Some expression values are NA. Replace each NA with a
    # zero vector matching the number of PA sites, such as "0,0".
    naid = which(is.na(exp))
    
    if (length(naid) > 0) {
      for (j in 1:length(naid)) {
        exp[naid[j]] = paste0(
          rep(0, times = nPA[naid[j]]),
          collapse = ','
        )
      }
    }
    
    exp = strsplit(exp, ',', fixed = TRUE)
    # nexp = unlist(lapply(exp, length))
    exp = unlist(exp)  # Expression value for each PA site
    
    ec2 = gsub(
      '_Separate_Exp',
      '',
      i,
      fixed = TRUE
    )
    
    exp = as.numeric(exp)
    d = cbind(d, exp)
    colnames(d)[ncol(d)] = ec2
  }
  
  cn = gsub('_Separate_Exp', '', ec, fixed = TRUE)
  d$Gene = NULL
  d$Loci = NULL
  
  colDataFile = as.data.frame(
    matrix(
      cn,
      ncol = 1,
      dimnames = list(cn, 'group')
    )
  )
  
  d = readPACds(d, colDataFile)
  
  return(d)
}

#' convert m6ADB data to PACdataset
#'
#' m6ADB2PACds convert m6ADB data to PACdataset.
#'
#' @param d A file name of m6ADB data.
#' Columns starting with `RPM.` is the sample columns.
#' coord=peak_start(+)/peak_end(-), add will add peak_score, UPA_start and UPA_end.
#' @return PACdataset.
#' @export
#' @examples
#' \dontrun{
#' d=m6ADB2PACds('human_GSE46705_HeLa_M3KO_Ctrl.csv')
#' }
#' @family m6Adata functions
m6ADB2PACds <- function(d) {
  if (is.character(d)) {
    if (grepl('.csv', d, fixed=TRUE)) {
      d=read.csv(d, header=TRUE)
    } else {
      d=read.table(d, header=TRUE)
    }
    
    cat('read m6ADB data:', nrow(d), 'rows\n')
  }
  
  cid=grep('^RPM', colnames(d))
  
  d=d[, c('chr', 'strand', 'chromStart', 'chromEnd', 'geneID', 'diff.log2FC', 'pvalue',  'fdr', 'score', colnames(d)[cid])]
  cid=(ncol(d)-length(cid)+1):ncol(d)
  d$coord=d$chromEnd
  d$coord[d$Strand=='-']=d$chromStart[d$Strand=='-']
  
  colnames(d)=c('chr','strand','peak_start','peak_end', 'gene', 'log2fc','pvalue','fdr','score', colnames(d)[cid], 'coord')
  
  colnames(d)[grep('^RPM', colnames(d))]=gsub('RPM.', '', colnames(d)[grep('^RPM', colnames(d))], fixed=TRUE)

  d$peak_score=d$score
  d$UPA_start=d$peak_start
  d$UPA_end=d$peak_end
  
  colDataFile=as.data.frame(matrix(colnames(d)[cid], ncol=1, dimnames=list(colnames(d)[cid], 'group')))
  d=readPACds(d, colDataFile, PAname='coord')
  return(d)
}

## column name mappings
COLMAPS=data.frame(old=c('seqnames', 'start', 'end', 'width',  'gene_name', 'geneId', 'score','geneStart','geneEnd', 'geneStrand'),
                   new=c('chr','peak_start','peak_end','peak_width', 'geneID', 'geneID', 'score','gene_start','gene_end', 'strand')
                   )

## Update the column name.
## If it is in `old` and `new` is not in d, then will update it
.updateColnames<-function(d, maps=NULL, verbose=TRUE) {
  
  if (is.null(maps)) maps=COLMAPS
  
  for (i in 1:nrow(maps)) {
    if (maps[i, 1] %in% colnames(d)) {
      if (!(maps[i, 2] %in% colnames(d))) {
        colnames(d)[which(colnames(d)==maps[i, 1])]=maps[i, 2]
        if (verbose) cat(maps[i, 1], '-->', maps[i, 2], '\n')
      } else {
        if (verbose) cat(maps[i, 1], '-->', maps[i, 2], 'failed,', maps[i, 2], 'already in d\n')
      }
    }
  }
  return(d)
}


#' Add raw counts to QAPA results
#'
#' addQAPARawCounts Add raw counts from QAPA's quant.sf files to QAPA's APA list.
#'
#' @param qfile A file name of the QAPA APA list.
#' @param sfs A vector of quant.sf files of QAPA's outputs which record raw counts. Names of `sfs` should be sample names in the `qfile`.
#' @param suffix Default is 'CNT'. Will add columns named `sample.CNT` in the `qfile`.
#' @param ofile If not NULL, then output the data.frame to ofile.
#' @return A file name or a data.frame. The function will double check the ID and TPM columns to ensure that each QAPA row receives a unique CNT.
#' @export
#' @examples
#' \dontrun{
#' qfile='sim3_qapa_hg38.txt'
#' sfs=c('QAPA_saf/qapa_SRR847370_quant.sf', 'QAPA_saf/qapa_SRR847371_quant.sf',
#'       'QAPA_saf/qapa_SRR847374_quant.sf', 'QAPA_saf/qapa_SRR847375_quant.sf')
#' names(sfs)=paste0('sample', 1:4)
#' ofile='sim3_qapa_hg38_CNT.txt'
#' q=addQAPARawCounts(qfile, sfs, ofile=ofile)
#' }
#' @family APAdata functions
#' @export
addQAPARawCounts <- function(qfile, sfs, suffix = "CNT", ofile = NULL) {
  q = read.table(qfile, header = T, comment.char = "#", sep = "\t")
  head(q)

  if (!all(paste0(names(sfs), ".TPM") %in% colnames(q))) {
    stop("sample cols in QAPA file not the same as the sf files!\n")
  }

  q$id = paste0(q$UTR3.Start, "_", q$UTR3.End)

  .extractID <- function(longid) {
    ss = unlist(strsplit(longid, split = "_"))
    n = length(ss)
    return(paste0(ss[n - 1], "_", gsub("\\(.*", "", ss[n])))
  }

  for (i in 1:length(sfs)) {
    sf = read.table(sfs[i], header = TRUE)
    sf$id = unlist(lapply(sf$Name, .extractID))

    count_col = paste0(names(sfs)[i], ".", suffix)

    sf[, count_col] = sf$NumReads
    sf$SFTPM = sf$TPM
    sf = sf[, c("id", count_col, "SFTPM")]

    if (!all(q$id %in% sf$id)) {
      notIn = q$id[!(q$id %in% sf$id)]
      stop(
        length(notIn),
        " QAPA UTRStart_UTREnd not in the sf-file, for example: ",
        notIn[1],
        "\n"
      )
    }

    q = merge(q, sf, by.x = "id", by.y = "id")

    if (sum(q[, paste0(names(sfs)[i], ".TPM")] - q$SFTPM) != 0) {
      notIn = q$id[
        q[, paste0(names(sfs)[i], ".TPM")] - q$SFTPM != 0
      ]

      stop(
        length(notIn),
        " QAPA TPM not the same as the sf-file TPM, for example: ",
        notIn[1],
        "\n"
      )
    }

    q$SFTPM = NULL

    cat(
      ">>> process sf file",
      i,
      ", add",
      count_col,
      "\n"
    )
  }

  q$id = NULL

  if (!is.null(ofile)) {
    write.table(
      q,
      file = ofile,
      col.names = TRUE,
      row.names = FALSE,
      sep = "\t",
      quote = F
    )

    return(ofile)
  }

  return(q)
}

#' Convert a m6A data frame to PACdataset.
#'
#' m6A2PACds convert a m6A data frame to PACdataset.
#'
#' @param d A data frame with columns: seqnames, geneStrand, start, end, width, geneId, score, geneStart, geneEnd.
#' @return A PACdataset object. The coord is peak_start(+) or peak_end(-), add will add UPA_start, UPA_end columns.
#' @export
#' @family m6Adata functions
m6A2PACds <- function(d) {
  d=d[, c('seqnames','geneStrand','start','end', 'width', 'geneId', 'score','geneStart','geneEnd')]
  d=.updateColnames(d)
  
  d$strand[d$strand==1]='+'
  d$strand[d$strand==2]='-'
  d$coord=d$peak_start
  d$coord[d$strand=='-']=d$peak_end[d$strand=='-']
  
  d$peak_score=d$score
  
  #Some genes have peaks in two directions, such as PMS2P5
  #In this situation, movAPA::readPACds cannot pass, so we need to change the gene to another column name, such as geneID
  ug=unique(paste0(d$geneID, d$strand))
  if (length(ug)!=length(unique(d$gene))) cat('Caution: there are genes in both strands!\n')
  
  d$UPA_start=d$peak_start
  d$UPA_end=d$peak_end
  
  colDataFile=as.data.frame(matrix(c('group1'), ncol=1, dimnames=list('score', 'group')))
  d=readPACds(d, colDataFile)
  d@anno$gene=d@anno$geneID
  return(d)
}



#' Convert m6AExpress's data frame to PACdataset.
#'
#' m6AExpressDf2PACds converts m6AExpress's data frame to PACdataset.
#'
#' @param d A m6AExpress's data frame, with columns: seqnames, start, end, width, strand, gene_name, and <IP1, IP2, Treated_IP1, Treated_IP2..>.
#' @param smpCols If not NULL then will extract `smpCols`, otherwise will auto detect sample columns as remaining columns.
#' @return A PACdataset object. The coord is peak_start(+) or peak_end(-), add will add UPA_start, UPA_end columns.
#' @export
#' @examples
#' \dontrun{
#' m6ads=m6AExpressDf2PACds(peakDfNorm,
#'                    smpCols=c("IP1", "IP2", "Treated_IP1", "Treated_IP2"))
#' }
#' @family m6Adata functions
m6AExpressDf2PACds <- function(d, smpCols=NULL) {
  d=.updateColnames(d)

  # auto detect smpCols
  if (is.null(smpCols)) smpCols=colnames(d)[-which(colnames(d) %in% c(COLMAPS$new, COLMAPS$old))]
  
  if (!is.null(smpCols)) {
    if (!all(smpCols %in% colnames(d))) stop('not all smpCols in colnames(d)!')
  }
  
  d$peak_start=as.numeric(d$peak_start)
  d$peak_end=as.numeric(d$peak_end)  
  d$peak_width=as.numeric(d$peak_width)    
  d$coord=d$peak_start
  d$coord[d$strand=='-']=d$peak_end[d$strand=='-']

  #Some genes have peaks in two directions, such as PMS2P5
  #In this situation, movAPA::readPACds cannot pass, so we need to change the gene to another column name, such as geneID
  ug=unique(paste0(d$geneID, d$strand))
  if (length(ug)!=length(unique(d$gene))) cat('Caution: there are genes in both strands!\n')
  
  d$UPA_start=d$peak_start
  d$UPA_end=d$peak_end
  
  rid=grep('chr', d$chr)
  if (length(rid)==0) {
    cat('add chr to chrnames: 1-->chr1\n')
    d$chr=paste0('chr', d$chr)
  }
  
  colDataFile=as.data.frame(matrix(smpCols, ncol=1, dimnames=list(smpCols, 'group')))
  d=readPACds(d, colDataFile)
  d@anno$gene=d@anno$geneID
  return(d)
}

#' Convert m6A expression data to a PACdataset
#'
#' Perform optional differential m6A analysis, normalize IP/Input signals,
#' and convert the result into a PACdataset object.
#'
#' @param peakDf A data.frame containing m6A peak count data.
#' @param ctrls Character vector of control sample names.
#' @param treats Optional character vector of treated sample names.
#' @param libSizes Optional numeric vector of library sizes.
#' @param doDE Logical value indicating whether differential analysis
#' should be performed.
#' @param filterDE Logical value indicating whether to retain genes
#' with at least one differential m6A peak.
#' @param DEcutoff Significance cutoff for differential analysis.
#' @param DEtype Either \code{"pval"} or \code{"padj"}.
#'
#' @return A PACdataset object.
#'
#' @export
m6AExpress2PACds <- function(
    peakDf,
    ctrls,
    treats = NULL,
    libSizes = NULL,
    doDE = TRUE,
    filterDE = FALSE,
    DEcutoff = 0.05,
    DEtype = c("pval", "padj")) {
  DEtype <- match.arg(DEtype)

  if (!is.data.frame(peakDf)) {
    stop("peakDf must be a data.frame.")
  }

  if (!is.character(ctrls) || length(ctrls) == 0L) {
    stop(
      "ctrls must be a non-empty character vector."
    )
  }

  if (is.null(treats) || length(treats) == 0L) {
    treats <- character(0)

    if (isTRUE(doDE)) {
      warning(
        "No treated m6A data were provided. ",
        "Setting doDE = FALSE."
      )
      doDE <- FALSE
    }

    if (isTRUE(filterDE)) {
      warning(
        "filterDE requires treated m6A data. ",
        "Setting filterDE = FALSE."
      )
      filterDE <- FALSE
    }
  } else {
    if (!is.character(treats)) {
      stop(
        "treats must be NULL or a character vector."
      )
    }

    if (length(ctrls) != length(treats)) {
      stop(
        "ctrls and treats must have the same length."
      )
    }
  }

  if (doDE) {
    cat("do DE for peakDf...\n")

    peakDf <- m6AIpInputDE(
      peakDf,
      ctrls = ctrls,
      treats = treats,
      libSizes = libSizes,
      prefix = "DE_"
    )

    if (filterDE) {
      cat("filter DE from peakDf...\n")

      if (DEtype == "pval") {
        isDE <- !is.na(peakDf$DE_pvalue) &
          peakDf$DE_pvalue < DEcutoff
      } else {
        isDE <- !is.na(peakDf$DE_padj) &
          peakDf$DE_padj < DEcutoff
      }

      if (!"gene_name" %in% colnames(peakDf)) {
        stop(
          "peakDf must contain a gene_name column ",
          "when filterDE = TRUE."
        )
      }

      de_genes <- unique(peakDf$gene_name[isDE])

      cat(
        length(de_genes),
        "genes have >=1 DE m6A.\n"
      )

      peakDf <- peakDf[
        peakDf$gene_name %in% de_genes,
        ,
        drop = FALSE
      ]
    }
  }

  cat(
    "normalize m6A by ",
    "log2(IP/lib(IP))-log2(Input/lib(Input))\n",
    sep = ""
  )

  peakDf <- .m6AIpInputNormalize(
    peakDf,
    ctrls = ctrls,
    treats = treats,
    libSizes = libSizes,
    rawPrefix = "raw_"
  )

  cat("convert peakDf to PACdataset...\n")

  smpCols <- .getIPnames(ctrls)

  if (length(treats) > 0L) {
    smpCols <- c(
      smpCols,
      .getIPnames(treats)
    )
  }

  m6AExpressDf2PACds(
    peakDf,
    smpCols = smpCols
  )
}



#' Annotate a BED file and convert it to a PACdataset.
#'
#' bed2PACds annotates a BED file and converts it to a PACdataset.
#'
#' @param d A BED file name.
#' @param txdb a txdb object.
#' @return A PACdataset, with the BED file annotated by ChIPseeker::annotatePeak.
#' @export
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
#' d='DM/diff_treat_vs_control_c3.0_cond1.bed'
#' p1.txdb=bed2PACds(d, txdb) #5615 entrezID
#' p2.gff=bed2PACds(d, txdb) #5585 genesymbol
#' ## #The gene_start/end obtained from the two txdbs may have some differences,
#' ## but the annotations are similar
#' p1.txdb[p1.txdb@anno$coord==62550816]
#' p2.gff[p2.gff@anno$coord==62550816]
#' }
#' @family APAdata functions
bed2PACds <-function(d, txdb) {
  
  if (!inherits(txdb, "TxDb"))  stop("txdb should be a TxDb object, please use makeTxDbFromGFF to get TxDb first!")
  if (is.character(d)) {
    if (!grepl('bed$|BED$', d)) stop('d is a character, should be .bed file name!')
    d=import(d)
  }
  d <- ChIPseeker::annotatePeak(d,
                       tssRegion = c(-100,-1),
                       TxDb = txdb,
                       level='gene',
                       assignGenomicAnnotation = TRUE,
                       genomicAnnotationPriority = c("3UTR", "5UTR", "Exon", "Intron", "Promoter",
                                                     "Downstream", "Intergenic"),
                       overlap="all",
                       addFlankGeneInfo = FALSE)
  d=as.data.frame(d)
  d$strand[d$geneStrand==2]='-'
  d$strand[d$geneStrand==1]='+'
  d$strand=as.character(d$strand)
  dn=paste0(d[, c('seqnames','start','end')])
  isdup=duplicated(dn)
  if (sum(isdup)>0) {
    cat(sum(isdup), 'duplicated rows with multiple annotations, duplicated ones will be removed!')
    d=d[!isdup, ]
  }
  
  d$coord=d$start
  d$coord[d$strand=='-']=d$end[d$strand=='-']
  d=d[, c('seqnames','strand','coord','start','end','width','geneId','geneStart','geneEnd','annotation','score')]

  
  colnames(d)=c('chr','strand','coord','peak_start','peak_end','peak_width','gene','gene_start','gene_end','annoftr','score')
  d$peak_score=d$score #score to @counts
  d$UPA_start=d$peak_start
  d$UPA_end=d$peak_end
  
  colDataFile=as.data.frame(matrix(c('group1'), ncol=1, dimnames=list('score', 'group')))
  
  d=readPACds(d, colDataFile=colDataFile)
  return(d)
}

## ----------- APAds -------------

#' Get 3'UTR DEAPA from a 3UTRAPApd-like PACdataset.
#'
#' get3UTRAPApdDE gets 3'UTR DEAPA from a 3UTRAPApd-like PACdataset. This function filters DEAPA by fisher.test for each replicates of each pair of conditions.
#' As long as one replicate of a pair of conditions satisfies pval<cutoff, then it is marked as DEAPA.
#' This function determines the number and name of conditions, replicates based on `scol`, and consider rep1 and rep2 for each pair of condition pairs in order.
#' @param pacds A PACdataset object with `pdWhich` column which is obtained by `get3UTRAPApd`.
#' @param pthd P-value cutoff of Fisher.test to filter DEAPA, default is 0.05.
#' @param scol The column denoting "group" in pacds@colData, default is 'group'.
#' @param DEcol Default is 'isDE', Will add isDE=TRUE/FALSE to pacds@anno.
#' @param pvcol Default is 'DEPval', Will add DEPval to pacds@anno.
#' @param filterDE Default is TRUE to filter isDE=TRUE rows.
#' @return A PACdataset with DE information or filtered.
#' @export
#' @examples
#' \dontrun{
#' apadsUTR_m3=movAPA::get3UTRAPAds(apads_m3)
#' apadsPD_m3=get3UTRAPApd(apadsUTR_m3, minDist=50, maxDist=5000, minRatio=0.05,
#'                        fixDistal=FALSE, addCols='pd')
#' apadsUTR_m3DE=get3UTRAPApdDE(apadsPD_m3, pthd=0.05, filterDE=TRUE)
#' }
#' @family APAds functions
get3UTRAPApdDE <- function(
  pacds,
  pthd = 0.05,
  scol = 'group',
  DEcol = 'isDE',
  pvcol = 'DEPval',
  filterDE = TRUE
) {
  if (!('pdWhich' %in% colnames(pacds@anno))) {
    stop(
      "pdWhich column not in pacds, ",
      "please run get3UTRAPApd first to get PD-APA!\n"
    )
  }
  
  conds = as.character(unique(pacds@colData[, scol]))
  nrep = nrow(pacds@colData) / length(conds)
  
  if (nrow(pacds@colData) %% length(conds)) {
    stop("Number of reps are not the same for each condition!\n")
  }
  
  # Sample names for all replicates in each condition
  condReps = lapply(
    conds,
    function(par) {
      rownames(pacds@colData)[pacds@colData[, scol] == par]
    }
  )
  
  # Sample names for each replicate across all conditions.
  # For example, reps[[1]] = c("ctrlM141", "treatM141")
  # represents all conditions for the first replicate.
  reps = list()
  
  for (i in 1:nrep) {
    reps[[i]] = unlist(lapply(condReps, '[', i))
  }
  
  cat(
    'Test DE for',
    length(condReps),
    'conditions',
    length(reps),
    'reps\n'
  )
  
  isDE = rep(FALSE, nrow(pacds@counts))
  DEPval = rep(1, nrow(pacds@counts))
  
  for (i in seq(1, nrow(pacds@counts), by = 2)) {
    minpv = 1
    
    for (j in 1:length(reps)) {
      # Process each replicate
      pairs = combn(reps[[j]], 2)
      
      # Generate all pairs of conditions
      for (k in 1:ncol(pairs)) {
        # Construct a contingency table
        m = pacds@counts[c(i, i + 1), pairs[, k]]
        pv = fisher.test(m)$p.value
        minpv = min(minpv, pv)
      }
    }
    
    if (minpv < pthd) {
      isDE[c(i, i + 1)] = TRUE
    }
    
    DEPval[c(i, i + 1)] = minpv
  }
  
  invisible(gc())
  
  pacds@anno[, DEcol] = isDE
  pacds@anno[, pvcol] = DEPval
  
  cat(
    'Add two DE columns to pacds@anno:',
    DEcol,
    pvcol,
    '\n'
  )
  
  cat('Total APA-gene#', length(isDE) / 2, '\n')
  cat('DE APA-gene#', sum(isDE) / 2, '\n')
  
  if (filterDE) {
    pacds = pacds[isDE]
  }
  
  return(pacds)
}


#' Get hexamers of polyA signals.

#' @param grams the grams vector or be V1/MM/MOUSE/MM10.
#' @return grams of uppercase.
#' @export
#' @examples
#' getVarGrams('v1')
#' getVarGrams('MM')
#' getVarGrams(c('aataaa','ggattc'))
getVarGrams<-function(grams) {
  if (!is.null(grams)) {
    grams=toupper(grams)
    if (length(grams)==1 & grams[1]=='V1')
      grams=c('AATAAA','TATAAA','CATAAA','GATAAA','ATTAAA','ACTAAA','AGTAAA','AAAAAA','AACAAA','AAGAAA','AATTAA','AATCAA','AATGAA','AATATA','AATACA','AATAGA','AATAAT','AATAAC','AATAAG')
    if (length(grams)==1 & (grams[1] %in% c('MM','MOUSE','MM10'))) {
      grams=c('AATAAA','ATTAAA','TATAAA','AGTAAA','AATACA','CATAAA','AATATA','GATAAA','AATGAA','AATAAT','AAGAAA','ACTAAA','AATAGA','ATTACA','AACAAA','ATTATA','AACAAG','AATAAG')
    }
  }
  return(grams)
}


#' Add polyA signal (PAS) to a PACdataset.
#'
#' addPAS2PACds adds polyA signal (PAS) to a PACdataset. This function adds PAS_dist/gram/coord columns to PACds@anno.
#' The first position of coord is the first A in AATAAA; dist is the distance from the last nt (upstream) or first nt (downstream) of AATAAA to the polyA site.
#' The polyA site position is 0, and the extracted sequence is `from+to+1` nt.
#' Unlike movAPA's annotateByPAS, in this function, dist has positive and negative values, and will add the coord column, and will change the NA to NOPAS.
#'
#' @param pacds A PACdataset object.
#' @param bsgenome chrmosome fasta files, an object of BSgenome or FaFile, see faFromPACds().
#' @param grams a character vector to specify a gram like AATAAA, or v1 (AATAAA's variants), or multiple grams. grams can be not equal length, like c('AATAAA','AAATTT','CCCT')
#' @param from to specify the range near PACs, PAC is the 0 position.
#' e.g., from=-50, to=-1 to subset 50 nt (PAC is the 0 or 51st position, upstream 50nt of PAC), see faFromPACds().
#' @param to similar to from.
#' @param priority a numeric vector to set the priority and subgroups of grams if grams has multiple elements, default is NULL.
#'               For example, if grams=c('AATAAA','ATTAAA','AAAAAA','TTTAT), priority=(1,2,3,3), then will first search for AATAAA, if not exists, then for ATTAAA, then the remaining AAAAAA/TTTAT.
#'               If priority=NULL, then will treat all elements in grams as the same group (no priority).
#' @param label label A character string specifying the output column prefix. The function adds one or two columns, including `label_gram` and `dist`.
#' The _gram column gives the gram that is closest to the PAC.
#' _dist is the start position of a gram to a PAC.
#' @param chrCheck if TRUE, then all chr in PACds should be in bsgenome, otherwise will ignore those non-consistent chr rows in PACds.
#' @return  A PACdataset with columns such as `label_gram` and `dist` added.
#' @export
#' @examples
#' \dontrun{
#' paAPA.AATAAA=addPAS2PACds(pdAPA, bsgenome, grams='AATAAA',
#'                          from=-100, to=10, label='AATAAA_')
#' }
#' @family APAds functions
addPAS2PACds <- function (pacds, bsgenome, grams, from, to, priority=1, label='PAS_', chrCheck=TRUE) {
  
  grams=toupper(grams)
  
  if (length(grams)==1) {
    if (!is.null(priority)) 
      if (!identical(priority, 1)) 
        cat("grams has only 1 element, priority is not applicable!\n")
    priority=rep(1, length(grams))
  }
  
  if (is.null(priority) | identical(priority, 1)) priority=rep(1, length(grams))
  
  if (length(priority)!=length(grams)) stop("priority must be the same length as grams!\n")
  
  if (is.null(label)) label='PAS_'
  
  grams=getVarGrams(grams)
  
  if (!.isChrConsistent(pacds, bsgenome, allin = FALSE)) {
    cat('PACds chr not all in seqnames of bsgenome\n')
    if (chrCheck) {
      stop('Please check chr names!\n')
    } else {
      cat("chrCheck=FALSE, Remove rows in PACds with chr name not in bsgenome\n")
      pacds=subsetPACds(pacds, chrs = GenomeInfoDb::seqnames(bsgenome), verbose=TRUE)
    }
  }
  
  #seq=faFromPACds(pacds[1:5], bsgenome, what='updn', fapre=NULL, up=-5, dn=5)
  seq=faFromPACds(pacds, bsgenome, what='updn', fapre=NULL, up=from, dn=to)
  
  if (length(seq)!=nrow(pacds@counts)) {
    cat(length(seq),'fasta seqs, but',nrow(pacds@counts),'PACs, not matched!\n')
    stop("This is probably because that the coord/from/to of pacds is out of range of chromosomes\n")
  }
  
  dcol=paste0(label, 'dist'); ccol=paste0(label, 'coord'); gcol=paste0(label, 'gram')
  cat('Will add three columns to pacds@anno:',gcol, dcol, ccol,'\n')
  
  dd=list()
  for (gram in grams) {
    vm=vmatchPattern(pattern=gram, subject=seq, fixed=TRUE)
    #CCCTTCCTTTATAACTAGTGTCGCAACAATAAAATTTGAGCTTTGATCA
    # vm$start=28, with 27 characters before AATAAA
    
    nr=elementNROWS(vm)
    hasGram=nr
    hasGram[nr!=0]=1
    idx=which(hasGram!=0)
    svm=lapply(vm[idx], function(par) {start(par)})
    
    # Convert to a zero-based relative position:
    # positive values indicate downstream distances,
    # whereas negative values indicate upstream distances.
    svm2=lapply(svm, function(par) {
      p=par+from-1
      md=min(abs(p))
      return(p[which(abs(p)==md)][1])
      #min(abs(par+from-1))
    })
    
    # Example:
    # > seq$`AASDHPPT:PA5:chr11;+;105969405` = 10
    # The first A of AATAAA is at position 1, and there are 10 nt
    # from this position to the end of the sequence.
    # 50-letter DNAString object
    # seq: TGTAATGGTGATATGAAAAACTTTGTCTTGTCATTATAATAATAAAAAAA
    
    gramDist=rep(NA, length(hasGram))
    gramDist[idx]=unlist(svm2)
    
    gramCoord=rep(NA, length(hasGram))
    sidx=intersect(idx, which(pacds@anno$strand=='+'))
    if (length(sidx)>0) gramCoord[sidx]=pacds@anno$coord[sidx]+gramDist[sidx]
    sidx=intersect(idx, which(pacds@anno$strand=='-'))
    if (length(sidx)>0) gramCoord[sidx]=pacds@anno$coord[sidx]-gramDist[sidx] 
    
    if (length(grams)==1) {
      pacds@anno[, gcol]=rep(NA, length(pacds))
      pacds@anno[idx, gcol]=grams
      pacds@anno[, dcol]=gramDist
      pacds@anno[, ccol]=gramCoord
    } else {
      dd[[gram]]=gramDist
    }
  }
  
  if (length(grams)==1) {
    pacds@anno[is.na(pacds@anno[, gcol]), gcol]='NOPAS'
    return(pacds)
  }
  
  ## dd is a list, like dd$AATAAA=[NA, 5, NA, NA,... nseq]
  
  # If there are multiple grams, get the gram with the minimum
  # distance to the PAC.
  dd2=as.matrix(as.data.frame(dd)) # row=PAs, col=grams
  pacds@anno[, gcol]=NA
  pacds@anno[, dcol]=NA
  
  levels=sort(unique(priority))
  for (l in levels) {
    
    # Process the grams at the current priority level.
    dd=dd2[, priority==l, drop=F]
    
    # Get the gram with the minimum distance.
    mi=apply(dd, 1, function(par) {
      if (sum(!is.na(par))==0) return(NA)
      m=min(abs(par), na.rm=TRUE)
      mi=which(abs(par)==m)[1]
      return(mi)
    })
    
    # Get the minimum distance.
    mind=apply(dd, 1, function(par) {
      if (sum(!is.na(par))==0) return(NA)
      m=min(abs(par), na.rm=TRUE)
      return(par[which(abs(par)==m)][1])
    })
    
    minGram=mi
    minGram[!is.na(mi)]=grams[priority==l][mi[!is.na(mi)]]
    
    # If no higher-priority gram was found, fill in the PAS using
    # the current priority level.
    idNA=is.na(pacds@anno[, gcol])
    pacds@anno[idNA, gcol]=minGram[idNA]
    pacds@anno[idNA, dcol]=mind[idNA]
  }
  
  
  pacds@anno[, ccol]=rep(NA, length(pacds))
  sidx=which(pacds@anno$strand=='+' & !is.na(pacds@anno[, dcol]))
  if (length(sidx)>0) pacds@anno[, ccol][sidx]=pacds@anno$coord[sidx]+pacds@anno[, dcol][sidx]
  sidx=which(pacds@anno$strand=='-' & !is.na(pacds@anno[, dcol]))
  if (length(sidx)>0) pacds@anno[, ccol][sidx]=pacds@anno$coord[sidx]-pacds@anno[, dcol][sidx]
  
  pacds@anno[is.na(pacds@anno[, gcol]), gcol]='NOPAS'
  return(pacds)
}



#' Add animal's polyA signal (PAS) to a PACdataset.
#'
#' addMMPAS2PACds adds animal's polyA signal (PAS) to a PACdataset.
#' This function considers the priority of different PAS and search animal PAS, AATAAA 100nt -> ATTAAA 100nt -> TGTA 100nt -> others 50nt.
#'
#' @param pacds A PACdataset object.
#' @param bsgenome chrmosome fasta files, an object of BSgenome or FaFile, see faFromPACds().
#' @param from to specify the range near PACs, PAC is the 0 position.
#' e.g., from=-50, to=-1 to subset 50 nt (PAC is the 0 or 51st position, upstream 50nt of PAC), see faFromPACds().
#' @param to similar to from.
#' @param label A character string specifying the output column prefix. The function adds one or two columns, including `label_gram` and `dist`.
#' The _gram column gives the gram that is closest to the PAC.
#' _dist is the start position of a gram to a PAC.
#' @return  A PACdataset with columns such as `label_gram` and `dist` added.
#' @export
#' @examples
#' \dontrun{
#' pdAPA.mut=addMMPAS2PACds(pdAPA.mut, bsgenome,
#'                         from=-100, to=10, label='PAS_')
#' }
#' @family APAds functions
addMMPAS2PACds<- function(pacds, bsgenome, from=-100, to=10, label='PAS_') {
  grams=c('AATAAA','ATTAAA','TGTA')
  priority=1:3
  pasds1=addPAS2PACds(pacds, bsgenome, grams=grams, from, to, priority=priority, label=label, chrCheck=TRUE) 
  
  gcol=paste0(label,'gram')
  dcol=paste0(label,'dist')
  
  print(table(pasds1@anno[,gcol]))
  
  if (sum(is.na(pasds1@anno[, dcol]))==0) return(pasds1)
  
  # 取NA的部分继续找，但在from的一半范围内
  grams=getVarGrams('mm')[-(1:2)]
  p2=pasds1; naid=is.na(p2@anno[, dcol])
  p2@anno=p2@anno[naid, ]
  p2@counts=p2@counts[naid, ]
  pasds2=addPAS2PACds(p2, bsgenome, grams=grams, from=round(from/2), to=to, priority=NULL, label=label, chrCheck=TRUE) 
  
  print(table(pasds2@anno[, gcol]))
  
  anno2=pasds2@anno[!is.na(pasds2@anno[, dcol]),]
  anno1=pasds1@anno[!is.na(pasds1@anno[, dcol]), ]
  naAnno=pasds2@anno[is.na(pasds2@anno[, dcol]), ]
  
  rn=rownames(pacds@anno)
  anno=rbind(anno1,anno2,naAnno)
  pacds@anno=anno
  pacds@anno=pacds@anno[rn, ]
  pacds@counts=pacds@counts[rn, ]
  return(pacds)
}


#' Get proximal and distal polyA sites in a PACdataset
#'
#' get3UTRAPApd gets proximal and distal polyA sites in a PACdataset
#' Selection criteria: among all PA pairs that meet the criteria of `dist` and `ratio`, select the pair with the highest ratio sum and the closest distance.
#' Note that there is no filtering of PA expression levels here, which can be done in movAPA::subsetPACds.
#'
#' @param pacds A PACdataset object. By default, all are 3'UTR polyA sites; If there is an `ftr` column in PACds, then will check if all ftr are 3'UTR.
#' @param minDist Default is 50. Require the distance between  proximal and distal polyA site to be > minDist and < maxDist.
#' @param maxDist Default is 1000. Require the distance between  proximal and distal polyA site to be > minDist and < maxDist.
#' @param minRatio Default is 0.05. If >0 then will filter out polyA sites that do not meet the requirements based on the ratio first, and then filter proximal and distal sites.
#' @param fixDistal Default is FALSE. True to fix the farthest PA, that is, only look for the proximal site.
#' @param addCols Default is `pd`. If not NULL, then will add columns naming <addCols>Which/Score/Ratio/Dist, like pdWhich=P or D.
#' @return  A PACdataset with derived gram and distance columns added.
#' @export
#' @examples
#' \dontrun{
#'pacds=APAtrap2PACds(paCtrl)
#'pd1=get3UTRAPApd(pacds, minDist=50, maxDist=1000, minRatio=0.05, fixDistal=F)
#'pd2=get3UTRAPApd(pacds, minDist=50, maxDist=1000, minRatio=0.05, fixDistal=T)
#'pacds@anno$gene[!(pacds@anno$gene %in% pd1@anno$gene)]
#'movAPA::subsetPACds(pd1, gene='CTC1')
#'movAPA::subsetPACds(pd2, gene='CTC1')
#' }
#' @family APAds functions
get3UTRAPApd <- function(pacds, minDist=50, maxDist=1000, minRatio=0.05, fixDistal=FALSE, addCols='pd') {
  
  if ('ftr' %in% colnames(pacds@anno)) {
    ftrs=unique(pacds@anno$ftr)
    if (length(ftrs)!=1 | ftrs!='3UTR') stop("ftr is in pacds@anno, there are non-3UTR PACs, please get3UTRAPAds first!")
  }
  

  nc=rowSums(pacds@counts)
  annos=cbind(pacds@anno[, c('chr','strand','coord', 'gene')], nc=nc)
  annos$PA=rownames(annos)
  
  # Calculate the total count and usage ratio
  gnc=tapply(annos$nc, as.factor(annos$gene), sum)
  gnp=tapply(annos$gene, as.factor(annos$gene), length) 
  g=data.frame(gene=names(gnc), gnc, gnp)
  
  annos=merge(annos, g, by='gene')
  annos$ratio=annos$nc/annos$gnc
  
  if (minRatio>0) {
    ridx=(annos$ratio>minRatio)
    if (sum(ridx)>0) {
      annos=annos[ridx, ]
      gnp=tapply(annos$gene, as.factor(annos$gene), length) 
      gpass=names(gnp)[gnp>=2]
      cat('filtering by minRatio: gene# before:', length(gnp), '; after:', length(gpass),'; remove:', length(gnp)-length(gpass), '\n')
      annos=annos[annos$gene %in% gpass, ]
    }
  }
  
  # Sort by gene, chromosome, strand, and coordinate,
  # with the distal PA site placed first
  .sortAnnos <- function(annos) {
    anno1=annos[annos$strand=='+', ]
    anno1=anno1[order(anno1$gene, anno1$chr, -anno1$coord), ]
    anno2=annos[annos$strand=='-', ]
    anno2=anno2[order(anno2$gene, anno2$chr, anno2$coord), ]
    annos=rbind(anno1,anno2)
    return(annos)
  }
  
  annos=.sortAnnos(annos)
  
  ####
  .getAPApd <- function(anno) {
    
    # Fix the most distal PA site as the first site
    if (fixDistal) {
      anno$dist=abs(anno$coord-anno$coord[1])
      idx=which(anno$dist>0 & anno$dist<=maxDist & anno$dist>=minDist)
      if (length(idx)>0) {
        # Among PA sites satisfying the distance requirements,
        # select the site with the highest expression; in the event
        # of a tie, select the one closest to the distal PA site
        mid=idx[anno$nc[idx]==max(anno$nc[idx])][1]
        return(anno$PA[c(1, mid)])
      } else {
        return(c())
      }
    }
    
    # When fixDistal is FALSE, calculate the distance and
    # combined usage ratio for every pair of PA sites
    anno$ratio=anno$nc/sum(anno$nc)
    d=c()
    r=c()
    p1=c()
    p2=c()
    for (i in 1:(nrow(anno)-1)) {
      for (j in (i+1):(nrow(anno))) {
        p1=c(p1, anno$PA[i]); p2=c(p2, anno$PA[j])
        d=c(d, abs(anno$coord[i]-anno$coord[j]))
        r=c(r, anno$ratio[i]+anno$ratio[j])
      }
    }
    d=data.frame(p1=p1, p2=p2, dist=d, ratio=r)
    idx1=which(d$dist>0 & d$dist<=maxDist & d$dist>=minDist)
    idx2=which(d$ratio==max(d$ratio))
    idx=intersect(idx1, idx2)
    if (length(idx)==1) {
      return(c(d$p1[idx], d$p2[idx]))
    } else if (length(idx)==0) {
      return(c())
    } else {
      # If multiple PA-site pairs have the same combined usage ratio,
      # return the pair with the shortest distance
      idx3=idx[which(d$dist[idx]==min(d$dist[idx]))]
      return(c(d$p1[idx3], d$p2[idx3]))
    }
    return(c())
  }
  ####
  
  sp=split(annos, annos$gene, drop=TRUE)
  
  pas=unlist(lapply(sp, .getAPApd))
  cat('filtering pd (dist between): gene# before:', length(sp), '; after:', length(pas)/2,'; remove:', length(sp)-length(pas)/2, '\n')
  
  if (length(pas)==0) {
    pacds@counts=pacds@counts[-(1:nrow(pacds@counts)), ]
    pacds@anno=pacds@anno[-(1:nrow(pacds@anno)), ]
    return(pacds)
  } else {
    
    pacds@anno=pacds@anno[pas, ]
    pacds@anno=.sortAnnos(pacds@anno)
    pacds@counts=pacds@counts[rownames(pacds@anno), ]
    
    if (!is.null(addCols)) {
      newcols=paste0(addCols[1], c('Which','Score','Ratio','Dist'))
      pacds@anno[, newcols[1]]=rep(c('D','P'),length.out=nrow(pacds@anno))
      
      nc=rowSums(pacds@counts)
      idD=seq(1, nrow(pacds@anno), by=2)
      idP=idD+1
      ncg=rep(nc[idD]+nc[idP], each=2)
      
      pacds@anno[, newcols[2]]=nc
      pacds@anno[, newcols[3]]=nc/ncg
      
      pacds@anno[, newcols[4]]=rep(abs(pacds@anno$coord[idD]-pacds@anno$coord[idP]), each=2)
    }
    return(pacds)
  }
}


#' Add the nearest upstream m6A site to each polyadenylation site.
#'
#' `addM6A2PACds()` adds the nearest upstream m6A site to each
#' polyadenylation site in a PACdataset.
#'
#' The function searches for the nearest upstream m6A site of each
#' polyadenylation site and records the distance and m6A coordinate.
#' It adds columns with the prefix specified by `label` to
#' `pacds@anno`, including overlap and distance information.
#'
#' @param pacds A PACdataset object containing polyadenylation sites.
#' @param m6ads A PACdataset object containing m6A sites.
#' @param peakcols A character vector of length two specifying the column
#'   names in `m6ads` that contain the m6A start and end coordinates.
#'   The default is `c("peak_start", "peak_en")`.
#' @param label A character string specifying the prefix of the output
#'   columns. The default is `"m6A_"`. For example, the function may
#'   create columns such as `m6A_ovp`, `m6A_dist`, and columns derived
#'   from `mcols`.
#' @param d A numeric value specifying the upstream distance used to
#'   extend the polyadenylation site region when determining overlapping
#'   m6A sites. The default is 2000.
#' @param mcols A character vector of additional metadata column names
#'   to add to `pacds@anno`. These columns use the prefix specified by
#'   `label`.
#' @param verbose Logical value indicating whether progress messages
#'   should be shown. The default is `TRUE`.
#' @return A PACdataset containing the added m6A overlap, distance, and
#'   metadata columns.
#' @export
#' @examples
#' \dontrun{
#' p1 <- findOvpPACds(
#'   qryPACds = pdAPA,
#'   sbjPACds = m6ads,
#'   qryMode = "point",
#'   sbjMode = "region",
#'   d = 1000
#' )
#'
#' p2 <- annotateByKnownPAC(
#'   pdAPA,
#'   knownPACdss = m6ads,
#'   labels = "m6A",
#'   d = 1000,
#'   verbose = TRUE
#' )
#'
#' p3 <- addM6A2PACds(
#'   pdAPA,
#'   m6ads,
#'   peakcols = c("peak_start", "peak_end"),
#'   label = "m6A_",
#'   d = 1000,
#'   mcols = c("score", "coord", "gene")
#' )
#' }
#' @family APAds functions
addM6A2PACds <- function (pacds, m6ads, peakcols=c('peak_start','peak_end'), label='m6A_', d=2000, mcols=NULL, verbose=TRUE) {
  
  if (!is.null(mcols)) {
    if (!AinB(mcols, colnames(m6ads@anno), TRUE)) stop(mcols, 'not all in m6ads@anno!')
    if (!AinB(peakcols, mcols, TRUE)) mcols=c(peakcols, mcols)
  }
  
  if (!AinB(c(peakcols,'coord'), colnames(m6ads@anno), TRUE)) stop(c(peakcols, 'coord'), 'not all in m6ads@anno!')
  
  gr1 <- with(pacds@anno, GRanges(seqnames = chr,
                                  ranges =IRanges(start=coord,
                                                  end=coord),
                                  strand = strand) )
  
  # Extend upstream from each PA site and search only for m6A peaks
  # located upstream of the PA site
  gr1=resize(gr1, width=d+1, fix="end", use.names=TRUE, ignore.strand=FALSE)
  
  
  gr2 <- GRanges(seqnames = m6ads@anno$chr, ranges =IRanges(start=m6ads@anno[, peakcols[1]],
                                                            end=m6ads@anno[, peakcols[2]]),
                 strand = m6ads@anno$strand)
  
  ov = findOverlaps(gr1, gr2,
                    maxgap=0L, minoverlap=0L,
                    type=c("any"), select='first',
                    ignore.strand=FALSE)
  isOvp=ov
  im=which(!is.na(ov))
  isOvp[-im]=0
  isOvp[im]=1
  minD=rep(NA, length(ov))
  minD[im]=abs(m6ads@anno$coord[ov[im]]-pacds@anno$coord[im])
  
  pacds@anno[, paste0(label,'ovp')]=isOvp
  pacds@anno[, paste0(label,'dist')]=minD
  
  if (!is.null(mcols)) {
    pacds@anno[, paste0(label, mcols)]=NA
    pacds@anno[im, paste0(label, mcols)]=m6ads@anno[ov[im], mcols]
  }
  
  if (verbose) {
    cat('query:',length(gr1),'\n')
    cat('known:',length(gr2),'(',label,')\n')
    cat(length(im), 'query PACs overlapping m6A peaks\n')
  }
  
  return(pacds)
}



## Statistically analyze the occurrence of m6A peak and PAS in the proportional and fractional PAs in PACds
preStat <- function(pacds) {
  res=NULL
  for (i in c('P','D','fixed')) {
    if (i=='fixed') {
      p=pacds@anno[pacds@anno$pdWhich=='D', ]
      
      # An m6A site truly belongs to the distal PA site only when
      # m6A_coord (i.e., peak_start) is located between the proximal
      # and distal PA sites. Otherwise, reset the distal m6A signal.
      p$m6Actrl_coord[is.na(p$m6Actrl_coord)]=0
      pcoord=pacds@anno[pacds@anno$pdWhich=='P', 'coord']
      isbtw=p$m6Actrl_ovp & p$m6Actrl_coord<=p$coord & p$m6Actrl_coord>=pcoord
      p$m6Actrl_ovp[!isbtw]=0
    } else {
      p=pacds@anno[pacds@anno$pdWhich==i, ]
    }
    
    nPA=nrow(p)
    AATAAA=sum(p$PAS_gram=='AATAAA')
    TGTA=sum(p$PAS_gram=='TGTA')
    ATTAAA=sum(p$PAS_gram=='ATTAAA')
    noPAS=sum(p$PAS_gram=='NOPAS')
    v1PAS=nPA-AATAAA-TGTA-ATTAAA-noPAS
    m6Actrl=sum(p$m6Actrl_ovp==1)
    m6Amut=sum(p$m6Amut_ovp==1)
    m6Actrl_TGTA=sum(p$m6Actrl_ovp==1 & p$PAS_gram=='TGTA')
    m6Amut_TGTA=sum(p$m6Amut_ovp==1 & p$PAS_gram=='TGTA')
    m6Actrl_AATAAA=sum(p$m6Actrl_ovp==1 & p$PAS_gram=='AATAAA')
    m6Amut_AATAAA=sum(p$m6Amut_ovp==1 & p$PAS_gram=='AATAAA')
    m6Actrl_ATTAAA=sum(p$m6Actrl_ovp==1 & p$PAS_gram=='ATTAAA')
    m6Amut_ATTAAA=sum(p$m6Amut_ovp==1 & p$PAS_gram=='ATTAAA')
    stats=c(nPA=nPA, AATAAA=AATAAA, TGTA=TGTA, ATTAAA=ATTAAA, v1PAS=v1PAS, noPAS=noPAS, 
            m6Actrl=m6Actrl, m6Actrl_TGTA=m6Actrl_TGTA, m6Actrl_AATAAA=m6Actrl_AATAAA, m6Actrl_ATTAAA=m6Actrl_ATTAAA, 
            m6Amut=m6Amut, m6Amut_TGTA=m6Amut_TGTA, m6Amut_AATAAA=m6Amut_AATAAA, m6Amut_ATTAAA=m6Amut_ATTAAA)
    
    if (is.null(res)) {
      res=data.frame(stats)
    } else {
      res=cbind(res, stats)
    }
    colnames(res)[ncol(res)]=i
    rownames(res)=names(stats)
  }
  return(res)
}


## -------- m6A peaks -----
#' Read peak information from an RDA file
#'
#' @param rdafile Path to an RDA file.
#' @return A list containing peak information, library sizes, and sample names.
#' @export
.getRdaData<-function(rdafile) {
  load(rdafile)
  #head(Get_peak_infor)
  peakDf = Get_peak_infor[[1]]
  libSizes = Get_peak_infor[[2]]
  names(libSizes)=c('IP1', 'IP2','Treated_IP1','Treated_IP2','Input1','Input2','Treated_Input1','Treated_Input2')

  ctrls=c('IP1','IP2', 'Input1','Input2')
  treats=paste0('Treated_', ctrls)
  libSizes=libSizes[c(ctrls, treats)] #Reorder libsize by ctrls and treatments
  return(list(peakDf, ctrls, treats, libSizes))
}



#' Check m6A IP/Input sample names
#'
#' Checks whether sample names are arranged as IP samples followed by
#' Input samples.
#'
#' @param smpCols Character vector of sample column names. The first half
#'   should be IP samples and the second half should be Input samples.
#'
#' @return Invisibly returns `NULL`. Stops if the number of sample columns
#'   is odd.
#'
#' @export
## check ctrls/treats colnames, should be like "..IP.. ..IP.. ..Input.. ..Input.."
.checkIPInputNames<-function(smpCols) {
  
  if (length(smpCols)%%2) stop('smpCols should be half IP half Input!')
  
  ips=smpCols[1:(length(smpCols)/2)]
  if (sum(!grepl('IP', ips))>0) warnings('It seems that left half of smpCols is not IP columns\n')
  
  ips=smpCols[-(1:(length(smpCols)/2))]
  if (sum(!grepl('Input', ips))>0) warnings('It seems that right half of smpCols is not Input columns\n')
}

#' Get m6A IP sample names
#'
#' Extracts the first half of the sample columns, which should contain
#' m6A IP samples.
#'
#' @param smpCols Character vector with IP sample names first and Input
#'   sample names second.
#'
#' @return A character vector containing the IP sample names.
#'
#' @export
.getIPnames<-function(smpCols) {
  return(smpCols[1:(length(smpCols)/2)])
}

#' Get m6A Input sample names
#'
#' Extracts the second half of the sample columns, which should contain
#' Input samples.
#'
#' @param smpCols Character vector with IP sample names first and Input
#'   sample names second.
#'
#' @return A character vector containing the Input sample names.
#'
#' @export
.getInputNames<-function(smpCols) {
  return(smpCols[-(1:(length(smpCols)/2))])
}


#' Use QNBtest to test DE m6A peak
#'
#' m6AIpInputDE  use QNBtest to test DE m6A peak
#' This function searches the nearest m6A upstream of each polyA site and record the distance, coordinate of m6A.
#' This function will add to pacds@anno the m6A_ovp/dist.. columns. Dist is from peak_start of m6A to coord of the polyA site.
#'
#' @param peakDf A data frame of peaks.
#' @param ctrls cond1's IP and Input column names or idx.
#' @param treats cond2's IP and Input column names or idx.
#' @param libSizes Total read counts for each ctrls/treats (same order as ctrls/treats); if NULL then use the colSums.
#' @param prefix Default is `DE_`, which will add <DE_> to columns names like DE_log2.OR, DE_pvalue, DE_padj.
#' @return The given peak data frame with three new cols log2.OR, pvalue, padj.
#' @export
#' @examples
#' \dontrun{
#' dat=.getRdaData('m3mut_cond1.rda')
#' m6adsM3=m6AExpress2PACds(peakDf=dat[[1]], ctrls=dat[[2]],
#'                         treats=dat[[3]], libSizes=dat[[4]],
#'                        doDE=TRUE, filterDE=FALSE)
#' peakDf=m6AIpInputDE(peakDf, ctrls=ctrls, treats=treats,
#'                     libSizes=libSizes, prefix='DE_')
#' }
#' @family m6Ads functions
m6AIpInputDE <- function(peakDf, ctrls, treats, libSizes=NULL, prefix='DE_') {
  
  if (!all(c(ctrls, treats) %in% colnames(peakDf))) stop(toString(ctrls), 'or', toString(treats),'not all in colnames(peakDf)!')
  
  .checkIPInputNames(ctrls)
  .checkIPInputNames(treats)  
  
  libSizes <- colSums(peakDf[, c(ctrls, treats)])
  
  if (length(libSizes)!=length(c(ctrls, treats))) stop("libSizes not the same length as ctrls+treats!")
  
  sf=as.numeric(libSizes/exp(mean(log(libSizes))))
  
  ctrlIP <- peakDf[, .getIPnames(ctrls)]
  ctrlInput <- peakDf[, .getInputNames(ctrls)]
  treatIP <- peakDf[, .getIPnames(treats)]
  treatInput <- peakDf[,  .getInputNames(treats)]

  DEres <- QNBtest(ctrlIP, treatIP, ctrlInput, treatInput, mode="per-condition")
  DEres=DEres[,c(4,5,7)] #log2.OR, pvalue, padj
  colnames(DEres)=paste0(prefix, colnames(DEres))
  
  DEres=cbind(peakDf, DEres)
  pcol=which(colnames(DEres)==paste0(prefix, 'pvalue'))
  cat(sum(!is.na(DEres[, pcol]) & DEres[, pcol] <0.05), 'from', nrow(DEres),'m6A peaks is DE (pval<0.05)\n')
  pcol=which(colnames(DEres)==paste0(prefix, 'padj'))
  cat(sum(!is.na(DEres[, pcol]) & DEres[, pcol] <0.05), 'from', nrow(DEres),'m6A peaks is DE (padj<0.05)\n')
  
  invisible(gc())
  return(DEres)
}




#' Normalize m6A levels
#'
#' m6AIpInputNormalize gets normalized m6A levels.
#' After standardizing the IP and Input columns according to libsize, the IP/Input of each rep is then taken as log2.
#' In m6Aexpress, M6A level=log2 (Cip/Sip) - log2 (Cinput/Input), where Sip and Input are normalization factors; Cip is the m6A level, and input is the background level.

#' @param peakDf A data frame of peaks.
#' @param ctrls cond1's IP and Input column names or idx.
#' @param treats cond2's IP and Input column names or idx.
#' @param libSizes Total read counts for each ctrls/treats (same order as ctrls/treats); if NULL then use the colSums.
#' @param rawPrefix Default is `raw_`. If not NULL to keep raw counts of IP and Input, but adding 'raw_'.
#' @return The given peak data frame with raw_ctrls/treats + (IP column names denoting normalized columns).
#' @export
#' @examples
#' \dontrun{
#' peakDfNorm=m6AIpInputNormalize(peakDf, ctrls=ctrls, treats=treats,
#'                               libSizes=libSizes, rawPrefix='raw_')
#' }
#' @family m6Ads functions
m6AIpInputNormalize <-  function (peakDf, ctrls, treats, libSizes=NULL, rawPrefix='raw_') {
  
  if (!all(c(ctrls, treats) %in% colnames(peakDf))) stop(ctrls, 'or', treats,'not all in colnames(peakDf)!')
  
  .checkIPInputNames(ctrls)
  .checkIPInputNames(treats)  
  
  if (is.null(libSizes)) libSizes <- colSums(peakDf[, c(ctrls, treats)])
  
  if (length(libSizes)!=length(c(ctrls, treats))) stop("libSizes not the same length as ctrls+treats!")
  
  sf=as.numeric(libSizes/exp(mean(log(libSizes))))
  
  # 标准化
  peakDf[, c(ctrls, treats)] <- as.data.frame(t(t(peakDf[, c(ctrls, treats)])/sf))
  
  # log2(IP/Input)
  ctrlIP <- peakDf[, .getIPnames(ctrls)]
  ctrlInput <- peakDf[, .getInputNames(ctrls)]
  treatIP <- peakDf[, .getIPnames(treats)]
  treatInput <- peakDf[,  .getInputNames(treats)]
  
  ctrlMethy=log2((ctrlIP+0.01)/(ctrlInput+0.01))
  treatMethy=log2((treatIP+0.01)/(treatInput+0.01))
  
  colnames(ctrlMethy)=colnames(ctrlIP)
  colnames(treatMethy)=colnames(treatIP)
  
  ctrlMethy[ctrlMethy<0]=0.001
  treatMethy[treatMethy<0]=0.001
  
  if (!is.null(rawPrefix)) {
    newcols=paste0(rawPrefix, c(ctrls, treats))
    colnames(peakDf)[match( c(ctrls, treats), colnames(peakDf))]=newcols
  } else {
    peakDf[, c(ctrls, treats)]=NULL
  }
  
  peakDf=cbind(peakDf, ctrlMethy, treatMethy)
  return(peakDf)
}



## -------- per gene -----
#Convert a PACds to dataframe [chr,strand, UPA_start/end, samples...], rownames=rownames(PACds@counts)
PACdsRanges2df <- function(PACds) {
  if (!AinB(c('chr','strand','UPA_start','UPA_end'), colnames(PACds@anno))) stop('chr/strand/UPA_start/UPA_end not in PACds@anno')
  df=cbind(PACds@anno[,c('chr','strand','UPA_start','UPA_end')], PACds@counts)
  return(df)
}


#Convert a PACds to dataframe [chr,strand,coord, samples...], rownames=rownames(PACds@counts)
PACds2PAdf <- function(PACds) {
  if (!AinB(c('chr','strand','coord'), colnames(PACds@anno))) stop('chr/strand/coord not in PACds@anno')
  df=cbind(PACds@anno[,c('chr','strand','coord')], PACds@counts)
  return(df)
}


#' Merge multiple PACdataset objects by cleavage-site coordinates
#'
#' mergePACdsCoords groups nearby poly(A) sites from one or multiple
#' `PACdataset` objects. This function was modified from
#' `movAPA::mergePACds` and aggregates pAs according to their `coord` values.
#'
#' This function is useful for grouping nearby cleavage sites into PACs. It can
#' also merge multiple PA or PAC datasets into a single `PACdataset` for
#' differential APA or other downstream analyses. After grouping or merging,
#' the returned object can be annotated using `annotatePAC()`.
#'
#' @param PACdsList A `PACdataset` object or a list of `PACdataset` objects.
#'   Each `PACds@anno` must contain `chr`, `strand`, and `coord` columns. If
#'   `colData` is absent, sample groups are assigned as `group1`, `group2`, and
#'   so on.
#' @param d Distance used to group nearby PACs. Default is 24 nt.
#'
#' @return A merged `PACdataset` object. The `counts` slot stores counts of
#'   merged samples. The `colData` slot stores group annotations derived from
#'   the first column of each input object's `colData`. The `anno` slot includes
#'   `chr`, `strand`, `coord`, `tottag`, `UPA_start`, `UPA_end`, `nPA`, and
#'   `maxtag`.
#'
#' @examples
#' \dontrun{
#' # Group nearby pAs into PACs.
#' data(PACds)
#' PACds@counts <- rbind(PACds@counts, PACds@counts)
#' PACds@anno <- rbind(PACds@anno, PACds@anno)
#' ds <- mergePACdsCoords(PACds, d = 24)
#'
#' # Merge two PACdataset objects.
#' pacds <- new("PACdataset", counts = PACds@counts, anno = PACds@anno)
#' PACdsList <- list(PACds, pacds)
#' ds <- mergePACdsCoords(PACdsList, d = 24)
#' }
#'
#' @family APAds functions
#' @export
mergePACdsCoords <- function(PACdsList, d = 24) {
  if (inherits(PACdsList, "PACdataset")) {
    PACdsList = list(PACdsList)
  }

  for (i in seq_along(PACdsList)) {
    if (!AinB(c("chr", "strand", "coord"),
              colnames(PACdsList[[i]]@anno))) {
      stop("chr/strand/coord not in PACdsList@anno")
    }

    if (nrow(PACdsList[[i]]@colData) == 0) {
      PACdsList[[i]]@colData = as.data.frame(matrix(
        rep(paste0("group", i), ncol(PACdsList[[i]]@counts)),
        ncol = 1,
        dimnames = list(colnames(PACdsList[[i]]@counts), "group")
      ))
    }
  }

  allpa = PACds2PAdf(PACdsList[[1]])

  if (length(PACdsList) >= 2) {
    for (j in 2:length(PACdsList)) {
      pa2 = PACds2PAdf(PACdsList[[j]])

      allpa = merge(
        allpa,
        pa2,
        all = TRUE,
        by.x = c("chr", "strand", "coord"),
        by.y = c("chr", "strand", "coord")
      )
    }
  }

  allpa[is.na(allpa)] = 0
  invisible(gc())

  gr <- with(
    allpa,
    GRanges(
      seqnames = chr,
      ranges = IRanges(start = coord, end = coord),
      strand = strand
    )
  )

  gr = GenomicRanges::resize(
    gr,
    width = d + 1,
    fix = "start",
    use.names = TRUE,
    ignore.strand = FALSE
  )

  itv = GenomicRanges::reduce(gr, drop.empty.ranges = TRUE)

  cat("Group PA to PACs\n")

  ov = GenomicRanges::findOverlaps(
    gr,
    itv,
    maxgap = -1L,
    minoverlap = 1L,
    type = "any",
    select = "all",
    ignore.strand = FALSE
  )

  ov = as.data.frame(ov)

  allpa$idx = seq_len(nrow(allpa))
  allpa = merge(allpa, ov, by.x = "idx", by.y = "queryHits")
  allpa$idx = NULL

  smpcols = colnames(allpa)[!(colnames(allpa) %in% c(
    "chr",
    "strand",
    "coord",
    "subjectHits"
  ))]

  cat("count tot tag for each sample within each PAC\n")

  allpa$tottag = rowSums(allpa[, smpcols, drop = FALSE])

  byItv <- dplyr::group_by(allpa, subjectHits)

  dots <- sapply(
    smpcols,
    function(x) substitute(sum(x), list(x = as.name(x)))
  )

  dots[["tottag"]] = substitute(sum(x), list(x = as.name("tottag")))

  pac = do.call(dplyr::summarise, c(list(.data = byItv), dots))

  cat("Annotate the range of each PAC\n")

  pacAnno = byItv %>% dplyr::summarise(
    nPA = n(),
    UPA_start = min(coord),
    UPA_end = max(coord),
    maxtag = max(tottag),
    coord = coord[which.max(tottag)],
    chr = chr[1],
    strand = strand[1]
  )

  pac = merge(pac, pacAnno, by = "subjectHits")

  smpcols = smpcols[smpcols != "tottag"]

  # Ensure that the PACdataset count matrix is numeric.
  counts_df = pac[, smpcols, drop = FALSE]

  counts = as.matrix(data.frame(
    lapply(counts_df, function(x) as.numeric(as.character(x))),
    check.names = FALSE
  ))

  rownames(counts) = rownames(pac)

  anno = pac[, c(
    "chr",
    "strand",
    "coord",
    "tottag",
    "UPA_start",
    "UPA_end",
    "nPA",
    "maxtag"
  ), drop = FALSE]

  colData = as.data.frame(matrix(
    unlist(lapply(
      PACdsList,
      function(ds) as.character(ds@colData[, 1])
    )),
    ncol = 1,
    dimnames = list(colnames(counts), "group")
  ))

  new(
    "PACdataset",
    counts = counts,
    colData = colData,
    anno = anno
  )
}




#' Merge multiple PACdataset
#'
#' mergePACdsRanges groups nearby PACs from single/multiple PACdataset objects.
#' This function was modified from movAPA::mergePACds, aggregating pAs according
#' to genomic ranges. Peaks are merged according to `UPA_start` and `UPA_end`,
#' where `d` is the allowed gap. The PA ID from each input PACdataset is recorded
#' for each merged peak.
#'
#' @param PACdsList A `PACdataset` object or a list of multiple `PACdataset`
#'   objects. `PACds@anno` must contain the `chr`, `strand`, `UPA_start`, and
#'   `UPA_end` columns. If names are provided, the corresponding PA IDs are
#'   recorded using these names; otherwise, names are assigned as `ds1`, `ds2`,
#'   and so on.
#' @param d Allowed gap between peaks in PACds. If `d = 0`, two peaks must
#'   overlap by at least 1 nt to be merged. Default is 24 nt.
#'
#' @return A merged `PACdataset` object with annotation columns including
#'   `chr`, `strand`, `tottag`, `UPA_start`, `UPA_end`, `nPA`, input-dataset
#'   PA IDs, and `coord`.
#'
#' @examples
#' \dontrun{
#' # Merge two PACdataset objects
#' pacds <- new("PACdataset", counts = PACds@counts, anno = PACds@anno)
#' PACdsList <- list(PACds, pacds)
#' ds <- mergePACdsRanges(PACdsList, d = 24)
#' }
#'
#' @family APAds functions
#' @export
mergePACdsRanges <- function(PACdsList, d = 0) {
  if (inherits(PACdsList, "PACdataset")) {
    PACdsList = list(PACdsList)
  }

  for (i in seq_along(PACdsList)) {
    if (!all(c("chr", "strand", "UPA_start", "UPA_end") %in%
             colnames(PACdsList[[i]]@anno))) {
      stop("chr/strand/UPA_start/UPA_end not in PACdsList@anno")
    }

    if (nrow(PACdsList[[i]]@colData) == 0) {
      PACdsList[[i]]@colData = as.data.frame(matrix(
        rep(paste0("group", i), ncol(PACdsList[[i]]@counts)),
        ncol = 1,
        dimnames = list(colnames(PACdsList[[i]]@counts), "group")
      ))
    }
  }

  if (is.null(names(PACdsList))) {
    names(PACdsList) = paste0("ds", seq_along(PACdsList))
  }

  allpa = PACdsRanges2df(PACdsList[[1]])
  allpa[, names(PACdsList)[1]] = rownames(PACdsList[[1]]@anno)

  if (length(PACdsList) >= 2) {
    for (j in 2:length(PACdsList)) {
      pa2 = PACdsRanges2df(PACdsList[[j]])
      pa2[, names(PACdsList)[j]] = rownames(PACdsList[[j]]@anno)

      allpa = merge(
        allpa,
        pa2,
        all = TRUE,
        by.x = c("chr", "strand", "UPA_start", "UPA_end"),
        by.y = c("chr", "strand", "UPA_start", "UPA_end")
      )
    }
  }

  allpa[is.na(allpa)] = 0
  invisible(gc())

  gr <- with(
    allpa,
    GRanges(
      seqnames = chr,
      ranges = IRanges(start = UPA_start, end = UPA_end),
      strand = strand
    )
  )

  mcols(gr) = allpa[, names(PACdsList), drop = FALSE]

  # Expand 3' end.
  if (d > 0) {
    gr = resize(
      gr,
      width = width(gr) + d,
      fix = "start",
      use.names = TRUE,
      ignore.strand = FALSE
    )
  }

  itv = GenomicRanges::reduce(gr, drop.empty.ranges = TRUE)

  cat("Group ranges\n")

  ov = GenomicRanges::findOverlaps(
    gr,
    itv,
    maxgap = -1L,
    minoverlap = 1L,
    type = "any",
    select = "all",
    ignore.strand = FALSE
  )

  ov = as.data.frame(ov)

  allpa$idx = seq_len(nrow(allpa))
  allpa = merge(allpa, ov, by.x = "idx", by.y = "queryHits")
  allpa$idx = NULL

  smpcols = colnames(allpa)[!(colnames(allpa) %in% c(
    "chr",
    "strand",
    "UPA_start",
    "UPA_end",
    "subjectHits",
    names(PACdsList)
  ))]

  cat("count tot tag for each sample within each range\n")

  allpa$tottag = rowSums(allpa[, smpcols, drop = FALSE])

  byItv <- dplyr::group_by(allpa, subjectHits)

  dots <- sapply(
    smpcols,
    function(x) substitute(sum(x), list(x = as.name(x)))
  )

  dots[["tottag"]] = substitute(sum(x), list(x = as.name("tottag")))

  pac = do.call(dplyr::summarise, c(list(.data = byItv), dots))

  pacAnno = byItv %>% dplyr::summarise(
    nPA = n(),
    UPA_start = min(UPA_start),
    UPA_end = max(UPA_end),
    chr = chr[1],
    strand = strand[1]
  )

  pac = merge(pac, pacAnno, by.x = "subjectHits", by.y = "subjectHits")

  smpcols = smpcols[smpcols != "tottag"]

  # Record the PA ID corresponding to each input PACdataset.
  paids = unique(allpa[, c("subjectHits", names(PACdsList))])

  .tostr <- function(par) {
    return(toString(par[par != 0]))
  }

  paids = dplyr::group_by(paids, subjectHits) %>%
    dplyr::summarise(dplyr::across(everything(), .tostr))

  pac = merge(pac, paids, by.x = "subjectHits", by.y = "subjectHits")

  # Ensure that the PACdataset count matrix is numeric.
  counts_df = pac[, smpcols, drop = FALSE]

  counts = as.matrix(data.frame(
    lapply(counts_df, function(x) as.numeric(as.character(x))),
    check.names = FALSE
  ))

  rownames(counts) = rownames(pac)

  anno = pac[, c(
    "chr",
    "strand",
    "tottag",
    "UPA_start",
    "UPA_end",
    "nPA",
    names(PACdsList)
  ), drop = FALSE]

  anno$coord = anno$UPA_start
  anno$coord[anno$strand == "-"] = anno$UPA_end[anno$strand == "-"]

  colData = as.data.frame(matrix(
    unlist(lapply(
      PACdsList,
      function(ds) as.character(ds@colData[, 1])
    )),
    ncol = 1,
    dimnames = list(colnames(counts), "group")
  ))

  ds = new(
    "PACdataset",
    counts = counts,
    colData = colData,
    anno = anno
  )

  return(ds)
}



#' Calculate gene level m6A score
#'
#' This function first searches the closest polyA site within the same gene for each m6A site,
#' and calculate the m6A score by distance between polyA site and m6A using M=exp (- d/d0).
#' d0 is the quantile of 75% of all nearest m6A and polyA site distance.
#' Then the gene level m6A score is the average m6A level of all m6As in a gene after distance weighting.
#'
#' @param pacds A PACdataset object of polyA sites.
#' Both m6ads and pacds are required to have a gene column.
#' We can first annotate both `ds` with 3' UTR extension (ensuring m6A and polyA sites are all within the same gene), and then call this function.
#' @param m6Ads A PACdataset object of m6A sites.
#' @param sufix Suffix for the columns names of m6A score, default is '_m6A'.
#' @return A data frame with columns like [gene, d0, nM6A, ...samples in m6Ads@counts with sufix _m6A].
#' @export
#' @examples
#' m6Ads=movAPA::makeExamplePACds()
#' m6Ads@anno$gene=sample(LETTERS, length(m6Ads), replace=TRUE)
#' pacds=movAPA::makeExamplePACds()
#' pacds@anno$gene=sample(LETTERS, length(pacds), replace=TRUE)
#' m=getM6AperGene(m6Ads, pacds)
#' head(m)
#' @family QBGLM functions
getM6AperGene <- function(m6Ads,
                          pacds,
                          sufix='_m6A') {

  if (!('gene' %in% colnames(m6Ads@anno)) | !('gene' %in% colnames(pacds@anno))) stop('getM6AperGene error: gene column not in @anno of m6Ads or pacds!\n')

  gm=unique(m6Ads@anno$gene)
  gp=unique(pacds@anno$gene)
  govp=intersect(gm, gp)
  if (length(govp)<length(gm)/10) stop('getM6AperGene error: gene name style is not the same between m6Ads and pacds!\n')

  if (length(govp) != length(gm))
    cat(length(govp), 'out of', length(gm), 'm6A genes have PACs, and will aggregate m6Ads@counts for these genes\n')

  m6Ads=movAPA::subsetPACds(m6Ads, genes=govp)

  pds=as.data.frame(matrix(nrow=0,ncol=3, dimnames = list(NULL, c('gene','mid','dist'))))

  for (g in govp) {
    mid=which(m6Ads@anno$gene==g)
    mc=m6Ads@anno$coord[mid]
    pc=pacds@anno$coord[pacds@anno$gene==g]

    # min dist to all PACs for each m6A site
    pd=lapply(mc, function(par) abs(par-pc))
    pd=unlist(lapply(pd, min))

    pds=rbind(pds, data.frame(gene=g, mid=rownames(m6Ads@anno)[mid], dist=pd))
  }

  d0 <- max(1, round(quantile( pds$dist, 0.75)))
  # decay: exp(-d/d0)
  pds$decay=exp(-(pds$dist)/d0)

  # each m6A score * decay
  counts=m6Ads@counts[pds$mid, ] #in order
  counts=counts*pds$decay
  #counts$gene=pds$gene

  # merge per gene
  sp=base::split(as.data.frame(counts), factor(pds$gene), drop=F)

  sp=lapply(sp, function(par) c(colMeans(par), nrow(par)))
  counts= t(as.data.frame(sp))
  colnames(counts)[-ncol(counts)]=paste0(colnames(counts)[-ncol(counts)], sufix)
  colnames(counts)[ncol(counts)]='nM6A'
  counts=as.data.frame(counts)
  counts$d0=d0
  counts$gene=names(sp)
  nc=ncol(counts)
  counts=counts[, c(nc, nc-1, nc-2, 1:(nc-3))]
  return(counts)
}

#' Calculate gene level RUD score
#'
#' This function calculates gene level RUD score.
#'
#' @param pacdsPD A PACdataset object with only proximal and distal polyA sites, which could be obtained by `get3UTRAPApd`.
#' @param sufix1 Suffix for the columns names of RUD score, default is '_RUD'.
#' @param sufix2 Suffix for the columns names of counts, default is '_weight'.
#' @return A data frame with columns like [pacdsPD@counts with sufix _RUD and _weight].
#' @export
#' @examples
#' \dontrun{
#' apadsUTR=movAPA::get3UTRAPAds(apads)
#' apadsPD=get3UTRAPApd(apadsUTR, minDist=50, maxDist=5000, minRatio=0.05,
#'                     fixDistal=FALSE, addCols='pd')
#' p=getRUDperGene(apadsPD)
#' }
#' @family QBGLM functions
getRUDperGene<-function(pacdsPD,
                        sufix1='_RUD',
                        sufix2='_weight') {

  if (!AinB(c('pdWhich','gene'), colnames(pacdsPD@anno))) stop("getRUDperGene error: please run get3UTRAPApd first!")

  #pacdsPD@counts=pacdsPD@counts+psudo
  #Do not add psudo values to PA with 0 expression.
  #If the gene is not expressed, then RUD is NaN;
  #Adding psudo in this way is also incorrect, as PA that is not expressed will become a ratio of 0.5!
  did=which(pacdsPD@anno$pdWhich=='D')
  pid=which(pacdsPD@anno$pdWhich=='P')

  rud=pacdsPD@counts[did, ]/(pacdsPD@counts[pid, ]+pacdsPD@counts[did, ])
  colnames(rud)=paste0(colnames(rud),sufix1)

  gene=pacdsPD@anno$gene[did]
  weights=round(pacdsPD@counts[pid, ]+pacdsPD@counts[did, ])
  colnames(weights)=paste0(colnames(weights),sufix2)

  res=as.data.frame(cbind(gene, rud, weights))
  res[, -1]=lapply(res[, -1], as.numeric)
  return(res)
}

#' Get highly variable genes by m6A and RUD scores
#'
#' This function calculates the CV of columns m6A and RUD, and filters highly variable genes.
#'
#' @param perGeneData A data frame of per gene data with m6A$ and RUD$ columns, which can be obtained by getRUDperGene and getM6AperGene.
#' @param m6Acutoff Cutoff for m6A score; if >0, then do filtering.
#' @param RUDcutoff Cutoff for RUD score; if >0, then do filtering.
#' @param totalWeights Cutoff for weights (total counts); if >0, then do filtering.
#' @param perWeights Cutoff for weights per sample; if >0, then do filtering. Default is 10.
#' @param nWeights Number of samples with weights>=perWeights, default is 6.
#' @param minAvgM6A Cutoff for per sample m6A score; if >0, then do filtering.
#' @return A filtered data frame.
#' @export
#' @examples
#' \dontrun{
#' mp1=filterHVPerGene(mp, m6Acutoff=0.3, RUDcutoff=0.3,
#'                     totalWeights=100, perWeights=10, nWeights=6)
#' }
filterHVPerGene<-function(perGeneData,
                          m6Acutoff=-1,
                          RUDcutoff=-1,
                          totalWeights=-1,
                          perWeights=10,
                          nWeights=6,
                          minAvgM6A=0.1) {

  keep1=rep(T, nrow(perGeneData))
  if (m6Acutoff>0) {
    cid=grep('m6A$', colnames(perGeneData))
    m=rowMeans(perGeneData[, cid])
    s=apply(perGeneData[, cid], 1, sd)
    cv=s/m
    keep1=(cv>m6Acutoff)
    cat(sum(keep1), 'm6A HVG\n')
  }

  keep2=rep(T, nrow(perGeneData))
  if (RUDcutoff>0) {
    cid=grep('RUD$', colnames(perGeneData))
    m=rowMeans(perGeneData[, cid])
    s=apply(perGeneData[, cid], 1, sd)
    cv=s/m
    keep2=(cv>RUDcutoff)
    cat(sum(keep2), 'RUD HVG\n')
  }

  keep3=rep(T, nrow(perGeneData))
  if (totalWeights>0 & perWeights>0 & nWeights>0) {
    cid=grep('weight', colnames(perGeneData))
    keep3=rowSums(perGeneData[, cid])>totalWeights & sum(perGeneData[, cid]>=perWeights)>=nWeights
    cat(sum(keep3), 'pass weights filtering\n')
  }

  keep4=rep(T, nrow(perGeneData))
  if (minAvgM6A>0) {
    cid=grep('m6A$', colnames(perGeneData))
    keep4=(rowMeans(abs(perGeneData[, cid]))>minAvgM6A)
    cat(sum(keep4), 'pass minAvgM6A filtering\n')
  }

  keep=(keep1 & keep2 & keep3 & keep4)
  cat(sum(keep), 'RUD&m6A&weights&minAvgM6A\n')
  return(perGeneData[keep, ])
}

## ------ GLM -----------
.FitModel <- setClass("FitModel",
                      slots = c(
                        type = "character", #glm, fitError
                        gene="character",
                        df.residual = "numeric",
                        dispersion = "numeric",
                        coefficients = "numeric",
                        raw.pvals= "numeric", #pvals without squeezeVar
                        cov.unscaled = "matrix",
                        var.post = "numeric",
                        df.post = "numeric"
                      )
)


#' @name FitModel
#' @title FitModel
#'
#' @description Function for contstructing a new `FitModel` object.
#'
#' @param type default set to fitError, can be a glm
#' @param gene gene
#' @param df.residual Numeric, residual of the fitted glm
#' @param dispersion Numeric, dispersion of the fitted glm
#' @param coefficients Numeric, coefficients of the fitted glm
#' @param raw.pvals Numeric, p-values of the fitted glm
#' @param cov.unscaled Numeric, unscaled variance of the fitted glm
#' @param var.post Numeric, posterior variance of the glm, default is NA
#' @param df.post Numeric, posterior degrees of freedom of the glm, default is NA
#' @return A FitModel object
#' @importFrom methods new
#' @export
#' @family QBGLM functions
FitModel <- function(type = "fitError",
                     gene = NA_character_,
                     df.residual = NA_real_,
                     dispersion = NA_real_,
                     coefficients = NA_real_,
                     raw.pvals = NA_real_,
                     cov.unscaled = matrix(NA,0,0),
                     var.post = NA_real_,
                     df.post = NA_real_) {
  out <- new("FitModel")
  out@type <- type
  out@gene <- gene
  out@df.residual <- df.residual
  out@dispersion <- dispersion
  out@coefficients <- coefficients
  out@raw.pvals <- raw.pvals
  out@cov.unscaled <- cov.unscaled
  out@var.post <- var.post
  out@df.post <- df.post
  return(out)
}



#' Remove fitError models in a model list
#'
#' This function removes fitError models in a model list.
#'
#' @param models A model list after fitting QBGLM.
#' @param verbose TRUE to show message.
#' @return A list with fitError model removed.
#' @export
#' @family QBGLM functions
removeFitErrorModels<-function(models, verbose=TRUE) {
  eid=which(unlist(lapply(models, function(m){m@type}))=='fitError')
  if (length(eid)>0) {
    if (verbose) cat("remove", length(eid), "fitError models\n")
    models=models[-eid]
  }
  return(models)
}



#' Fit QBGLM for genes
#'
#' This function fits QBGLM for a data frame of genes with m6A+RUD+Weights, and returns a model list with each element a fitted GLM for a gene.
#'
#' @param mpMat A data frame of genes with m6A+RUD+Weights, there should be columns ending with "_m6A", "_RUD", "_weight".
#' @param noFitError TRUE to remove fitError models.
#' @return A list of fitted models with each element a fitted GLM for a gene.
#' @export
#' @examples
#' \dontrun{
#' p=getRUDperGene(apadsPD)
#' m=getM6AperGene(m6ads, apadsPD)
#' mp=merge(m, p, by.x='gene', by.y='gene')
#' models=fitQuasiGLM(mp)
#' table(unlist(lapply(models, function(par) par@type)))
#' }
#' @family QBGLM functions
fitQuasiGLM<-function(mpMat,
                      noFitError=TRUE) {

  # m6A, RUD, weights
  mid=grep('_m6A', colnames(mpMat))
  rid=grep('_RUD', colnames(mpMat))
  wid=grep('_weight', colnames(mpMat))

  if (length(mid)!=length(rid) | length(rid) !=length(wid))
    stop("Number of columns with _m6A, _RUD, _weight not equal!")

  .fitOne <- function(aRow) {

    if (any(is.na(aRow))) {
      .out <- FitModel(type = "fitError")
      return(.out)
    }

    n=length(aRow)/3
    d=data.frame(m6A=aRow[1:n], rud=aRow[(n+1):(2*n)], weights=aRow[(2*n+1):(3*n)])

    model <- try(glm(rud ~ m6A, family = "quasibinomial", data=d, weights=d$weights))
    #summary(model)

    if (inherits(model, "try-error")) {
      .out <- FitModel(type = "fitError")
      return(.out)
    } else {
      if (any(is.na(model$coefficients))) {
        .out <- FitModel(type = "fitError")
        return(.out)
      }

      type <- "glm"
      .out <- .FitModel(
        type = type,
        coefficients=model$coefficients,
        raw.pvals=summary(model)$coefficients[, 4],
        df.residual=model$df.residual,
        dispersion=summary(model)$dispersion,
        cov.unscaled=summary(model)$cov.unscaled,
        var.post = as.numeric(NA),
        df.post = as.numeric(NA)
      )

      # The result may have NA, which also considered as FitError
      if (is.na(.out@dispersion) | anyNA(.out@coefficients) ) .out@type <- "fitError"

      return(.out)
    }
  }

  mpDat=mpMat[, c(mid, rid, wid)]
  models=apply(mpDat, 1, .fitOne)
  names(models)=mpMat$gene

  for (i in 1:length(models)) models[[i]]@gene=names(models)[i]
  if (noFitError) models=removeFitErrorModels(models)

  return(models)
}

#' Add posterior var to fitted models
#'
#' This function adds posterior var to fitted models by calling limma::squeezeVar.
#'
#' @param models A model list.
#' @param robust Default is FALSE. See limma::squeezeVar's robust.
#' @return A model list with updated var.post and df.post.
#' @export
#' @family QBGLM functions
addPost2FitModels<-function(models, robust=FALSE) {

  # Squeeze a set of sample variances together
  # by computing empirical Bayes posterior means
  hlp <- limma::squeezeVar(
    var = unlist(lapply(models, function(m){m@dispersion} )),
    df = unlist(lapply(models, function(m){m@df.residual} )),  #all df is 2
    robust = robust
  )

 # d=data.frame(cbind(x=unlist(lapply(models, function(m) m@dispersion)), y=hlp$var.post))
 # head(d)
 # write.csv(d, file='var.csv')
 # summary(unlist(lapply(models, function(m) m@dispersion)) - hlp$var.post)

  # put variance and degrees of freedom in appropriate slots
  for (i in seq_along(models)) {
    models[[i]]@var.post <- as.numeric(hlp$var.post[i])

    if (length(hlp$df.prior)>1) #robust=T
      models[[i]]@df.post <- hlp$df.prior[i] + models[[i]]@df.residual
    else
      models[[i]]@df.post <- hlp$df.prior + models[[i]]@df.residual

    if (length(models[[i]]@df.post)==0) models[[i]]@df.post=NA # replace numeric(0) with NA
    models[[i]]@type <- 'glm.post'
  }

  return(models)
}

.getFitModelSummary <- function(m) {
  estimates=m@coefficients
  se=sqrt(diag(m@cov.unscaled * m@var.post))
  #rawse=sqrt(diag(m@cov.unscaled * m@dispersion))
  df=m@df.post
  rawpv=m@raw.pvals
  if (length(se)!=2) stop(m@gene)
  return(c(estimates, se, df, rawpv))
}


#' Test the significance of beta for a fitted model list.
#'
#' This function tests the significance of beta for a fitted model list.
#'
#' @param models A model list after calling fitQuasiGLM and addPost2FitModels.
#' @param qval TRUE to adjust the pvalue by qvalue.
#' @param FDR TRUE to adjust the pvalue by BH adjustment, default is FALSE.
#' @param empirical TRUE to adjust the pvalue by empirical pvalue adjustment of locfdr, default is FALSE.
#' @return A data frame with columns: b0, b1, stdError0, stdError1, df, rawPval0, rawPval1,  tstat0, tstat1, pval0, pval1, and [qval0, qval1, padj0, padj1, emp_pval0, emp_padj0, ..].
#' @export
#' @family QBGLM functions
testBetas <- function(models,
                      qval=TRUE,
                      FDR=FALSE,
                      empirical=FALSE
                      ) {

  mtype = unique(unlist(lapply(models, function(m){m@type} )))
  if (!('glm.post' %in% mtype)) stop("testBetas: please fitQuasiGLM and addPost2FitModels first!")

  models=removeFitErrorModels(models)

  stats=unlist(lapply(models, .getFitModelSummary))
  ncol=length(stats)/length(models)
  if (ncol!=7) stop("error in .getFitModelSummary, not 7 values!\n")
  stats=matrix(stats, ncol=7, byrow=TRUE)

  ## b0, b1
  estimates = stats[, 1:2]
  ## SE of b0, b1
  se=stats[, 3:4]
  ## df.post
  df=stats[, 5]
  ## T statistics of b0, b1
  t <- estimates / se
  ## Pval of b0, b1
  pval <- pt(-abs(t), df) * 2

  res=cbind(stats, t, pval)
  colnames(res)=c('b0','b1','stdError0','stdError1','df','rawPval0','rawPval1', 'tstat0','tstat1','pval0','pval1')
  res=as.data.frame(res)
  res=cbind(gene=unlist(lapply(models, function(m) m@gene)), res)

  if (qval) {
    res$qval0=qvalue::qvalue(res$pval0)$qvalues
    res$qval1=qvalue::qvalue(res$pval1)$qvalues
  }

  if (FDR) {
    res$padj0=stats::p.adjust(res$pval0)
    res$padj1=stats::p.adjust(res$pval1)
  }

  if (empirical) {
    empirical <- p.adjust_empirical(res$pval0,
                                    res$tstat0,
                                    diagplot1 = F,
                                    diagplot2 = F)
    res$emp_pval0 <- empirical$pval
    res$emp_padj0 <- empirical$FDR

    empirical <- p.adjust_empirical(res$pval1,
                                    res$tstat1,
                                    diagplot1 = F,
                                    diagplot2 = F)
    res$emp_pval1 <- empirical$pval
    res$emp_padj1 <- empirical$FDR
  }

  return(res)
}


## ---------- GLM factor (temp) ---------
## Add sample factors to the GLM model
fitQuasiGLM_sample<-function(mpMat, noFitError=TRUE, grp) {

  # m6A, RUD, weights
  mid=grep('_m6A', colnames(mpMat))
  rid=grep('_RUD', colnames(mpMat))
  wid=grep('_weight', colnames(mpMat))

  .fitOne <- function(aRow) {
    if (any(is.na(aRow))) {
      .out <- FitModel(type = "fitError")
      return(.out)
    }

    n=length(aRow)/3
    d=data.frame(m6A=aRow[1:n], rud=aRow[(n+1):(2*n)], weights=aRow[(2*n+1):(3*n)], sample=grp)

    model <- try(glm(rud ~ m6A+sample, family = "quasibinomial", data=d, weights=d$weights))
    #summary(model)

    if (inherits(model, "try-error")) {
      .out <- FitModel(type = "fitError")
      return(.out)
    } else {
      if (any(is.na(model$coefficients))) {
        .out <- FitModel(type = "fitError")
        return(.out)
      }
      type <- "glm"
      .out <- .FitModel(
        type = type,
        coefficients=model$coefficients[1:2],
        raw.pvals=summary(model)$coefficients[1:2, 4], #only retain b0 and b1
        df.residual=model$df.residual,
        dispersion=summary(model)$dispersion,
        cov.unscaled=summary(model)$cov.unscaled, #more than two rows
        var.post = as.numeric(NA),
        df.post = as.numeric(NA)
      )

      if (is.na(.out@dispersion) | anyNA(.out@coefficients) ) .out@type <- "fitError"

      return(.out)
    }
  }

  mpDat=mpMat[, c(mid, rid, wid)]

  models=list()
  for (i in 1:nrow(mpDat)) {
    #cat(mpMat$gene[i],',')
    models[[i]]=.fitOne(unlist(mpDat[i, ]))
  }

  #models=apply(mpDat, 1, .fitOne)
  names(models)=mpMat$gene

  #aRow=unlist(mpDat[mpMat$gene=='105371267', ])

  for (i in 1:length(models)) models[[i]]@gene=names(models)[i]
  if (noFitError) models=removeFitErrorModels(models)

  return(models)
}


## ---------- Output ---------

getFitSummary_sample <- function(models) {

  .getOneSummay<-function(m) {
    estimates=m@coefficients
    rawpv=m@raw.pvals
    return(c(estimates, rawpv))
  }

  stats=unlist(lapply(models, .getOneSummay))
  ncol=length(stats)/length(models)
  if (ncol!=4) stop("error in .getFitModelSummary, not 7 values!\n")
  res=matrix(stats, ncol=4, byrow=TRUE)
  colnames(res)=c('b0','b1','rawPval0','rawPval1')
  res=as.data.frame(res)
  res=cbind(gene=unlist(lapply(models, function(m) m@gene)), res)
  res$rawPval1[is.na(res$rawPval1)]=1
  return(res)
}


#' Get summary of the GLM-fitted results
#'
#' This function count number of genes of statistical significance, with pvalue or adjusted pvalue or qvalue < cutoff.
#'
#' @param resMp A data frame from `testBetas`, with rawPval1, pval1 and (qval1, padj1...).
#' @param pval1 cutoff for pval1, default is 0.05.
#' @param qval1 cutoff for qval1, default is 0.2.
#' @param padj1 cutoff for qadj1, default is 0.1.
#' @return NULL, print statistics on the screen.
#' @export
#' @family output functions
statFitRes<-function(resMp, pval1=0.05, qval1=0.2, padj1=0.1) {
  resMp$rawPval1[is.na(resMp$rawPval1)]=1

  nall=nrow(resMp)

  pvals=resMp$pval1
  np=sum(pvals<pval1)

  npRaw=sum(resMp$rawPval1<pval1)

  npBoth=sum(resMp$rawPval1<pval1 & pvals<pval1)

  qv=c(nall, npRaw, np, npBoth)
  names(qv)=c('N', paste0('rawPval1<',pval1), paste0('pval1<',pval1) , paste0('rawPval1 & pval1<',pval1) )

  if ('qval1' %in% colnames(resMp) & qval1<1) {
    nq1=sum(resMp$qval1<qval1)
    qv[paste0('qval1<',qval1)]=nq1
    cat("min qval1", min(resMp$qval1), '\n')
  }


  if ('padj1' %in% colnames(resMp) & padj1<1) {
    bh1=sum(resMp$padj1<padj1)
    qv[paste0('padj1<', padj1)]=bh1
    cat("min padj1", min(resMp$padj1), '\n')
  }

 # cat(sprintf("Cutoff: %s<%.2f; %s<%.2f; %s<%.2f\n", 'pval1',pthd,'qval1',qthd,'padj1',bthd))

  for (i in 1:length(qv)) cat(qv[i], names(qv)[i],  '\n')
  return(invisible(NULL))
}



#' Get genes' GLM-fitted results
#'
#' getFitDataByGenes gets the m6A/RUD/weight data frame of given gene names  (which can be the gene, gene assembly, gene symbol column).
#'
#' @param fitRes A data frame from `testBetas`, with rawPval1, pval1 and (qval1, padj1...).
#' @param genes genes to be subset.
#' @return A data frame.
#' @examples
#' \dontrun{
#' getFitDataByGenes(res, genes=c(7114, 2107))
#' }
#' @export
#' @family output functions
getFitDataByGenes<-function(fitRes, genes) {
  # m6A, RUD, weights
  mid=grep('_m6A', colnames(fitRes))
  rid=grep('_RUD', colnames(fitRes))
  wid=grep('_weight', colnames(fitRes))

  gid=which(fitRes$gene %in% genes)
  if (length(gid)==0) {
    if ('gene_ensembl' %in% colnames(fitRes))
      gid=which(fitRes$gene_ensembl %in% genes)
  }
  if (length(gid)==0) {
    if ('gene_symbol' %in% colnames(fitRes))
      gid=which(fitRes$gene_symbol %in% genes)
  }
  if (length(gid)==0) {
    return(NULL)
  }

  gids=gid
  dats=data.frame()
  for(gid in gids) {
    gns=unlist(fitRes[gid, grep('gene', colnames(fitRes))])
    dat=data.frame(m6A=unlist(fitRes[gid, mid]), RUD=unlist(fitRes[gid, rid]), weight=unlist(fitRes[gid, wid]))
    cat(gns,'\n')
    print(dat)
    dat=cbind(gene=gns[1], sample=rownames(dat), dat)
    if(nrow(dats)==0) {
      dats=dat
    } else {
      dats=rbind(dats, dat)
    }
  }
  rownames(dats)=NULL
  return(dats)
}


#' Subset GLM-fitted results
#'
#' This function subsets genes of statistical significance, with pvalue or adjusted pvalue or qvalue < cutoff, and simplifies columns.
#'
#' @param resMp A data frame from `testBetas`, with rawPval1, pval1 and (qval1, padj1...).
#' @param rawPval1 cutoff for rawPval1, default is 1.1 (not filtering).
#' @param pval1 cutoff for pval1, default is 1.1 (not filtering).
#' @param qval1 cutoff for qval1, default is 1.1 (not filtering).
#' @param padj1 cutoff for qadj1, default is 1.1 (not filtering).
#' @return A data frame.
#' @export
#' @family output functions
subsetFitRes<-function(resMp, rawPval1=1.1, pval1=1.1, qval1=1.1, padj1=1.1) {

  n1=nrow(resMp)

  resMp$rawPval1[is.na(resMp$rawPval1)]=1
  resMp$pval1[is.na(resMp$pval1)]=1
  resMp=resMp[resMp$rawPval1<rawPval1, ]
  resMp=resMp[resMp$pval1<pval1, ]

  if ('qval1' %in% colnames(resMp)) {
    resMp$qval1[is.na(resMp$qval1)]=1
    resMp=resMp[resMp$qval1<qval1, ]
  }

  if ('padj1' %in% colnames(resMp)) {
    resMp$padj1[is.na(resMp$padj1)]=1
    resMp=resMp[resMp$padj1<padj1, ]
  }

  cols0=grep('0$|d0|nM6A|stdError|df|tstat', colnames(resMp))
  if (length(cols0)>0) resMp=resMp[, -cols0]

  if (('qval1' %in% colnames(resMp))) {
    resMp=resMp[order(resMp$pval1, resMp$rawPval1, decreasing = F), ]
  } else {
    resMp=resMp[order(resMp$qval1, resMp$pval1, resMp$rawPval1, decreasing = F), ]
  }

  if(nrow(resMp)!=n1)
    cat('Sig. number:', n1 , '-->', nrow(resMp),'\n')
  return(resMp)
}


#' Get genomic regions of all genes
#'
#' getGeneAnnos retrieves genomic ranges of all genes from different genome annotation sources.
#'
#' @param annoObj An OrganismDb or Mart object.
#' @return A data frame of genomic ranges of all genes, with fixed columns of chr/strand/start/end and other columns starting with 'gene_' like gene_entrez.
#' @examples
#' \dontrun{
#' library(Homo.sapiens)
#' orgdb=Homo.sapiens
#' geneAnnos=getGeneAnnos(orgdb)
#' }
#' @name getGeneAnnos
#' @export
getGeneAnnos <- function(annoObj) {

  if (inherits(annoObj, "OrganismDb")) {
    genes=suppressMessages(genes(annoObj, columns=c("SYMBOL","ENSEMBL","GENEID")))
    genes=as.data.frame(genes)
    colnames(genes)=c('chr','start','end','width','strand','gene_entrezid','gene_ensembl','gene_symbol')
    for (g in c('gene_entrezid','gene_ensembl','gene_symbol')) {
      genes[, g]=unlist(lapply(genes[, g], '[', 1))
    }
    return(genes[, c('chr','strand','start','end', 'gene_entrezid','gene_ensembl','gene_symbol')])
  }

  if (inherits(annoObj, "Mart")) {
    #listAttributes(bm)[grep('gene_id', listAttributes(bm)$name), ]
    genes=biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id", "chromosome_name",'strand','start_position','end_position'), mart = annoObj)
    colnames(genes)=c('gene_ensembl','gene_symbol','gene_entrezid', 'chr','strand','start','end')
    genes$strand[genes$strand==1]='+'
    genes$strand[genes$strand==-1]='-'
    return(genes[, c('chr','strand','start','end', 'gene_entrezid','gene_ensembl','gene_symbol')])
  }

  cat(
  "The function for annotation type",
  paste(class(annoObj), collapse = ", "),
  "has not yet been realized\n")
  return(data.frame())
}


## ---------- Wrapper funs ---------

#' Wrapper function to convert APA's PACdataset to PACdataset for QBGLM.
#'
#' APAds2GLMds is a wrapper function to convert APA's PACdataset to PACdataset for QBGLM. It annotates polyA sites with TXDB and get proximal and distal sites of 3'UTR-APA genes. Finally, it filters DE-APA by fisher.test of replicates.
#'
#' @param APAds A pacdataset of APA sites.
#' @param PDminDist Default is 50. Require the distance between  proximal and distal polyA site to be > minDist and < maxDist.
#' @param PDmaxDist Default is 1000. Require the distance between  proximal and distal polyA site to be > minDist and < maxDist.
#' @param PDminRatio Default is 0.05. If >0 then will filter out polyA sites that do not meet the requirements based on the ratio first, and then filter proximal and distal sites.
#' @param PDfixDistal Default is FALSE. True to fix the farthest PA, that is, only look for the proximal site.
#' @param DEpval cutoff to filter DE-APA, default is 0.05. If =1, then not filter DE.
#' @param txdb An TXDB object for genome annotation.
#' @param extUTRlen Length to extend 3'UTR to include more near downstream polyA sites.
#' @return  A PACdataset with only proximal and distal, and DE-APA (if DEpval<1) sites.
#' @export
#' @examples
#' \dontrun{
#' ## first, get a APAds from QAPA results
#' # add raw counts from .sf files to QAPA list
#' qfile='M3KD.qapa.txt'
#' sfs=c('qapa_SRR847370_quant.sf', 'qapa_SRR847371_quant.sf',
#'       'qapa_SRR847374_quant.sf', 'qapa_SRR847375_quant.sf')
#' names(sfs)=paste0('sample', 1:4)
#' q=addQAPARawCounts(qfile, sfs, suffix='CNT', ofile=NULL)
#' head(q)
#'
#' # convert QAPA list to PACdataset
#' apads=QAPA2PACds(q, vcol='CNT')
#'
#' # change sample names: ctrlM3 and treatM3 and 1/2 for reps
#' apads=setPACdsSmpInfo(apads,
#'                 smpInfo=cbind(old=paste0('sample', 1:4),
#'                               new=c('ctrlM31','ctrlM32','treatM31','treatM32'),
#'                               group=c('ctrlM3','ctrlM3','treatM3','treatM3')))
#' apads@anno$gene_QAPA=apads@anno$gene
#' movAPA::summary(apads)
#'
#' ## get GLMds
#' apads=APAds2GLMds(APAds=apads,
#'           PDminDist=50, PDmaxDist=5000, PDminRatio=0.05, PDfixDistal=FALSE,
#'           DEpval=0.05,
#'           txdb=txdb,
#'           extUTRlen=2000)
#' }
#' @family APAds functions
APAds2GLMds<-function(APAds,
                     PDminDist=50, PDmaxDist=5000, PDminRatio=0.05, PDfixDistal=FALSE,
                     DEpval=0.05,
                     txdb,
                     extUTRlen=2000) {

  ## annotate apads and m6ads
  apads=movAPA::annotatePAC(APAds, txdb)

  ## Extend 3'UTR to include more APA and m6A sites in downstream 3'UTR
  apads=movAPA::ext3UTRPACds(apads, extUTRlen)

  # Select two "best" pAs in 3'UTR as proximal and distal pA
  apadsUTR=movAPA::get3UTRAPAds(apads)
  apadsPD=get3UTRAPApd(apadsUTR, minDist=PDminDist, maxDist=PDmaxDist, minRatio=PDminRatio, fixDistal=PDfixDistal, addCols='pd')
  cat(sprintf("PD genes# %d\n",length(apadsPD) ) )

  if (DEpval<1) {
    # calculate DEAPA by fisher.test for each replicate
    apadsUTRDE=get3UTRAPApdDE(apadsPD, pthd=DEpval, filterDE=TRUE)
    cat(sprintf("DE-APA PD sites# %d\n",length(apadsUTRDE) ) )
    return(apadsUTRDE)
  } else {
    return(apadsPD)
  }
}



#' Wrapper function to load QAPA results to PACdataset for QBGLM.
#'
#' QAPA2GLMds is a wrapper function to load QAPA results to PACdataset for QBGLM. It first get raw counts of polyA sites by reading qapafile and quant.sf files.
#' Then it annotates polyA sites with TXDB and get proximal and distal sites of 3'UTR-APA genes. Finally, it filters DE-APA by fisher.test of replicates.
#'
#' @param qapaFile A file name of the QAPA APA list.
#' @param sfFiles A vector of quant.sf files of QAPA's outputs which record raw counts. Names of `sfFiles` should be sample names in the `qapaFile`.
#' @param newSampleNames new sample names for the sfFiles, e.g., c('ctrlM31','ctrlM32','treatM31','treatM32').
#' @param newSampleGroups new group names for the sfFiles, e.g., c('ctrlM3','ctrlM3','treatM3','treatM3').
#' @param PDminDist Default is 50. Require the distance between  proximal and distal polyA site to be > minDist and < maxDist.
#' @param PDmaxDist Default is 1000. Require the distance between  proximal and distal polyA site to be > minDist and < maxDist.
#' @param PDminRatio Default is 0.05. If >0 then will filter out polyA sites that do not meet the requirements based on the ratio first, and then filter proximal and distal sites.
#' @param PDfixDistal Default is FALSE. True to fix the farthest PA, that is, only look for the proximal site.
#' @param DEpval cutoff to filter DE-APA, default is 0.05. If =1, then not filter DE.
#' @param txdb An TXDB object for genome annotation.
#' @param extUTRlen Length to extend 3'UTR to include more near downstream polyA sites.
#' @return  A PACdataset with only proximal and distal, and DE-APA (if DEpval<1) sites.
#' @export
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb=TxDb.Hsapiens.UCSC.hg38.knownGene
#' qapaFile='M3KD.qapa.txt'
#' sfFiles=c('qapa_SRR847370_quant.sf', 'qapa_SRR847371_quant.sf',
#'          'qapa_SRR847374_quant.sf', 'qapa_SRR847375_quant.sf')
#' names(sfFiles)=paste0('sample', 1:4)
#' newSampleNames=c('ctrlM31','ctrlM32','treatM31','treatM32')
#' newSampleGroups=c('ctrlM3','ctrlM3','treatM3','treatM3')
#' apads=QAPA2GLMds(qapaFile, sfFiles,
#'           newSampleNames, newSampleGroups,
#'           PDminDist=50, PDmaxDist=5000, PDminRatio=0.05, PDfixDistal=FALSE,
#'           DEpval=0.05,
#'           txdb=txdb,
#'           extUTRlen=2000)
#' }
#' @family APAds functions
QAPA2GLMds<-function(qapaFile, sfFiles,
                     newSampleNames, newSampleGroups,
                     PDminDist=50, PDmaxDist=5000, PDminRatio=0.05, PDfixDistal=FALSE,
                     DEpval=0.05,
                     txdb,
                     extUTRlen=2000) {

  if (is.null(names(sfFiles))) stop("Please set names of sfFiles as the sample names in qapaFile!")

  # add raw counts from .sf files to QAPA list
  q=addQAPARawCounts(qapaFile, sfFiles, suffix='CNT', ofile=NULL)

  # convert QAPA list to PACdataset
  apads=QAPA2PACds(q, vcol='CNT')

  # change sample names: ctrlM3 and treatM3 and 1/2 for reps
  apads=setPACdsSmpInfo(apads,
                        smpInfo=cbind(old=names(sfFiles),
                                      new=newSampleNames,
                                      group=newSampleGroups))
  apads@anno$gene_QAPA=apads@anno$gene
  #movAPA::summary(apads)

  apads=APAds2GLMds(APAds=apads,
                    PDminDist=PDminDist, PDmaxDist=PDmaxDist, PDminRatio=PDminRatio, PDfixDistal=PDfixDistal,
                    DEpval=DEpval,
                    txdb=txdb,
                    extUTRlen=extUTRlen)
  return(apads)
}

#' Wrapper function to load m6A's PACdataset to PACdataset for QBGLM.
#'
#' m6Ads2GLMds is a wrapper function to load m6A's PACdataset to PACdataset for QBGLM. It annotates m6A sites with TXDB. Finally, it filters DE-m6As.
#'
#' @param m6Ads A pacdataset of m6A peaks.
#' @param DEpval cutoff to filter DE-m6A, default is 0.05. If =1, then not filter DE.
#' @param txdb An TXDB object for genome annotation.
#' @param extUTRlen Length to extend 3'UTR to include more near downstream m6A sites.
#' @return A PACdataset of m6A peaks or DE-m6A peaks (if DEpval<1), with @anno have raw_<IP, Input count columns>, DE_log.OR/pvalue/padj.
#' @export
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb=TxDb.Hsapiens.UCSC.hg38.knownGene
# A table of Get_peak_infor of m6Aexpress.
#' m6Afile='M3KD.m6A.txt'
#'
#' # libsize of m6Aexpress, with two columns: sample and libsize.
#' libSizes=read.table('M3KD.m6A.libsize.txt', header=TRUE)
#'
#' # new sample and group names for samples in libSizes
#' newSampleNames=c('ctrlM31','ctrlM32','treatM31','treatM32')
#' newSampleGroups=c('ctrlM3','ctrlM3','treatM3','treatM3')
#'
#' # get m6A for QBGLM
#' m6ads=m6A2GLMds(m6Afile, libSizes,
#'                 newSampleNames, newSampleGroups,
#'                 DEpval=0.05,
#'                 txdb=txdb,
#'                 extUTRlen=2000)
#'
#' movAPA::summary(m6ads)
#' m6Ads2GLMds(m6Ads=m6ads, DEpval=0.05, txdb, extUTRlen=2000)
#' }
#' @family m6Adata functions
m6Ads2GLMds<-function(m6Ads,
                      DEpval=0.05,
                      txdb,
                      extUTRlen=2000) {

  # After doDE, we can filter DE-m6A genes by DE_pvalue
  if (DEpval<1) {
    DEMgenes=m6Ads@anno$gene[m6Ads@anno$DE_pvalue<DEpval]
    m6ads_DEM=movAPA::subsetPACds(m6Ads, genes=DEMgenes)
    cat(sprintf("m6A:%d\nDE-m6A:%d\n",length(m6Ads), length(m6ads_DEM) ) )

    ## annotate m6ads_DEM
    m6ads_DEM=movAPA::annotatePAC(m6ads_DEM, txdb)

    ## Extend 3'UTR to include more APA and m6A sites in downstream 3'UTR
    m6ads_DEM=movAPA::ext3UTRPACds(m6ads_DEM, extUTRlen)

    if (sum(is.na(m6ads_DEM@anno$gene))>0)
      m6ads_DEM=m6ads_DEM[!is.na(m6ads_DEM@anno$gene)]

    return(m6ads_DEM)
  } else {
    # annotate m6A with txdb
    m6Ads=movAPA::annotatePAC(m6Ads, txdb)

    ## Extend 3'UTR to include more APA and m6A sites in downstream 3'UTR
    m6Ads=movAPA::ext3UTRPACds(m6Ads, extUTRlen)

    ## for some reason, some m6A cannot be annotated to any genes, which will be deleted
    if (sum(is.na(m6Ads@anno$gene))>0)
      m6Ads=m6Ads[!is.na(m6Ads@anno$gene)]

    return(m6Ads)
  }

}



#' Wrapper function to load m6Aexpress' results to PACdataset for QBGLM.
#'
#' m6Apeak2GLMds is a wrapper function to load m6Aexpress' results to PACdataset for QBGLM. It first get m6A scores by reading m6Afile and libSizes.
#' Then it annotates m6A sites with TXDB. Finally, it filters DE-m6As.
#'
#' @param m6Afile A m6AExpress's m6A peak file, with columns seqnames, start, end, width, strand, gene_name, and <IP1, IP2, Treated_IP1, Treated_IP2..>.
#' @param libSizes A file name or a data frame with 2 columns: sample and counts. The first column is the sample column, with the order -- ctrls: IP, Input; treats: IP, Input.
#' The second columns is the total read counts for each sample.
#' @param newSampleNames new sample names, e.g., c('ctrlM31','ctrlM32','treatM31','treatM32'), should be as the order of libSizes.
#' @param newSampleGroups new group name, e.g., c('ctrlM3','ctrlM3','treatM3','treatM3'), should be as the order of libSizes.
#' @param DEpval cutoff to filter DE-m6A, default is 0.05. If =1, then not filter DE.
#' @param txdb An TXDB object for genome annotation.
#' @param extUTRlen Length to extend 3'UTR to include more near downstream m6A sites.
#' @return A PACdataset of m6A peaks or DE-m6A peaks (if DEpval<1), with @anno have raw_<IP, Input count columns>, DE_log.OR/pvalue/padj.
#' @export
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb=TxDb.Hsapiens.UCSC.hg38.knownGene
#' m6Afile='M3KD.m6A.txt'
#' libSizes=read.table('M3KD.m6A.libsize.txt', header=TRUE)
#' newSampleNames=c('ctrlM31','ctrlM32','treatM31','treatM32')
#' newSampleGroups=c('ctrlM3','ctrlM3','treatM3','treatM3')
#' m6A2GLMds(m6Afile, libSizes,
#'           newSampleNames, newSampleGroups,
#'           DEpval=0.05,
#'           txdb=txdb,
#'           extUTRlen=2000)
#' }
#' @family m6Adata functions
m6Apeak2GLMds<-function(m6Afile, libSizes,
                    newSampleNames, newSampleGroups,
                    DEpval=0.05,
                    txdb,
                    extUTRlen=2000) {

  peakDf=read.table(m6Afile, header = TRUE)

  if (is.character(libSizes)) libSizes=read.table(libSizes, header = TRUE)

  ctrls=libSizes[,1][1:(nrow(libSizes)/2)] # "IP1"    "IP2"    "Input1" "Input2"
  treats=libSizes[,1][((nrow(libSizes)/2)+1) : nrow(libSizes)] # "Treated_IP1"    "Treated_IP2"    "Treated_Input1" "Treated_Input2"

  if (!all(c(ctrls, treats) %in% names(peakDf))) stop("Sample names in libSizes[1, ] should be all in column names of m6Afile!")

  # convert to PACdataset and detect DE-m6A, but without filtering
  m6ads=m6AExpress2PACds(peakDf=peakDf, ctrls=ctrls, treats=treats, libSizes=libSizes[, 2], doDE=TRUE, filterDE=FALSE)

  # change sample names: ctrlM3 and treatM3 and 1/2 for reps
  old=c(ctrls[1:(length(ctrls)/2)], treats[1:(length(treats)/2)])       #"IP1", "IP2", "Treated_IP1", "Treated_IP2"
  m6ads=setPACdsSmpInfo(m6ads,
                        smpInfo=cbind(old=old,
                                      new=newSampleNames,
                                      group=newSampleGroups))

  m6ads=m6Ads2GLMds(m6Ads=m6ads,
                    DEpval=DEpval,
                    txdb=txdb,
                    extUTRlen=extUTRlen)
  return(m6ads)

}


#' Wrapper function of QBGLM
#'
#' QBGLM is a wrapper function of QBGLM using APA from `QAPA2GLMds` and m6A `m6A2GLMds`. It first links APAds and m6Ads to get per gene data.
#' And then do QBGLM to test each gene and get betas.
#'
#' @param APAds A PACdataset from QAPA2GLMds, txdb-annotated DE-PD-APA dataset.
#' @param m6Ads A PACdataset from m6A2GLMds, txdb-annotated DE-m6A dataset.
#' @return A data frame with each row a gene.
#' @export
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb=TxDb.Hsapiens.UCSC.hg38.knownGene
#'
#' # output APA list from QAPA, with TPM for each sample
#' qapaFile='M3KD.qapa.txt'
#'
#' # quant.sf: raw counts of salmon from QAPA's output
#' sfFiles=c('qapa_SRR847370_quant.sf', 'qapa_SRR847371_quant.sf',
#'           'qapa_SRR847374_quant.sf', 'qapa_SRR847375_quant.sf')
#' names(sfFiles)=paste0('sample', 1:4)
#'
#' # new sample and group names for sfFiles
#' newSampleNames=c('ctrlM31','ctrlM32','treatM31','treatM32')
#' newSampleGroups=c('ctrlM3','ctrlM3','treatM3','treatM3')
#'
#' # get APA for QBGLM
#' apads=QAPA2GLMds(qapaFile=qapaFile, sfFiles=sfFiles,
#'                  newSampleNames=newSampleNames, newSampleGroups=newSampleGroups,
#'                  PDminDist=50, PDmaxDist=5000, PDminRatio=0.05, PDfixDistal=FALSE,
#'                  DEpval=0.05,
#'                  txdb=txdb,
#'                  extUTRlen=2000)
#' movAPA::summary(apads)
#'
#'
#' # A table of Get_peak_infor of m6Aexpress.
#' m6Afile='M3KD.m6A.txt'
#'
#' # libsize of m6Aexpress, with two columns: sample and libsize.
#' libSizes=read.table('M3KD.m6A.libsize.txt', header=TRUE)
#'
#' # new sample and group names for samples in libSizes
#' newSampleNames=c('ctrlM31','ctrlM32','treatM31','treatM32')
#' newSampleGroups=c('ctrlM3','ctrlM3','treatM3','treatM3')
#'
#' # get m6A for QBGLM
#' m6ads=m6A2GLMds(m6Afile, libSizes,
#'                 newSampleNames, newSampleGroups,
#'                 DEpval=0.05,
#'                 txdb=txdb,
#'                 extUTRlen=2000)
#'
#' movAPA::summary(m6ads)
#'
#' # QBGLM
#' resMp=QBGLM(APAds=apads, m6Ads=m6ads)
#' nrow(resMp)
#' # just simplify columns without filtering
#' resMp=subsetFitRes(resMp, qval1=2)
#'
#' ## Add gene symbols and entrezids
#' library(Homo.sapiens)
#' orgdb=Homo.sapiens
#' geneAnnos=getAnnoGenes(orgdb)
#' resMp=merge(resMp,
#'            geneAnnos[, c('gene_entrezid','gene_ensembl','gene_symbol')],
#'            all.x=TRUE, all.y=FALSE, by.x='gene', by.y='gene_entrezid')
#' head(resMp)
#' write.csv(resMp, file='M3KD_HeLa1.resMp.csv')
#'
#' ## subset genes with qval<0.2
#' resMpSig=subsetFitRes(resMp, qval1=0.2)
#' write.csv(resMpSig, file='M3KD_HeLa1.m6APAreg.csv')
#' }
#' @family QBGLM functions
QBGLM<-function(APAds, m6Ads) {

  ## Using DE-APA and DE-m6A to calculate m6A scores and RUD per gene

  # Calculate RUD (distal PA usage) for each gene
  p=getRUDperGene(APAds)

  # Obtain the m6A level of each gene in each sample based on decoy function by distance from the nearest APA site.
  # Then combine all m6A levels under one gene to form a single m6A score for each gene.
  m=getM6AperGene(m6Ads, APAds)

  # get a data frame, with each row denoting one gene recording m6A, RUD, weights per gene
  mpAll=merge(m, p, by.x='gene', by.y='gene')

  ## Fitting a GLM for each gene
  modelsFit=fitQuasiGLM(mpAll)
  models=addPost2FitModels(modelsFit, robust=F)
  res=testBetas(models, qval=TRUE, FDR=FALSE, empirical = FALSE)

  ## add gene+APA+m6A info to the pvalue list
  resMp=merge(res, mpAll, by.x='gene', by.y='gene')
  # just simplify columns without filtering
  resMp=subsetFitRes(resMp, qval1=2)

  ## count number of sig. genes
  statFitRes(resMp, pval1=0.05, qval1=0.2, padj1=1)

  return(resMp)
}


#' Wrapper function of QBGLM for multiple samples
#'
#' QBGLMmultiDs is a wrapper function of QBGLM for multiple samples, using APA from `QAPA2GLMds` and m6A `m6A2GLMds`.
#' It first combines multiple APAds and multiple m6Ads. And then links APAds and m6Ads to get per gene data.
#' And then do QBGLM to test each gene and get betas.
#'
#' @param APAdsList A PACdataset list of APAds from functions like `QAPA2PACds`, without annotation and DE calculation.
#' @param m6AdsList A PACdataset list of m6Ads from `m6AExpress2PACds`, which have done DE.
#' @param distAPA distance to group nearby polyA sites in APAds, default is 24 nt.
#' @param distM6A distance to group nearby m6A sites in m6Ads, default is 50 nt.
#' @param txdb An TXDB object for genome annotation.
#' @param extUTRlen Length to extend 3'UTR to include more near downstream polyA sites.
#' @param DEm6Apval cutoff to filter DE-m6A, default is 0.05.
#' @param PDminDist Default is 50. Require the distance between  proximal and distal polyA site to be > minDist and < maxDist.
#' @param PDmaxDist Default is 1000. Require the distance between  proximal and distal polyA site to be > minDist and < maxDist.
#' @param PDminRatio Default is 0.05. If >0 then will filter out polyA sites that do not meet the requirements based on the ratio first, and then filter proximal and distal sites.
#' @param PDfixDistal Default is FALSE. True to fix the farthest PA, that is, only look for the proximal site.
#' @param DEAPApval cutoff to filter DE-APA, default is 0.05.
#' @return A data frame with each row a gene.
#' @export
#' @family QBGLM functions
QBGLMmultiDs<-function(APAdsList, m6AdsList,
                       distAPA=24, distM6A=50,
                       txdb, extUTRlen=2000,
                       DEm6Apval=0.05,
                       PDminDist=50, PDmaxDist=5000,
                       PDminRatio=0.05, PDfixDistal=FALSE,
                       DEAPApval=0.05 ) {

  ### APA ###
  ## merge APAds list to one APAds
  apadsm=mergePACdsCoords(APAdsList, d=distAPA)

  ## annotate APAds
  apadsm=movAPA::annotatePAC(apadsm, txdb)
  apadsm=movAPA::ext3UTRPACds(apadsm, extUTRlen)

  ## get proximal and distal
  apadsmUTR=movAPA::get3UTRAPAds(apadsm)
  apadsmPD=get3UTRAPApd(apadsmUTR, minDist=PDminDist, maxDist=PDmaxDist, minRatio=PDminRatio, fixDistal=PDfixDistal, addCols='pd')

  # get DE-APA (slow)
  cat("Get DE-APA for merged APAds by pair-replicate fisher.test, this may take several minutes...\n")
  apadsmDE=get3UTRAPApdDE(apadsmPD, pthd=DEAPApval, filterDE=TRUE)

  # get RUD
  RUD=getRUDperGene(apadsmDE)

  ### m6A ###

  ## get DEm6A-gene names for each m6Ads
  DEm6Agenes=c()
  for (i in 1:length(m6AdsList)) {
    g=unique(m6AdsList[[i]]@anno$gene[m6AdsList[[i]]@anno$DE_pvalue<DEm6Apval])
    DEm6Agenes=c(DEm6Agenes, g)
    cat(sprintf("%s DE-m6A genes# %d\n", names(m6AdsList)[i], length(g)))
  }
  DEm6Agenes=unique(DEm6Agenes)
  cat(sprintf("Total DE-m6A genes# %d\n", length(DEm6Agenes)))

  for (i in 1:length(m6AdsList)) {
    m6AdsList[[i]]=movAPA::subsetPACds(m6AdsList[[i]], genes=DEm6Agenes)
  }

  ## merge m6Ads list to one m6Ads
  m6adsmDE=mergePACdsRanges(PACdsList=m6AdsList, d=distM6A)
  #table(m6adsm@anno$nPA)

  ## subset merged m6A genes to get only DE-m6As
  #m6adsmDE=movAPA::subsetPACds(m6adsm, genes=DEm6Agenes)

  ## annotate m6A to the same gene system as APA
  ## it is not important which region m6A located
  m6adsmDE=movAPA::annotatePAC(m6adsmDE, txdb)
  if (sum(is.na(m6adsmDE@anno$gene))>0)
    m6adsmDE=m6adsmDE[!is.na(m6adsmDE@anno$gene)]
  m6adsmDE=movAPA::ext3UTRPACds(m6adsmDE, extUTRlen)

  # get m6A per gene by decoy function
  m6A=getM6AperGene(m6adsmDE, apadsmDE)

  # get m6A, RUD weights per gene
  mpAll=merge(m6A, RUD, by.x='gene', by.y='gene')

  ## Fitting a GLM for each gene
  modelsFit=fitQuasiGLM(mpAll)
  models=addPost2FitModels(modelsFit, robust=F)
  res=testBetas(models, qval=TRUE, FDR=FALSE, empirical = FALSE)

  ## add gene+APA+m6A info to the pvalue list
  resMp=merge(res, mpAll, by.x='gene', by.y='gene')
  # just simplify columns without filtering
  resMp=subsetFitRes(resMp, qval1=2)

  ## count number of sig. genes
  statFitRes(resMp, pval1=0.05, qval1=0.2, padj1=1)

  return(resMp)

}


## ---------- temp ---------

#' Summarize p-value and q-value results
#'
#' @param resMp A data.frame containing p-value columns.
#' @param pthd The p-value threshold. Default is 0.05.
#' @param qthd The q-value threshold. Default is 0.2.
#'
#' @return This function prints summary statistics and returns invisibly.
#'
#' @export
.statPQ <- function(resMp, pthd=0.05, qthd=0.2) {

  resMp$rawPval1[is.na(resMp$rawPval1)]=1
  
  nall=nrow(resMp)
  
  pvals=resMp$pval1
  np=sum(pvals<pthd)  
  
  npRaw=sum(resMp$rawPval1<pthd)  
  
  npBoth=sum(resMp$rawPval1<pthd & pvals<pthd) 
  
  ## Use the minimum p-value to calculate q-values
  #pvals=ifelse(resMp$rawPval1>resMp$pval1, resMp$pval1,resMp$rawPval1)
  #hist(pvals)
  
  #pvals=ifelse(resMp$pval1<0.05, resMp$pval1*0.01, resMp$pval1)
  #par(mfrow=c(1,2))
  #hist(resMp$pval1, breaks=40)
  #hist(pvals, breaks=40)
  
  ## Method 1: Calculate q-values without using covariates
  q1=qvalue::qvalue(pvals)
  nq1=sum(q1$qvalues<qthd) 
  
  ## Method 2: Benjamini-Hochberg FDR adjustment;
  ## this method yields fewer significant results
  bh1=p.adjust(pvals, method = "BH")
  nbh=sum(bh1<qthd) #6
  
  # ## Calculate q-values using covariates
  # ### Use nM6A as the only covariate
  # q2 <- swfdr::lm_qvalue(pvals, X=resMp$nM6A)
  # nq2=sum(q2$qvalues<qthd) 
  # 
  # ### Use both nM6A and weights as covariates
  # weights=rowSums(resMp[, grep('weight', colnames(resMp))])
  # q3 <- swfdr::lm_qvalue(pvals, X=cbind(resMp$nM6A, weights))
  # nq3=sum(q3$qvalues<qthd)
  
  #pval_N_q3=sum(pvals<qthd & q3$qvalues<qthd)
  
  # qv=c('N','rawPv', 'pval','pval_n_rawPv', 'bhFDR','qvalue1','qvalue2','qvalue3', 'pval_N_q3')
  # N=c(nall, npRaw, np, npBoth, nbh, nq1, nq2, nq3, pval_N_q3)
  
  qv=c('N','rawPv', 'pval','pval_n_rawPv', 'bhFDR','qvalue1')
  N=c(nall, npRaw, np, npBoth, nbh, nq1)
  
  cat('pthd=', pthd,'; qthd=',qthd, '\nmin qvalue1=', min(q1$qvalues), '\n')
  for (i in 1:length(qv)) cat(qv[i], N[i], '\n')
}


#' Extract gene-level m6A, RUD, and weight data
#'
#' Given a gene identifier, this function searches the gene, gene_ensembl,
#' and gene_symbol columns and extracts gene-level m6A, RUD, and weight data.
#'
#' @param resMp A data.frame containing gene identifiers and result columns.
#' @param gene A gene identifier.
#'
#' @return A data.frame containing m6A, RUD, and weight values.
#' Returns NULL if the gene cannot be found.
#'
#' @export
.getGeneData <- function(resMp, gene) {
  # m6A, RUD, and weights
  mid=grep('_m6A', colnames(resMp))
  rid=grep('_RUD', colnames(resMp))
  wid=grep('_weight', colnames(resMp))
  
  gid=which(resMp$gene==gene)
  if (length(gid)==0) {
    if ('gene_ensembl' %in% colnames(resMp))
      gid=which(resMp$gene_ensembl==gene)
  }
  if (length(gid)==0) {
    if ('gene_symbol' %in% colnames(resMp))
      gid=which(resMp$gene_symbol==gene)
  }
  if (length(gid)==0) {
    return(NULL)
  }
  
  dat=data.frame(
    m6A=unlist(resMp[gid, mid]),
    RUD=unlist(resMp[gid, rid]),
    weight=unlist(resMp[gid, wid])
  )
  cat(unlist(resMp[gid, grep('gene', colnames(resMp))]),'\n')
  print(dat)
}


#' Compare significant genes between two result objects
#'
#' Compare the overlapping genes between two result data.frames.
#' If both pval1 and rawPval1 are available, their intersection is used.
#' Otherwise, rawPval1 is used.
#'
#' @param res1 The first result data.frame.
#' @param res2 The second result data.frame.
#'
#' @return A character vector containing overlapping gene identifiers.
#'
#' @export
.compareTwoResMp <- function(res1, res2) {
  has_pval1 <- all(
    c("pval1", "rawPval1") %in% colnames(res1),
    c("pval1", "rawPval1") %in% colnames(res2)
  )

  if (has_pval1) {
    res1$rawPval1[is.na(res1$rawPval1)] <- 1
    res2$rawPval1[is.na(res2$rawPval1)] <- 1

    res1$pval1[is.na(res1$pval1)] <- 1
    res2$pval1[is.na(res2$pval1)] <- 1

    g1 <- res1$gene[
      res1$rawPval1 < 0.05 &
        res1$pval1 < 0.05
    ]

    g2 <- res2$gene[
      res2$rawPval1 < 0.05 &
        res2$pval1 < 0.05
    ]

    cat("rawPval1 < 0.05 & pval1 < 0.05\n")
  } else {
    if (!all(c("rawPval1", "gene") %in% colnames(res1)) ||
        !all(c("rawPval1", "gene") %in% colnames(res2))) {
      stop(
        "Both result objects must contain gene and rawPval1 columns."
      )
    }

    res1$rawPval1[is.na(res1$rawPval1)] <- 1
    res2$rawPval1[is.na(res2$rawPval1)] <- 1

    g1 <- res1$gene[res1$rawPval1 < 0.05]
    g2 <- res2$gene[res2$rawPval1 < 0.05]

    cat("rawPval1 < 0.05\n")
  }

  overlap <- intersect(g1, g2)

  cat(
    "gene#: ",
    length(g1),
    " ",
    length(g2),
    "\n",
    sep = ""
  )

  cat(
    "ovp#: ",
    length(overlap),
    "\n",
    sep = ""
  )

  invisible(overlap)
}

#' Subset results using p-value and q-value thresholds
#'
#' Filter a result data.frame using raw p-value, adjusted p-value,
#' and q-value thresholds.
#'
#' @param resMp A data.frame containing result statistics.
#' @param rawPval1 Raw p-value threshold. Default is 1.
#' @param pval1 Adjusted or moderated p-value threshold. Default is 1.
#' @param qval1 Q-value threshold. Default is 1.
#'
#' @return A filtered and sorted data.frame.
#'
#' @export
.subsetResMp<-function(resMp, rawPval1=1, pval1=1, qval1=1) {
  resMp$rawPval1[is.na(resMp$rawPval1)]=1
  resMp$pval1[is.na(resMp$pval1)]=1
  
  if (!('qval1' %in% colnames(resMp))) resMp$qval1=qvalue::qvalue(resMp$pval1)$qvalues
  resMp$qval1[is.na(resMp$qval1)]=1
  
  n1=nrow(resMp)
  
  if (qval1<1) {
    resMp=resMp[resMp$qval1<qval1, ]
  }
  
  if (rawPval1<1) {
    resMp=resMp[resMp$rawPval1<rawPval1, ]
  }
  
  if (pval1<1) {
    resMp=resMp[resMp$pval1<pval1, ]
  }
  
  cols0=grep('0$|d0|nM6A|stdError|df|tstat', colnames(resMp))
  if (length(cols0)>0) resMp=resMp[, -cols0]
  resMp=resMp[order(resMp$pval1, resMp$rawPval1, decreasing = F), ]
  cat('After filtering pval/qval:', n1 , '-->', nrow(resMp),'\n')
  return(resMp)
}

#' Convert an m6A expression data.frame to a PACdataset
#'
#' @param peakDf A data.frame containing normalized m6A expression values.
#' @param smpCols Character vector of sample columns.
#'
#' @return A PACdataset object.
#'
#' @keywords internal
.m6AExpressDf2PACds <- function(peakDf, smpCols) {
  if (!is.data.frame(peakDf)) {
    stop("peakDf must be a data.frame.")
  }

  if (!all(smpCols %in% colnames(peakDf))) {
    missing_cols <- setdiff(smpCols, colnames(peakDf))

    stop(
      "The following sample columns are missing from peakDf: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  m6AExpressDf2PACds(
    peakDf,
    smpCols = smpCols
  )
}


#' Normalize m6A IP/Input signals
#'
#' Normalize IP and Input count columns by library size and calculate
#' log2-normalized IP/Input values.
#'
#' @param peakDf A data.frame containing IP and Input count columns.
#' @param ctrls Character vector of control sample names.
#' @param treats Optional character vector of treated sample names.
#' @param libSizes Optional numeric vector of library sizes.
#' @param rawPrefix Prefix added to the original count columns.
#'
#' @return A data.frame containing normalized m6A values.
#'
#' @keywords internal
.m6AIpInputNormalize <- function(
    peakDf,
    ctrls,
    treats = NULL,
    libSizes = NULL,
    rawPrefix = "raw_") {
  if (!is.data.frame(peakDf)) {
    stop("peakDf must be a data.frame.")
  }

  if (!is.character(ctrls) || length(ctrls) == 0L) {
    stop(
      "ctrls must be a non-empty character vector."
    )
  }

  if (is.null(treats)) {
    treats <- character(0)
  }

  if (!is.character(treats)) {
    stop(
      "treats must be NULL or a character vector."
    )
  }

  if (length(treats) > 0L &&
      length(ctrls) != length(treats)) {
    stop(
      "ctrls and treats must have the same length."
    )
  }

  sample_names <- c(ctrls, treats)

  if (!all(sample_names %in% colnames(peakDf))) {
    missing_cols <- setdiff(
      sample_names,
      colnames(peakDf)
    )

    stop(
      "The following sample columns are missing from peakDf: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  .checkIPInputNames(ctrls)

  if (length(treats) > 0L) {
    .checkIPInputNames(treats)
  }

  if (is.null(libSizes)) {
    libSizes <- colSums(
      peakDf[, sample_names, drop = FALSE],
      na.rm = TRUE
    )
  }

  libSizes <- as.numeric(libSizes)

  if (length(libSizes) != length(sample_names)) {
    stop(
      "libSizes must have the same length as the available sample columns."
    )
  }

  if (any(!is.finite(libSizes)) ||
      any(libSizes <= 0)) {
    stop(
      "libSizes must contain positive finite values."
    )
  }

  sf <- libSizes / exp(mean(log(libSizes)))

  peakDf[, sample_names] <- as.data.frame(
    t(
      t(
        peakDf[, sample_names, drop = FALSE]
      ) / sf
    )
  )

  ip_cols <- .getIPnames(ctrls)
  input_cols <- .getInputNames(ctrls)

  if (length(treats) > 0L) {
    ip_cols <- c(
      ip_cols,
      .getIPnames(treats)
    )

    input_cols <- c(
      input_cols,
      .getInputNames(treats)
    )
  }

  required_cols <- c(ip_cols, input_cols)

  if (!all(required_cols %in% colnames(peakDf))) {
    missing_cols <- setdiff(
      required_cols,
      colnames(peakDf)
    )

    stop(
      "The following IP/Input columns are missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  ip_input <- log2(
    (
      peakDf[, ip_cols, drop = FALSE] + 0.01
    ) / (
      peakDf[, input_cols, drop = FALSE] + 0.01
    )
  )

  ip_input[ip_input < 0] <- 0.001
  colnames(ip_input) <- ip_cols

  if (!is.null(rawPrefix)) {
    raw_cols <- paste0(
      rawPrefix,
      sample_names
    )

    colnames(peakDf)[
      match(sample_names, colnames(peakDf))
    ] <- raw_cols
  } else {
    peakDf[, sample_names] <- NULL
  }

  cbind(
    peakDf,
    ip_input
  )
}
