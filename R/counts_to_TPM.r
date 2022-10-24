##process normalised counts to TPM values

#' Title
#'
#' @param annot full annotation from GTF model (created with prepare_annotation.r)
#' @param dds
#'
#' @return
#' @export
#'
#' @examples
#' load("D:/Box Sync/Projects/CHM/RNA-seq/analysis/CHMvsWT/deseq2-HF-CHMvsRefWT-t_all-WALDp0.05/deseq2-HF-CHMvsRefWT_Rdata.rdata")
#' annot <- read.csv("D:/annotation.v104_danio_rerio/annotationDanio_rerio.GRCz11.104.gtf.annotation.full.from.GTF.csv")
#' a<-counts_to_TPM(dds, annot)

counts_to_TPM<-function(dds, annot){
  countsNorm<-DESeq2::counts(dds, normalized=TRUE)
  countsSampleNos<-dim(countsNorm)[2]

  # process annotation to get average tx length
  txData<-annot %>%
    dplyr::select(gene_name, tx_len) %>%
    group_by(gene_name) %>%
    summarise(avg=median(tx_len))
  names(txData)[2]<-"tx_len_av"

  countsNorm<-as.data.frame(countsNorm)
  countsNorm$gene_id<-rownames(countsNorm)

  # dfList <- list(countsNorm,  txData)
  # countsNormAnnot<-dfList %>%
  #   reduce(full_join, by='gene_id')
  countsNorm<-merge(countsNorm, annot, by="gene_id")
  df.counts.norm.annotated <- merge(countsNorm, txData, by="gene_name", all.x=TRUE)

    tpms <- apply(df.counts.norm.annotated[,3:(countsSampleNos+2)],2 , function(x)
    tpm(x, df.counts.norm.annotated$tx_len_av))
    tpms<-as.data.frame(tpms)

  # SANITY CHECK
  tpms.rowsums<-data.frame(colname=names(tpms), colSums_total=colSums(tpms))
  if ((sum(tpms.rowsums$colSums_total)/10^6) != countsSampleNos){
    stop("Error in SANITY check for TPMS rowsums..\n")
  }
  # Annotate tpms
  annotDim<-dim(df.counts.norm.annotated)[2]
  tpms.annotated<<-cbind(tpms, df.counts.norm.annotated [,c(1:2,(countsSampleNos+3):annotDim)])
  return(tpms.annotated)
}


##Function for TPM calculation from normalised counts
#' Title
#'
#' @param counts
#' @param lengths
#'
#' @return
#' @export
#'
#' @examples
tpm<-function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

