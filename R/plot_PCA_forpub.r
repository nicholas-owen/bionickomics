

#' Plot PCA for publication format
#'
#' @param rlogdata - object of transformed data from DESeq2
#' @param samplesToRemove  - list object of samples to remove from dataset, default="NONE")
#' @param output  - filename of the output in cwd
#' @param outFormat  - output format of the plot (default = SVG)
#'
#' @return
#' @export
#'
#' @examples
plot_PCA_forpub <- function(rlogdata, samplesToRemove = "NONE",
                            output, outFormat = "svg") {
  require("PCAtools")
  require("DESeq2")
  require("ggplot2")
  require("ggalt")
  #input will be the rld object from transformed DDS object

  #reordering based on alpha
  rlogcounts<-assay(rlogdata)
  #rlogcounts<-rlogcounts[,order(colnames(rlogcounts))]
  rlogcounts<-rlogcounts[, !colnames(rlogcounts) %in% samplesToRemove]
  dfMeta<-colData(rld)
  dfMeta<-dfMeta[!(row.names(dfMeta) %in% samplesToRemove), ]

  p <- pca(rlogcounts, metadata = dfMeta, removeVar = 0.1)
  fig <- biplot_omics(
    p,
    colby = 'condition',
    colLegendTitle = 'Condition',
    # encircle config
    encircle = TRUE,
    encircleFill = TRUE,
    hline = 0,
    vline = c(-25, 0, 25),
    legendPosition = 'bottom',
    legendLabSize = 12,
    legendIconSize = 7.0
  )
message("Saving output PCA..")
  ggsave(
    filename = paste0(output, ".", outFormat),
    plot = fig,
    device=outFormat,
    width = 16,
    height = 12,
    dpi = 600
  )





}
