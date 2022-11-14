#' Plot sample to sample distance heatmap for samples in DDS object
#'
#' @param rld - transformed dds object from DESeq2
#' @param outputFilename - filename of output image
#' @param outputFormat  - format of output image
#'
#' @return
#' @export
#'
#' @examples
plot_Sample2Sample<-function(rld, outputFilename, outputFormat="svg"){
  require("RColorBrewer")
  require("gplots")
  require("DESeq2")
  require("ggplot2")

  dfSupport<-colData(rld)
  condition<-dfSupport$condition
  mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))]
  sampleDists <- as.matrix(dist(t(assay(rld))))
  fig<-heatmap.2(as.matrix(sampleDists), key=F, trace="none",
            col=colorpanel(100, "black", "white"),
            ColSideColors=mycols[condition], RowSideColors=mycols[condition],
            margin=c(10, 10), main="Sample Distance Matrix")
  message("Saving output PCA..")

  svg(paste0(outputFilename, ".", outputFormat), w=10, h=10, pointsize=20)
  heatmap.2(as.matrix(sampleDists), key=F, trace="none",
            col=colorpanel(100, "black", "white"),
            ColSideColors=mycols[condition], RowSideColors=mycols[condition],
            margin=c(10, 10), main="Sample Distance Matrix")

  dev.off()

}
