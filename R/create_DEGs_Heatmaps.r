

#' Create Heatmaps of Differentially Expressed Genes
#'
#' input can be from pairwise analysis of grouped - see keycode:extra for Differential Expression Analysis markdown document.
#'
#' @param dataFolder Location of .RData file from Differential Gene Expression analysis (absolute path)
#' @param geneList a list object of genes of interest, either Ensembl ID or Hugo gene names, or "" or blank
#' @param outputFileName  string of characters for file output
#'
#' @return outputs files as PDF for important as needed
#' @export
#'
#' @examples dataFolder<-"D:/Cloud/GoogleDrive/Work UCL/Documents/Manuscripts/Shared_Pathways_MAC/data/BM-PX-WT_t35/"
#' @examples outputFileName<-"testing"
#' @examples
#' @examples create_DEGs_Heatmaps(dataFolder=dataFolder, geneList = "", outputFileName = outputFileName)
#'


create_DEGs_Heatmaps<-function(dataFolder, geneList="", outputFileName){

  `%notin%` <- Negate(`%in%`)

  dataFile<-list.files(dataFolder, pattern = "*.rdata", recursive = TRUE)
  proteins_only<-"Y"

  message("Loading RData objects.")
  load(paste0(dataFolder, dataFile))
  message("Loading R libraries")
  library(bionickomics)
  load_pack(readr)
  load_pack(dplyr)
  load_pack(ggplot2)
  load_pack(DESeq2)
  load_pack(NMF)
  load_pack(dendextend)
  load_pack(RColorBrewer)
  load_pack(dichromat)
  load_pack(DESeq2)
  load_pack(readxl)
  load_pack(pheatmap)

  message("DDS object includes:")
  cat(resultsNames(dds))

  # dfConditions<-expand.grid(test=conditions.to.keep, ref=conditions.to.keep)
  # dfConditions<-crossing( conditions.to.keep[1:(length(conditions.to.keep))-1], conditions.to.keep)
  if (length(resultsNames(dds))>2){
    dfConditions<-combn(conditions.to.keep, 2, simplify = TRUE)
    row.names(dfConditions)<-c("Test", "Ref")

    for (i in 1:dim(dfConditions)[2]){

      #get first character of each group condition name:
      condition.matrix.reference<-as.character(dfConditions[2,i])
      condition.matrix.test<-as.character(dfConditions[1,i])


      if (!(condition.matrix.reference == condition.matrix.test) & condition.matrix.test == "WT"){

        condition.matrix.ref.temp<-(dfConditions[2,i])
        dfConditions[2,i]<-dfConditions[1,i]
        dfConditions[1,i]<-condition.matrix.ref.temp

      }

    }
    #dfConditions<-as.data.frame(t(dfConditions))
    conditions.to.test.total<-dim(dfConditions)[[2]]
  } else {
        dfConditions<-t(data.frame(Test=(levels(condition)[1]), Ref=(levels(condition)[2])))
        conditions.to.test.total<-dim(dfConditions)[[2]]
      }



  dfConditions<-as.data.frame(t(dfConditions))

  for (i in 1:dim(dfConditions)[1]){
    loopTest<-dfConditions[i,1]
    loopRef<-dfConditions[i,2]
    cat("Results: ", loopTest, "-",loopRef, "\n")

    res<-results(dds, contrast = c("condition",loopTest,loopRef))
    results<-as.data.frame(res)
    # Order by adjusted p-value
    results <- results[order(results$padj), ]
    ## Merge with normalized count data
    resdata <- merge(as.data.frame(results), as.data.frame(counts(dds, normalized=TRUE)),
                     by="row.names", sort=FALSE)
    names(resdata)[1] <- "gene_id"
    message("Head of the results data object.")
    head(resdata, n=1)
    ##merge with results
    results<-merge(resdata, annotation, by='gene_id', all.x = TRUE)
    message("Filtering results.")
    results<-results %>%
      select(gene_id, gene_name, gene_biotype, everything()) %>%
      arrange(padj)

    filt.padj<-0.05
    filt.lfc<-1
    results.filtered<-subset(results, (results$baseMean >=10) & (results$padj<=filt.padj) &
                               (abs(results$log2FoldChange))>=filt.lfc)

    if (proteins_only=="Y"){
      results.filtered<-subset(results, results$gene_biotype=="protein_coding")
      message("Keeping biotype=protein only.")
    }

    if (geneList!=""){
    results.filtered<-results.filtered[results.filtered$gene_name %in% geneList,]
    message("Filtering on user specified gene list.")
    top.hits<-dim(results.filtered)[[1]]
    message("Number of observations: ", top.hits)
    } else {

      top.hits<-100
      message("No gene list supplied. Fixing top DEG list to 100")
    }



    dds.rld.transformed<-rlog(dds, blind=FALSE) ###BLIND tested for TRUE also

    sample.names<-colnames(dds)

    heatmap.data.rld.assay<-assay(dds.rld.transformed)[arrange(results.filtered, padj, pvalue)$gene_id[1:top.hits],]
    samples.to.keep<-grep(paste0(loopTest,"|", loopRef), colnames(heatmap.data.rld.assay))

    heatmap.annotation.columns<-cbind.data.frame(Condition = colData(dds)$condition)
    heatmap.data.rld.assay.subset<-heatmap.data.rld.assay[, samples.to.keep]
    heatmap.data.rld.assay.subset<-as.matrix(heatmap.data.rld.assay.subset)

    heatmap.annotation.columns.subset<-droplevels(heatmap.annotation.columns[samples.to.keep,])

    ntd <- normTransform(dds)

    # select <- order(rowMeans(counts(dds,normalized=TRUE)),
    #                 decreasing=TRUE)[1:20]
    # df <- as.data.frame(colData(dds)[,c("condition")])
    # pheatmap(assay(rld)[select,], cluster_rows=TRUE, show_rownames=TRUE,
    #          cluster_cols=TRUE,
    #          scale = "row",
    #          #annotation_row=toupper(arrange(results.filtered, padj,pvalue)$gene_name[1:20]),
    #          #annotation_col = sample.names
    #          cutree_rows = 10,
    #          cutree_cols = 3
    # )


    # load_pack(RColorBrewer)
    colrs=colorRampPalette(rev(brewer.pal(3,"Set1")))(100)
    colrs=colorRampPalette(c( "blue","white", "orange"))(100)

    filt.padj<-0.05
    filt.lfc<-1

    colrs = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255)

    colrs=colorRampPalette(rev(brewer.pal(3,"Set1")))(100)
    colrs=colorRampPalette(c( "blue","white", "orange"))(100)

    heatmap.data.rld.assay<-assay(dds.rld.transformed)[arrange(results.filtered, padj, pvalue)$gene_id[1:top.hits],]
    heatmap.annotation.columns<-cbind.data.frame(Tissue = colData(dds)$condition)

    samples.to.plot<-colnames(dds)
    samples.to.keep<-grep(paste0(loopTest,"|", loopRef), colnames(heatmap.data.rld.assay))
    sample.names<-colnames(dds)[samples.to.keep]
    head(heatmap.data.rld.assay[,samples.to.keep])

    heatmap.data.rld.assay.subset<-heatmap.data.rld.assay[, samples.to.keep]
    heatmap.data.rld.assay.subset<-as.matrix(heatmap.data.rld.assay.subset)
    heatmap.annotation.columns.subset<-droplevels(heatmap.annotation.columns[samples.to.keep,])
    # class(heatmap.data.rld.assay.subset[arrange(results.filtered, padj, pvalue)$gene_id[1:top.hits], ])
    heatmap.data<-heatmap.data.rld.assay.subset[arrange(results.filtered, padj, pvalue)$gene_id[1:top.hits], ]
    rownames(heatmap.data)<-tolower(arrange(results.filtered, padj,pvalue)$gene_name[1:top.hits])
    my.sample.cols<-data.frame(Tissue = heatmap.annotation.columns.subset)
    rownames(my.sample.cols)<-sample.names
    my.sample.cols$Tissue<-as.factor(my.sample.cols$Tissue)

    #my_colour = list(
    #  Tissue = c(WT = "#000000", AN = "#555555", PX= "#808080")#
    #
    #)

    my_colour = list(
      Tissue = c(a = "#000000", b = "#555555", c="#EEEEEE", d="#DDDDDD", e="#AAAAAA", f="#BBBBBB")#
    )

    my_colour$Tissue<-sample(my_colour$Tissue, size = 2, replace = FALSE)
    names(my_colour$Tissue)<-c(loopRef, loopTest)

    fontsize_row = 10 - nrow(heatmap.data) / 15


    #set unique filename
    file_ext <- ".pdf"
    base_path<-getwd()
    file_unique <- format(Sys.time(), format = "%Y-%m-%d-T%H%M%S")
    file_output <- file.path(dataFolder, paste0(file_unique,"_", loopTest,
                                                 "vs", loopRef,"_", outputFileName, file_ext))


    graphics.off()
    message("Creating output PDF of heatmap: ", loopTest, " vs ", loopRef)
    message("Saving to file: ", file_output)
    pdf(file_output, paper="special", width=pnas.width.onehalf.column.inches, height=pnas.height.inches)

    rownames(heatmap.data)<-toupper(rownames(heatmap.data))
    pheatmap(heatmap.data, cluster_rows=TRUE, show_rownames=TRUE,
             cluster_cols=TRUE,
             annotation_colors = my_colour,
             scale = "row",
             #annotation_row=toupper(arrange(results.filtered, padj,pvalue)$gene_name[1:20]),
             #annotation_col = sample.names
             #cutree_rows = 10,
             # cutree_cols = 3,
             color=colrs,
             annotation_col = my.sample.cols,
             fontsize_row=fontsize_row,
             border_color = NA
    )


    dev.off()
    graphics.off()


  }

}



