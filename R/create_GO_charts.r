#' Create Gene Ontology (GO) Bar Charts
#'
#' @param df  data.frame out Webgestalt output
#' @param ont Ontology to be assessed (MF/BP/CC etc)
#' @param fill RColorBrewer colour palette for example "Reds" for further information check http://applied-r.com/rcolorbrewer-palettes/
#'
#' @return ggplot2 object of the bar chart
#' @export
#'
#' @examples

create_GO_barchart<-function(df, ont, fill){
  # df - dataframe Webgestalt output
  # ont Gene Ontology GO source , BP/MF/CC
  # fill colour of barplot using RBrewer scales check https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/

  load_pack(qvalue)
  load_pack(RColorBrewer)
  load_pack(dplyr)
  load_pack(rrvgo)
  load_pack(ggplot2)
  load_pack(stringr)

  message("Creating qvalues from Webgestalt output data.frame")
  df$qvalue <- qvalue(p = df$pValue, fdr.level = 0.05, pi0 = 1)$qvalues
  colnames(df)[1]<-"term_id"
  limitEnrich<-round(max(df$enrichmentRatio)+2)
  if (limitEnrich <=6){
    limitStep=1
  } else {
    limitStep=2
  }

  if(ont=="BP"){
    plotXTitle<-"Biological Process"
  } else if (ont == "CC"){
    plotXTitle<-"Cellular Compartment"
  } else if (ont =="MF"){
    plotXTitle<-"Molecular Function"
  } else if (ont =="KE"){
    plotXTitle<-"KEGG"
  } else if (ont =="Pa"){
    plotXTitle<-"PANTHER"
  } else if (ont =="RE"){
    plotXTitle<-"REACTOME"
  }


  df<-df %>%
    mutate(Bins = cut(pValue,
                      breaks = seq(
                        from = 0,
                        to = 0.05,
                        length.out = (dim(df)[1])+1
                      )))

  colPal<-brewer.pal(8, name = fill)
  colPal<-colorRampPalette(colPal)((dim(df)[1])*2)
  colPal<-colPal[ (2*(dim(df)[1])):(dim(df)[1]+1)]
  names(colPal) <- levels(df$Bins)
  df$Bins[is.na(df$Bins)] <- names(colPal)[1]

  if (ont =="BP" | ont =="CC" | ont == "MF"){

  simMatrix <- calculateSimMatrix(df$term_id,
                                  orgdb="org.Hs.eg.db",
                                  ont=ont,
                                  method="Rel")

  scores <- setNames(-log10(df$qvalue), df$term_id)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")
  colnames(reducedTerms)[1]<-"term_id"
  df<-merge(reducedTerms,df, by="term_id")
  # p_heatmap<-heatmapPlot(simMatrix,
  #             reducedTerms,
  #             annotateParent=TRUE,
  #             annotationLabel="parentTerm",
  #             fontsize=6)
  # p_scatter<-scatterPlot(simMatrix, reducedTerms)
  # p_treemap<-treemapPlot(reducedTerms)
  # p_wordcloud<-wordcloudPlot(reducedTerms, min.freq=1, colors="black")
  }

  p_barchart<-ggplot(df, aes(x=reorder(description, enrichmentRatio), enrichmentRatio, fill=Bins)) +
    geom_bar(stat='identity') +
    coord_flip() +
    xlab(plotXTitle) +
    ylab("Enrichment Ratio") +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14), legend.text = element_text(size=14),
          legend.title = element_text(size=14)) +
    expand_limits(x = 0, y = 0) +
    scale_y_continuous(name="Enrichment Ratio", limits=c(0, limitEnrich),
                       breaks=seq(0, limitEnrich, limitStep)) +
    scale_fill_manual(values = colPal)

  p_barchart<-p_barchart+ scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "),
                                                                         width = 30))
  # result<-list(dataframe=df, plot_heatmap=p_heatmap, plot_scatter=p_scatter, plot_treemap = p_treemap,
  #                plot_wordcloud = p_wordcloud )

  result<-p_barchart

  return(result)
}
