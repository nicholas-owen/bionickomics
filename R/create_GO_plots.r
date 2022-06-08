


#' Cleans Gene Ontology outputs and combines with ReviGO reduced data
#'
#' @param results
#' @param revigo
#'
#' @return df
#' @export
#'
#' @examples
clean_GO_output <- function(results, revigo) {
  colnames(results)[1] <- "geneSet"
  df <- merge(results, revigo)
  df <- subset(df, df$Eliminated == "False")
  df$overlapscaled <- (df$overlap) / (max(df$overlap))
  top.GO.number <- dim(df)[1]
  #set max to 16 for -log10pval
  df$logp <- pmin((-(log10(df$pValue))), 16)
  df$description <- str_wrap(df$description, width = 20)
  return(df)
}





#' Creates Gene Ontology plots for Enrichment
#'
#' @param df combined output from clean_GO_output
#' @param ont Ontology type BP/CC/MF/KEGG etc
#' @param fill Plot colour
#'
#' @return p
#' @export
#'
#' @examples
make_GO_plots <- function(df, ont, fill) {
  p <-
    ggplot(df, aes(
      x = enrichmentRatio,
      y = -log10(pValue),
      label = description
    )) +
    geom_point(
      stat = 'identity',
      size = 1 + (15 * df$overlapscaled),
      color = "black",
      pch = 21,
      fill = fill
    ) +
    xlab(expression(" Enrichment Ratio")) +
    ylab(expression(-Log["10"](p - value))) +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14)
    ) +
    #expand_limits(x = 0, y = 0) +
    geom_text_repel(size = 5,
                    #geom_label_repel if boxes required
                    data = df[1:top.GO.number, ],
                    ##use 1:topgenenumber of dataframe
                    #fontface = "italic",
                    #arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "last"),
                    force = 5)

  p <-
    p + scale_y_continuous(
      limits = c(0, max(df$logp)+1),
      expand = c(0, 0),
      breaks = seq(0, max(df$logp)+1, by = 2)
    )
  p <-
    p + scale_x_continuous(limits = c(
      0,
      round(max(df$enrichmentRatio) + 1),
      expand = c(0, 0),
      breaks = seq(0, round(max(
        df$enrichmentRatio
      ) + 1),
      by =
        1)
    ))
  p <-
    p + theme_bw() + theme(legend.position = "none")

  p

  #set unique filename
  file_ext <- ".pdf"
  base_path <- getwd()
  file_unique <-
    format(Sys.time(), format = "%Y-%m-%d-T%H%M%S")
  file_output <-
    file.path(base_path,
              paste0(file_unique, "_GO_Volcano_", ont
                     , file_ext))
  pnas.height.inches <- 225 / 25.4
  pnas.width.single.column.inches <-
    87 / 25.4
  pnas.width.double.column.inches <-
    178 / 25.4
  pnas.width.onehalf.column.inches <-
    114 / 25.4



  ggsave(
    filename = paste0(substr(file_output, 1, nchar(file_output) - 4),
                      ".svg"),
    plot = p,
    width = pnas.width.double.column.inches,
    height = pnas.height.inches / 1.5,
    dpi = 600
  )


  ###without labels for publication

  p2 <-
    ggplot(df, aes(
      x = enrichmentRatio,
      y = -log10(pValue),
      label = description
    )) +
    geom_point(
      stat = 'identity',
      size = 1 + (15 * df$overlapscaled),
      color = "black",
      pch = 21,
      fill = fill
    ) +
    xlab(expression(" Enrichment Ratio")) +
    ylab(expression(-Log["10"](p - value))) +
    theme(
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 18),
      legend.text = element_text(size = 14),
      legend.title = element_text(size =
                                    14)
    )# +
  #expand_limits(x = 0, y = 0) +
  #geom_text_repel( size=5, #geom_label_repel if boxes required
  #                data=df.final[1:top.GO.number,], ##use 1:topgenenumber of dataframe
  #fontface = "italic",
  #arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "last"),
  #                force = 5
  #)


  p2 <-
    p2 + scale_y_continuous(
      limits = c(0, max(df$logp)+1),
      expand = c(0, 0),
      breaks = seq(0, max(df$logp)+1, by = 2)
    )
  p <-
    p + scale_x_continuous(limits = c(
      0,
      round(max(df$enrichmentRatio) + 1),
      expand = c(0, 0),
      breaks = seq(0, round(max(
        df$enrichmentRatio
      ) + 1),
      by =
        1)
    ))
  p2 <-
    p2 + theme_bw() + theme(legend.position = "none")


  ggsave(
    filename = paste0(substr(file_output, 1, nchar(file_output) - 4),
                      "_nolabels.svg"),
    plot = p2,
    width = pnas.width.double.column.inches,
    height = pnas.height.inches / 1.5,
    dpi = 600
  )





  return(p)

}
