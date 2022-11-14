#' PCA biplot
#'
#' @param pcaobj
#' @param x
#' @param y
#' @param showLoadings
#' @param ntopLoadings
#' @param showLoadingsNames
#' @param colLoadingsNames
#' @param sizeLoadingsNames
#' @param boxedLoadingsNames
#' @param fillBoxedLoadings
#' @param drawConnectorsLoadings
#' @param widthConnectorsLoadings
#' @param colConnectorsLoadings
#' @param lengthLoadingsArrowsFactor
#' @param colLoadingsArrows
#' @param widthLoadingsArrows
#' @param alphaLoadingsArrow
#' @param colby
#' @param colkey
#' @param colLegendTitle
#' @param singlecol
#' @param shape
#' @param shapekey
#' @param shapeLegendTitle
#' @param pointSize
#' @param legendPosition
#' @param legendLabSize
#' @param legendTitleSize
#' @param legendIconSize
#' @param encircle
#' @param encircleFill
#' @param encircleFillKey
#' @param encircleAlpha
#' @param encircleLineSize
#' @param encircleLineCol
#' @param ellipse
#' @param ellipseConf
#' @param ellipseFill
#' @param ellipseFillKey
#' @param ellipseAlpha
#' @param ellipseLineSize
#' @param ellipseLineCol
#' @param xlim
#' @param ylim
#' @param lab
#' @param labSize
#' @param labhjust
#' @param labvjust
#' @param boxedLabels
#' @param selectLab
#' @param drawConnectors
#' @param widthConnectors
#' @param colConnectors
#' @param xlab
#' @param xlabAngle
#' @param xlabhjust
#' @param xlabvjust
#' @param ylab
#' @param ylabAngle
#' @param ylabhjust
#' @param ylabvjust
#' @param axisLabSize
#' @param title
#' @param subtitle
#' @param caption
#' @param titleLabSize
#' @param subtitleLabSize
#' @param captionLabSize
#' @param hline
#' @param hlineType
#' @param hlineCol
#' @param hlineWidth
#' @param vline
#' @param vlineType
#' @param vlineCol
#' @param vlineWidth
#' @param gridlines.major
#' @param gridlines.minor
#' @param borderWidth
#' @param borderColour
#' @param returnPlot
#'
#' @return
#' @export
#'
#' @examples
biplot_omics<-function (pcaobj, x = "PC1", y = "PC2", showLoadings = FALSE,
                       ntopLoadings = 5, showLoadingsNames = if (showLoadings) TRUE else FALSE,
                       colLoadingsNames = "black", sizeLoadingsNames = 3,
                       boxedLoadingsNames = TRUE, fillBoxedLoadings = alpha("white",
                                                                            1/4), drawConnectorsLoadings = TRUE, widthConnectorsLoadings = 0.5,
                       colConnectorsLoadings = "grey50", lengthLoadingsArrowsFactor = 1.5,
                       colLoadingsArrows = "black", widthLoadingsArrows = 0.5,
                       alphaLoadingsArrow = 1, colby = NULL, colkey = NULL, colLegendTitle = if (!is.null(colby)) colby else NULL,
                       singlecol = NULL, shape = NULL, shapekey = NULL, shapeLegendTitle = if (!is.null(shape)) shape else NULL,
                       pointSize = 3, legendPosition = "none", legendLabSize = 12,
                       legendTitleSize = 14, legendIconSize = 5, encircle = FALSE,
                       encircleFill = TRUE, encircleFillKey = NULL, encircleAlpha = 1/4,
                       encircleLineSize = 0.25, encircleLineCol = NULL, ellipse = FALSE,
                       ellipseConf = 0.95, ellipseFill = TRUE, ellipseFillKey = NULL,
                       ellipseAlpha = 1/4, ellipseLineSize = 0.25, ellipseLineCol = NULL,
                       xlim = if (showLoadings) c(min(pcaobj$rotated[, x]) - 5,
                                                  max(pcaobj$rotated[, x]) + 5) else c(min(pcaobj$rotated[,
                                                                                                          x]) - 1, max(pcaobj$rotated[, x]) + 1), ylim = if (showLoadings) c(min(pcaobj$rotated[,
                                                                                                                                                                                                y]) - 5, max(pcaobj$rotated[, y]) + 5) else c(min(pcaobj$rotated[,
                                                                                                                                                                                                                                                                 y]) - 1, max(pcaobj$rotated[, y]) + 1), lab = rownames(pcaobj$metadata),
                       labSize = 3, labhjust = 1.5, labvjust = 0, boxedLabels = FALSE,
                       selectLab = NULL, drawConnectors = TRUE, widthConnectors = 0.5,
                       colConnectors = "grey50", xlab = paste0(x, ", ",
                                                               round(pcaobj$variance[x], digits = 2), "% variation"),
                       xlabAngle = 0, xlabhjust = 0.5, xlabvjust = 0.5, ylab = paste0(y,
                                                                                      ", ", round(pcaobj$variance[y], digits = 2), "% variation"),
                       ylabAngle = 0, ylabhjust = 0.5, ylabvjust = 0.5, axisLabSize = 16,
                       title = "", subtitle = "", caption = "",
                       titleLabSize = 16, subtitleLabSize = 12, captionLabSize = 12,
                       hline = NULL, hlineType = "longdash", hlineCol = "black",
                       hlineWidth = 0.4, vline = NULL, vlineType = "longdash",
                       vlineCol = "black", vlineWidth = 0.4, gridlines.major = TRUE,
                       gridlines.minor = TRUE, borderWidth = 0.8, borderColour = "black",
                       returnPlot = TRUE)
{
  labFun <- xidx <- yidx <- NULL
  th <- theme_bw(base_size = 24) + theme(legend.background = element_rect(),
                                         plot.title = element_text(angle = 0, size = titleLabSize,
                                                                   face = "bold", vjust = 1), plot.subtitle = element_text(angle = 0,
                                                                                                                           size = subtitleLabSize, face = "plain", vjust = 1),
                                         plot.caption = element_text(angle = 0, size = captionLabSize,
                                                                     face = "plain", vjust = 1), axis.text.x = element_text(angle = xlabAngle,
                                                                                                                            size = axisLabSize, hjust = xlabhjust, vjust = xlabvjust),
                                         axis.text.y = element_text(angle = ylabAngle, size = axisLabSize,
                                                                    hjust = ylabhjust, vjust = ylabvjust), axis.title = element_text(size = axisLabSize),
                                         legend.position = legendPosition, legend.key = element_blank(),
                                         legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = legendLabSize),
                                         title = element_text(size = legendLabSize), legend.title = element_text(size = legendTitleSize))
  plotobj <- NULL
  plotobj$x <- pcaobj$rotated[, x]
  plotobj$y <- pcaobj$rotated[, y]
  if (!is.null(lab)) {
    plotobj$lab <- lab
  }
  plotobj <- as.data.frame(plotobj, stringsAsFactors = FALSE)
  if (!is.null(selectLab)) {
    if (is.null(lab)) {
      stop(paste0("You have specified lab as NULL ",
                  "- no labels can be selected!"))
    }
    else {
      names.new <- rep(NA, length(plotobj$lab))
      indices <- which(plotobj$lab %in% selectLab)
      names.new[indices] <- plotobj$lab[indices]
      plotobj$lab <- names.new
    }
  }
  if (is.null(colby)) {
    if (!is.null(lab)) {
      plotobj$col <- lab
    }
    else {
      plotobj$col <- seq_len(length(pcaobj$yvars))
    }
  }
  else {
    plotobj$col <- pcaobj$metadata[, colby]
  }
  if (!is.null(shape)) {
    plotobj$shape <- pcaobj$metadata[, shape]
  }
  plot <- ggplot(plotobj, aes(x = x, y = y)) + th + guides(fill = guide_legend(),
                                                           shape = guide_legend(), colour = guide_legend(override.aes = list(size = legendIconSize)))
  if (is.null(singlecol)) {
    if (!is.null(shape)) {
      plot <- plot + geom_point(aes(color = col, shape = shape),
                                size = pointSize)
    }
    else {
      plot <- plot + geom_point(aes(color = col), size = pointSize)
    }
  }
  else if (!is.null(singlecol)) {
    if (!is.null(shape)) {
      plot <- plot + geom_point(aes(color = singlecol,
                                    shape = shape), size = pointSize)
    }
    else {
      plot <- plot + geom_point(aes(color = singlecol),
                                size = pointSize)
    }
  }
  if (!is.null(colkey)) {
    plot <- plot + scale_colour_discrete("") + scale_color_manual(values = colkey)
  }
  if (!is.null(shapekey)) {
    plot <- plot + scale_shape_manual(values = shapekey)
  }
  if (showLoadings) {
    xidx <- order(abs(pcaobj$loadings[, x]), decreasing = TRUE)
    yidx <- order(abs(pcaobj$loadings[, y]), decreasing = TRUE)
    vars <- unique(c(rownames(pcaobj$loadings)[xidx][seq_len(ntopLoadings)],
                     rownames(pcaobj$loadings)[yidx][seq_len(ntopLoadings)]))
    r <- min((max(pcaobj$rotated[, x]) - min(pcaobj$rotated[,
                                                            x])/(max(pcaobj$loadings[, x]) - min(pcaobj$loadings[,
                                                                                                                 x]))), (max(pcaobj$rotated[, y]) - min(pcaobj$rotated[,
                                                                                                                                                                       y])/(max(pcaobj$loadings[, y]) - min(pcaobj$loadings[,
                                                                                                                                                                                                                            y]))))
    plot <- plot + geom_segment(data = pcaobj$loadings[vars,
    ], aes(x = 0, y = 0, xend = pcaobj$loadings[vars,
                                                x] * r * lengthLoadingsArrowsFactor, yend = pcaobj$loadings[vars,
                                                                                                            y] * r * lengthLoadingsArrowsFactor), arrow = arrow(length = unit(1/2,
                                                                                                                                                                              "picas"), ends = "last"), color = colLoadingsArrows,
    size = widthLoadingsArrows, alpha = alphaLoadingsArrow,
    show.legend = NA)
    if (showLoadingsNames) {
      if (drawConnectorsLoadings) {
        if (boxedLoadingsNames) {
          plot <- plot + coord_equal() + geom_label_repel(max.overlaps=Inf, data = pcaobj$loadings[vars,
          ], aes(label = vars, x = pcaobj$loadings[vars,
                                                   x] * r * lengthLoadingsArrowsFactor, y = pcaobj$loadings[vars,
                                                                                                            y] * r * lengthLoadingsArrowsFactor, hjust = 0),
          color = colLoadingsNames, size = sizeLoadingsNames,
          fill = fillBoxedLoadings, segment.color = colConnectorsLoadings,
          segment.size = widthConnectorsLoadings)
        }
        else {
          plot <- plot + coord_equal() + geom_text_repel(max.overlaps=Inf, data = pcaobj$loadings[vars,
          ], aes(label = vars, x = pcaobj$loadings[vars,
                                                   x] * r * lengthLoadingsArrowsFactor, y = pcaobj$loadings[vars,
                                                                                                            y] * r * lengthLoadingsArrowsFactor, hjust = 0),
          color = colLoadingsNames, size = sizeLoadingsNames,
          segment.color = colConnectorsLoadings, segment.size = widthConnectorsLoadings)
        }
      }
      else {
        if (boxedLoadingsNames) {
          plot <- plot + coord_equal() + geom_label(data = pcaobj$loadings[vars,
          ], aes(label = vars, x = pcaobj$loadings[vars,
                                                   x] * r * lengthLoadingsArrowsFactor, y = pcaobj$loadings[vars,
                                                                                                            y] * r * lengthLoadingsArrowsFactor, hjust = 0),
          color = colLoadingsNames, size = sizeLoadingsNames,
          fill = NA)
        }
        else {
          plot <- plot + coord_equal() + geom_text(data = pcaobj$loadings[vars,
          ], aes(label = vars, x = pcaobj$loadings[vars,
                                                   x] * r * lengthLoadingsArrowsFactor, y = pcaobj$loadings[vars,
                                                                                                            y] * r * lengthLoadingsArrowsFactor, hjust = 0),
          color = colLoadingsNames, size = sizeLoadingsNames,
          check_overlap = TRUE)
        }
      }
    }
  }
  plot <- plot + xlab(xlab) + ylab(ylab)
  if (!is.null(xlim)) {
    plot <- plot + xlim(xlim[1], xlim[2])
  }
  if (!is.null(ylim)) {
    plot <- plot + ylim(ylim[1], ylim[2])
  }
  plot <- plot + labs(title = title, subtitle = subtitle, caption = caption,
                      fill = "", colour = colLegendTitle, shape = shapeLegendTitle)
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline, linetype = vlineType,
                              colour = vlineCol, size = vlineWidth)
  }
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = hline, linetype = hlineType,
                              colour = hlineCol, size = hlineWidth)
  }
  plot <- plot + theme(panel.border = element_rect(colour = borderColour,
                                                   fill = NA, size = borderWidth))
  if (gridlines.major == TRUE) {
    plot <- plot + theme(panel.grid.major = element_line())
  }
  else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor == TRUE) {
    plot <- plot + theme(panel.grid.minor = element_line())
  }
  else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }
  if (boxedLabels) {
    if (drawConnectors) {
      labFun <- function(...) geom_label_repel(max.overlaps=Inf,...)
    }
    else {
      labFun <- function(...) geom_label(...)
    }
  }
  else {
    if (drawConnectors) {
      labFun <- function(...) geom_text_repel(max.overlaps=Inf,...)
    }
    else {
      labFun <- function(...) geom_text(...)
    }
  }
  if (!is.null(lab)) {
    if (drawConnectors && is.null(selectLab)) {
      plot <- plot + labFun(data = plotobj, aes(label = lab),
                            size = labSize, segment.color = colConnectors,
                            segment.size = widthConnectors, hjust = labhjust,
                            vjust = labvjust)
    }
    else if (drawConnectors && !is.null(selectLab)) {
      plot <- plot + labFun(data = subset(plotobj, !is.na(plotobj[,
                                                                  "lab"])), aes(label = lab), size = labSize,
                            segment.color = colConnectors, segment.size = widthConnectors,
                            hjust = labhjust, vjust = labvjust)
    }
    else if (!drawConnectors && !is.null(selectLab)) {
      if (boxedLabels) {
        plot <- plot + labFun(data = subset(plotobj,
                                            !is.na(plotobj[, "lab"])), aes(label = lab),
                              size = labSize, hjust = labhjust, vjust = labvjust)
      }
      else {
        plot <- plot + labFun(data = subset(plotobj,
                                            !is.na(plotobj[, "lab"])), aes(label = lab),
                              size = labSize, check_overlap = TRUE, hjust = labhjust,
                              vjust = labvjust)
      }
    }
    else if (!drawConnectors && is.null(selectLab)) {
      if (boxedLabels) {
        plot <- plot + labFun(data = plotobj, aes(label = lab),
                              size = labSize, check_overlap = TRUE, hjust = labhjust,
                              vjust = labvjust)
      }
      else {
        plot <- plot + labFun(data = plotobj, aes(label = lab),
                              size = labSize, check_overlap = TRUE, hjust = labhjust,
                              vjust = labvjust)
      }
    }
  }
  if (encircle) {
    if (encircleFill) {
      if (is.null(encircleLineCol)) {
        plot <- plot + ggalt::geom_encircle(aes(group = col,
                                                fill = col, colour = col), alpha = encircleAlpha,
                                            size = encircleLineSize, show.legend = FALSE,
                                            na.rm = TRUE)
      }
      else {
        plot <- plot + ggalt::geom_encircle(aes(group = col,
                                                fill = col), colour = encircleLineCol, alpha = encircleAlpha,
                                            size = encircleLineSize, show.legend = FALSE,
                                            na.rm = TRUE)
      }
    }
    else {
      if (is.null(encircleLineCol)) {
        plot <- plot + ggalt::geom_encircle(aes(group = col,
                                                colour = col), fill = NA, alpha = encircleAlpha,
                                            size = encircleLineSize, show.legend = FALSE,
                                            na.rm = TRUE)
      }
      else {
        plot <- plot + ggalt::geom_encircle(aes(group = col),
                                            colour = encircleLineCol, fill = NA, alpha = encircleAlpha,
                                            size = encircleLineSize, show.legend = FALSE,
                                            na.rm = TRUE)
      }
    }
    if (encircleFill) {
      if (is.null(encircleFillKey)) {
        if (!is.null(colkey)) {
          plot <- plot + scale_fill_manual(values = colkey)
        }
      }
      else {
        plot <- plot + scale_fill_manual(values = encircleFillKey)
      }
    }
  }
  if (ellipse) {
    if (ellipseFill) {
      if (is.null(ellipseLineCol)) {
        plot <- plot + stat_ellipse(aes(group = col,
                                        fill = col, colour = col), geom = "polygon",
                                    level = ellipseConf, alpha = ellipseAlpha,
                                    size = ellipseLineSize, show.legend = FALSE,
                                    na.rm = TRUE)
      }
      else {
        plot <- plot + stat_ellipse(aes(group = col,
                                        fill = col), colour = ellipseLineCol, geom = "polygon",
                                    level = ellipseConf, alpha = ellipseAlpha,
                                    size = ellipseLineSize, show.legend = FALSE,
                                    na.rm = TRUE)
      }
    }
    else {
      if (is.null(ellipseLineCol)) {
        plot <- plot + stat_ellipse(aes(group = col,
                                        colour = col), fill = NA, geom = "polygon",
                                    level = ellipseConf, alpha = ellipseAlpha,
                                    size = ellipseLineSize, show.legend = FALSE,
                                    na.rm = TRUE)
      }
      else {
        plot <- plot + stat_ellipse(aes(group = col),
                                    colour = ellipseLineCol, fill = NA, geom = "polygon",
                                    level = ellipseConf, alpha = ellipseAlpha,
                                    size = ellipseLineSize, show.legend = FALSE,
                                    na.rm = TRUE)
      }
    }
    if (ellipseFill) {
      if (is.null(ellipseFillKey)) {
        if (!is.null(colkey)) {
          plot <- plot + scale_fill_manual(values = colkey)
        }
      }
      else {
        plot <- plot + scale_fill_manual(values = ellipseFillKey)
      }
    }
  }
  if (returnPlot) {
    return(plot)
  }
  else if (!returnPlot) {
    plot
  }
}
