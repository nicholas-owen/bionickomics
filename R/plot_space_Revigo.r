
#' Plots 2D semantic space of ReviGO results for gene ontology analysis
#' data format:
#' (df format TermID,	Name,	Value,	LogSize,	Frequency,	Uniqueness,	Dispensability,	PC_0,	PC_1,	Representative)
#'
#' @param revigoData - dataframe of the output table from Revigo analysis
#' @param dispens - numeric filter threshold for dispensibility
#' @param reduce - reduces the points plotted in line with dispensibility
#' @param output - filename with full path of output plot (.svg)
#' @param scale - scaling of Value; as-is (default) vs 'log'
#' @param ...

#'
#' @return
#' @export
#'
#' @examples
plot_space_Revigo<-function(revigoData, dispens, reduce=TRUE, output, scale="as-is", ...){
  library(bionickomics)
  require(ggplot2)
  require(scales)
  require(stringr)


  df <- data.frame(revigoData);
  df <- df [(df$PC_0 != "null" & df$PC_1 != "null"), ];
  df$PC_0 <- as.numeric( as.character(df$PC_0) );
  df$PC_1 <- as.numeric( as.character(df$PC_1) );
  df$LogSize <- as.numeric( as.character(df$LogSize) );
  df$Value <- as.numeric( as.character(df$Value) );
  df$Frequency <- as.numeric( as.character(df$Frequency) );
  df$Uniqueness <- as.numeric( as.character(df$Uniqueness) );
  df$Dispensability <- as.numeric( as.character(df$Dispensability) );

  if(reduce==TRUE){
   df<-subset(df, df$Dispensability<=dispens)
   message("Reducing dataframe.")
  }

  df$Name <- str_wrap(df$Name, width = 20)
  if (scale=="log"){
  df$Value<-log(-(df$Value))
  message("Rescaling to log(Value) ")
  }

  p1 <- ggplot(data = df)

  p1 <-
    p1 + geom_point(aes(PC_0, PC_1, colour = -Value, size = LogSize), alpha = I(0.6))

  if (scale=="log"){
    p1 <- p1 + scale_colour_gradientn(colours = rainbow(6), limits = c(-(max(df$Value)), 0))
  } else {
    p1 <- p1 + scale_colour_gradientn(colours = rainbow(6), limits = c(-(max(df$Value)), -(min(df$Value))))
  }




  p1 <-
    p1 + geom_point(
      aes(PC_0, PC_1, size = LogSize),
      shape = 21,
      fill = "transparent",
      colour = I (alpha ("black", 0.6))
    )
  p1 <-
    p1 + scale_size(range = c(5, 30)) + theme_bw()
  # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
  ex <- df [df$Dispensability < dispens,]

  p1 <-
    p1 + geom_text(
      data = ex,
      aes(PC_0, PC_1, label = Name),
      colour = I(alpha("black", 0.85)),
      size = 3
    )

  p1 <- p1 + labs (y = "semantic space y", x = "semantic space x")

  p1 <- p1 + theme(legend.key = element_blank())

  one.x_range = max(df$PC_0) - min(df$PC_0)

  one.y_range = max(df$PC_1) - min(df$PC_1)

  p1 <-
    p1 + xlim(min(df$PC_0) - one.x_range / 10,
              max(df$PC_0) + one.x_range / 10)

  p1 <-
    p1 + ylim(min(df$PC_1) - one.y_range / 10,
              max(df$PC_1) + one.y_range / 10)


  # --------------------------------------------------------------------------
  # Output the plot to screen



  # Uncomment the line below to also save the plot to a file.
  # The file type depends on the extension (default=pdf).

  p1<-p1 + theme_gray()

  p1+ theme( axis.text.x = element_text(color="black"),
             axis.ticks = element_line(color = "black"))

  p1<-p1+theme(axis.line.x.bottom=element_line(color="black"), axis.line.y.left=element_line(color="black"),
               axis.line.x.top=element_line(color="black"), axis.line.y.right=element_line(color="black"))

  p1<-p1+geom_hline(yintercept=00, linetype="dashed",
                    color = "#333333", size=0.5)

  p1<-p1+geom_vline(xintercept=00, linetype="dashed",
                    color = "#333333", size=0.5)

  p1

  #ggsave("./revigo-plotCC.pdf")

  pnas.height.inches<- 225/25.4
  pnas.width.single.column.inches<-87/25.4
  pnas.width.double.column.inches<-178/25.4
  pnas.width.onehalf.column.inches<-114/25.4



  #set unique filename
  # file_ext <- ".pdf"
  # base_path<-getwd()
  # file_unique <- format(Sys.time(), format = "%Y-%m-%d-T%H%M%S")
  # file_output <- file.path(paste0(base_path,"/", file_unique,"HF11_ALL_GO_Volcano_CCRevigo_LFC2"
  #                                 , file_ext))
  file_output<-output
  # ggsave(filename = paste0(substr(file_output, 1, nchar(file_output)-4),
  #                          ".svg"),plot = p1,
  #        width = pnas.width.double.column.inches, height = pnas.height.inches/1.5, dpi=600)

  ggsave(filename = output,plot = p1,
         width = pnas.width.double.column.inches*1.25, height = (pnas.height.inches/1.5)*1.25, dpi=600)
  message("Plot saved to", output)

  print(p1)


}
