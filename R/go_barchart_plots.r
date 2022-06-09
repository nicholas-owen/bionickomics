
#' Go make those barchart plots for Webgestalt output
#'
#' @param fileList to loop over
#' @param filePath to find the files at
#'
#' @return
#' @export
#'
#' @examples


go_barchart_plots<-function(fileList, filePath){
  #create plot dimensions
  pnas.height.inches <- 225 / 25.4
  pnas.width.single.column.inches <-   87 / 25.4
  pnas.width.double.column.inches <-  178 / 25.4
  pnas.width.onehalf.column.inches <-  114 / 25.4


  for (i in 1:length(fileList)){
    dataFile<-file.path(filePath, fileList[i])
    df<-read.delim(dataFile)
    ont<-strtrim(fileList[i],width = 2)
    if (ont == "BP"){
      plotColor<-"Blues"
    } else if (ont == "MF"){
      plotColor<-"Reds"
    } else if (ont == "CC"){
      plotColor<-"Greens"
    } else if (ont == "KE"){
      plotColor<-"Oranges"
    } else if (ont == "RE"){
      plotColor<-"RdPu"
    } else if (ont == "Pa"){
      plotColor<-"PuRd"
    } else {
      plotColor<-"Greys"
    }

    p<-create_GO_barchart(df, ont, plotColor)

    #set unique filename
    file_ext <- ".svg"
    base_path <- getwd()
    dataName<-strsplit(filePath, "/")
    dataName<-unlist(dataName)[grep("_t\\d\\d$", unlist(dataName))]

    file_unique <-
      format(Sys.time(), format = "%Y-%m-%d-T%H%M%S")
    file_output <-
      file.path(base_path,
                paste0(file_unique,"_", dataName, "_GO_barchart_", ont
                       , file_ext))


    ggsave(
      filename = paste0(substr(file_output, 1, nchar(file_output) - 4),
                        ".svg"),
      plot = p,
      width = pnas.width.double.column.inches*2,
      height = pnas.height.inches ,
      dpi = 600
    )



  }

}
