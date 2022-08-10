#' Plot Organoid Expression GOI: Gene of Interest
#'
#' @param goi - gene of interest as string
#' @param timepoint - day time point of data, options 'All', 'd00', 'd20', 'd35'
#' @param conditions_to_keep - list of conditions to keep in the plot ('All' or others)
#'
#' @return
#' @export
#'
#' @examples plot_organoid_goi("CRB1", "All", "PX,WT,AN,BM,UM,MY")
#'
#'

plot_organoid_goi<-function(goi, timepoint, conditions_to_keep){

  # Packages ----
  load_pack(ggplot2)
  load_pack(dplyr)
  load_pack(reshape2)


  ## OPTIONS TO SET FOR ANALYSIS
  select_time<- timepoint #"All" # options are: All, d00, d20, d35
  #goi<-"CRB1" # your gene of interest, HGNC name
  opt.cexaxis<-20
  opt.cex<-1
  if (conditions_to_keep == "All"){
    conditions_to_keep<-c("PX", "WT", "AN", "BM", "UM", "MY")
  } else {
    conditions_to_keep<-(strsplit(conditions_to_keep, split = ","))
    }

  message("Searching organoid data for : ", goi, " at timepoint : ", select_time, " Conditions : ", conditions_to_keep)
  # setup any outliers
  samples_to_remove<-c("WT02_20_1,PX01_20_4,PX01_35_4,AN02_00_3") # updated 4/8/21 for all outliers from PX01_35_4,AN02_00_3

  # Loading Data ----
  df_tpms <-  bionickomics:::tpms #read.delim("./RNA-seq-Organoid_tpms_all.tsv")
  df_tpms<-df_tpms %>%
    # select(c(-tx_len_av, -(paste0(samples_to_remove))))
    dplyr::select(-(tx_len_av))

  df_tpms_melt<-melt(df_tpms)
  df_tpms_melt$time<-as.factor(paste0("d",substring(df_tpms_melt$variable, 6,7)))
  df_tpms_melt$variable<-strtrim(df_tpms_melt$variable,2)


  #keep only certain conditions
  df_tpms_melt<- df_tpms_melt[df_tpms_melt$variable %in% conditions_to_keep[[1]], ]

  p <-ggplot(
    {if (select_time == "All") {
      data=df_tpms_melt[df_tpms_melt$gene_name == goi,]
    } else if (select_time == "d00") {
      data=subset(df_tpms_melt[df_tpms_melt$gene_name == goi,], df_tpms_melt[df_tpms_melt$gene_name == goi,]$time=="d00")
    } else if (select_time == "d20") {
      data=subset(df_tpms_melt[df_tpms_melt$gene_name == goi,],df_tpms_melt[df_tpms_melt$gene_name == goi,]$time=="d20")
    } else if (select_time == "d35") {
      data=subset(df_tpms_melt[df_tpms_melt$gene_name == goi,], df_tpms_melt[df_tpms_melt$gene_name == goi,]$time=="d35")
    }}
    , aes(x = variable, fill=time, y=value)) +
    geom_boxplot(colour = "black") + geom_jitter(width=0.2, size=opt.cex, color="black") +
    xlab("Condition") +
    ylab("TPM") +
    ggtitle(paste0(goi," Expression plot")) +
    expand_limits(y = 0) +
    theme(legend.position="bottom") +
    theme(text = element_text(size = opt.cexaxis))


  fileOuput<-paste0("TPM_exp_", goi, "_time-", select_time, ".svg")
  message("Saving : ", fileOuput)
  ggsave(filename = fileOuput, plot = p, width = 16, height = 12)
  message("Returning plot..")
  print(p)

}


