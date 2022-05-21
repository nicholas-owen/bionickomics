
# https://github.com/skurscheid/deepToolsUtils/blob/master/R/preProcessingEnsembl.R
#' Check Ensembl for available releases
#'
#' @param url
#'
#' @return
#' @export
#'
#' @examples
CheckAvailableEnsemblReleases <- function(url = "ftp://ftp.ensembl.org/pub/"){
  releases <- getURL(url = url, dirlistonly = TRUE)
  releases <- strsplit(releases, "\r*\n")[[1]]
  releases <- releases[grep("release", releases)]
  return(releases)
}
#' Check Ensembl for available Organisms
#'
#' @param url
#' @param release
#'
#' @return
#' @export
#'
#' @examples
CheckAvailableEnsemblOrganisms <- function(url = "ftp://ftp.ensembl.org/pub/", release = NULL){
  if (is.null(release)) stop("Ensembl release number missing. Please check with CheckAvailableEnsemblReleases() and specify in function call.")
  tryCatch(match.arg(release, choices = CheckAvailableEnsemblReleases()))
  organisms <- getURL(url = paste(url, release, "/fasta/", sep = ""), dirlistonly = TRUE)
  organisms <- strsplit(organisms, "\r*\n")[[1]]
  return(organisms)
}

#' Download Ensembl GTF
#'
#' @param version
#' @param url
#' @param temp_dir
#' @param organism
#' @param data_type
#'
#' @return
#' @export
#'
#' @examples

DownloadEnsemblGTF <- function(version = NULL,
                               url = "ftp.ensembl.org/pub/",
                               temp_dir = NULL,
                               organism = NULL,
                               data_type = "genes"){
  match.arg(organism, choices = CheckAvailableEnsemblOrganisms(release = version))
  #match.arg(data_type, choices = c())

  url<-paste0(url, version,"/gtf/", organism, "/")
  avail_files<- data.frame(strsplit(getURL(url, dirlistonly=TRUE),"\r\n", fixed=T))
  gtf_file<-avail_files[grepl("[0-9].gtf.gz", avail_files[,1]),1]



  source_file = paste0(url, gtf_file)
  dest_file = paste0(temp_dir, gtf_file)

  if (! file.exists(dest_file)){
    tryCatch(curl::curl_download(url = source_file, destfile = dest_file, mode = "wb"))
    return(dest_file)
  } else {
    return(dest_file)
  }
}








###create annotation for Kallisto transcript to gene id information from gtf file


# https://github.com/skurscheid/deepToolsUtils/blob/master/R/preProcessingEnsembl.R


##EXAMPLE : load_pack(stringr)

suppressPackageStartupMessages(load_pack(pacman))
suppressPackageStartupMessages(load_pack(tidyverse))
suppressPackageStartupMessages(load_pack(tidyr))
suppressPackageStartupMessages(load_pack(RCurl))
suppressPackageStartupMessages(load_pack(optparse))

#devtools::install_github("bartongroup/rats", ref="master")
suppressPackageStartupMessages(load_pack(rats))
#BiocManager::install("COMBINE-lab/wasabi")
suppressPackageStartupMessages(load_pack(GenomicFeatures))
suppressPackageStartupMessages(load_pack(rtracklayer))
suppressPackageStartupMessages(load_pack(devtools))
suppressPackageStartupMessages(load_pack(stringr))

#CheckAvailableEnsemblReleases()
#CheckAvailableEnsemblOrganisms(release = "release-98")

#pacman::p_load("tidyverse", "tidyr", rats, wasabi, GenomicFeatures, rtracklayer)

#source("https://bioconductor.org/biocLite.R")
#if (!require("ggsci")) biocLite("ggsci")


# # input species and GTF location
# option_list<-list(
#   make_option(c("--species"), help="Taxon species full name; e.g. mus_musculus"),
#   make_option(c("--version"), help="Ensembl release version e.g. 102"),
#   make_option(c("--loc"), help="Location of the reference genome files")
# )

#get arguments
# option.parser<- OptionParser(option_list=option_list)
# opt<-parse_args(option.parser)

# to test help
#  parse_args(option.parser, args=c("--help"))
# to test mouse, version 98, downloaded locally
# opt<-parse_args(option.parser, args=c("--species=mus_musculus",
#                                      "--version=98",
#                                      "--loc=D:/Box Sync/Projects/ref_genomes/mm/GRCm38/v98/"))






#' Prepare annotation for RNA-seq data analysis from Ensembl
#'
#' @param species valid Ensembl species name
#' @param ver version number
#' @param loc location of where to store annotation files
#'
#' @return
#' @export
#'
#' @examples

prepare_annotation<-function(species, ver, loc){

  options(warn=-1)
organism<-species
version<-paste0("release-", ver)
ref_dir<-loc





# Capture input errors
if (is.null(organism)){
  stop("No organism specified.\n")
}
if (is.null(version)){
  stop("Please specify the version number (e.g. 102)\n")
}
if (is.null(ref_dir)){
  stop("No reference directory specified.\n")
}

cat("Downloading ", organism, " GTF file, ", version, " to ", ref_dir, "\n\n")
annotation_file<-DownloadEnsemblGTF(version = version,
                                    temp_dir=ref_dir,
                                    organism =  organism)


ref_genome_dir<-gsub("(.*)/.*","\\1",annotation_file)

annotation_dir<-paste0(ref_genome_dir, "/annotation.v",  sub(".*-", "", version), "_", organism)
dir.create(annotation_dir)
cat("Creating output dir: ", annotation_dir)
setwd(annotation_dir)

GTF.file.local<-list.files(path = "../", pattern="*.gtf.gz")[1]
GTF.file.local<-file.path("..", GTF.file.local, fsep="/")
GTF.name<-substr(GTF.file.local,1,nchar(GTF.file.local)-3)


###RATS analysis of transcriopt level expression
###https://github.com/bartongroup/RATS

# Extract transcript ID to gene ID index from a GTF annotation.
myannot <- annot2ids(GTF.file.local)
tx2gene<-myannot
names(tx2gene)<-c("tx_id",	"gene_id")

GTF.file<-import(GTF.file.local)
Gene.Model<-as.data.frame(unique(mcols(GTF.file)[,
                                                 c( "transcript_id","gene_id", "gene_name", "gene_biotype")])) # adding biotype and genename does not increase NAs
#rownames(Gene.Model)<-GTF.file$transcript_id
Gene.Model<-Gene.Model[complete.cases(Gene.Model),]
names(Gene.Model)<-c("tx_id", "gene_id", "gene_name", "gene_biotype")
unique(Gene.Model$gene_biotype) #test all levels of biotype present

write.csv(Gene.Model, file=paste0(str_split(GTF.name, "/")[[1]][2],".Gene.Model.from.GTF.csv"), row.names = FALSE, quote = FALSE)
##to annotate fully with strand and exon positions
GTF.annotated<-GTF.file
GTF.annotated$chr <- as.character(seqnames(GTF.file))
GTF.annotated$start_pos <- start(GTF.file)
GTF.annotated$end_pos <- end(GTF.file)
GTF.annotated$chr_strand<- as.character(strand(GTF.file))
annotation<-as.data.frame(unique(mcols(GTF.annotated)[,
                                                      c( "transcript_id","gene_id", "gene_name",
                                                         "gene_biotype", "chr", "exon_number",
                                                         "start_pos", "end_pos",
                                                         "chr_strand")]))
# keep only those columns with NA in transcript_id
annotation.to.keep<-is.na(annotation$transcript_id)
annotation.by.gene <- annotation[annotation.to.keep,]
##load annotation model and process out duplicate entries due to having more than one transcript
annotation.file<-paste0(str_split(GTF.name, "/")[[1]][2],".Gene.Model.from.GTF.csv")
annotation<-read.csv(annotation.file)
annotation<-annotation[,2:4]
annotation<-unique(annotation)
#####
tx_organism<-gsub("(^)([[:alpha:]])", "\\1\\U\\2", organism, perl=TRUE)
tx_organism<-gsub("_", " ", tx_organism)
txdb<-makeTxDbFromGFF(GTF.file.local, organism=tx_organism)
txdbByGene<-transcriptsBy(txdb, "gene")
df.txdbByGene<-as.data.frame(txdbByGene)
names(df.txdbByGene)[2]<-"gene_id"
gene.ids<-Gene.Model[,c("gene_id", "gene_name")]
gene.ids<-gene.ids[!duplicated(gene.ids),]
df.annotation<-merge(df.txdbByGene, gene.ids, by="gene_id")
###tidy annotation setup
data.order<-c("group", "seqnames", "strand", "gene_id", "gene_name",
              "tx_id", "tx_name", "start", "end", "width")
##note width is not the length of the transcript!
##to calculate transcript length
##transcript length
df.txlen <- transcriptLengths(txdb, with.utr5_len=TRUE, with.utr3_len=TRUE)
df.txlen<-df.txlen[,c(2,4,5,6,7)]
df.annotation<-merge(df.annotation, df.txlen, by="tx_name")
data.order<-c(data.order, "tx_len", "nexon", "utr5_len", "utr3_len")
df.annotation<-df.annotation[,data.order]
df.annotation<-df.annotation[order(df.annotation[2]),]
##tidy R memory
rm(myannot)
rm(data.order)
rm(txdb)
rm(df.txdbByGene)
##output necessary files for annotation
write.csv(Gene.Model, file=paste0(str_split(GTF.name, "/")[[1]][2],".Gene.Model.from.GTF.csv"), row.names = FALSE, quote = FALSE)
write.csv(annotation, file=paste0(str_split(GTF.name, "/")[[1]][2], ".annotation.exons.from.GTF.csv"), row.names = FALSE, quote = FALSE)
write.csv(Gene.Model, "tx2gene.csv", quote=F, row.names = F) # save tx2gene id conversion table created from gtf for future use
write.csv(df.annotation, file=paste0(str_split(GTF.name, "/")[[1]][2],".annotation.full.from.GTF.csv"), row.names = FALSE, quote = FALSE)
save.image(file=paste0(str_split(GTF.name, "/")[[1]][2], "_annotation_data.RData"), compress = "gzip")
cat("Annotation successfully created.\n\n")
list.files(annotation_dir)
options(warn=0)
}
