% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_DEGs_Heatmaps.r
\name{create_DEGs_Heatmaps}
\alias{create_DEGs_Heatmaps}
\title{Create Heatmaps of Differentially Expressed Genes}
\usage{
create_DEGs_Heatmaps(dataFolder, geneList = "", outputFileName)
}
\arguments{
\item{dataFolder}{Location of .RData file from Differential Gene Expression analysis (absolute path)}

\item{geneList}{a list object of genes of interest, either Ensembl ID or Hugo gene names, or "" or blank}

\item{outputFileName}{string of characters for file output}
}
\value{
outputs files as PDF for important as needed
}
\description{
input can be from pairwise analysis of grouped - see keycode:extra for Differential Expression Analysis markdown document.
}
\examples{
dataFolder<-"D:/Cloud/GoogleDrive/Work UCL/Documents/Manuscripts/Shared_Pathways_MAC/data/BM-PX-WT_t35/"
outputFileName<-"testing"
create_DEGs_Heatmaps(dataFolder=dataFolder, geneList = "", outputFileName = outputFileName)

}
