% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_space_Revigo.r
\name{plot_space_Revigo}
\alias{plot_space_Revigo}
\title{Plots 2D semantic space of ReviGO results for gene ontology analysis
data format:
(df format TermID,	Name,	Value,	LogSize,	Frequency,	Uniqueness,	Dispensability,	PC_0,	PC_1,	Representative)}
\usage{
plot_space_Revigo(
  revigoData,
  dispens,
  reduce = TRUE,
  output,
  scale = "as-is",
  ...
)
}
\arguments{
\item{revigoData}{\itemize{
\item dataframe of the output table from Revigo analysis
}}

\item{dispens}{\itemize{
\item numeric filter threshold for dispensibility
}}

\item{reduce}{\itemize{
\item reduces the points plotted in line with dispensibility
}}

\item{output}{\itemize{
\item filename with full path of output plot (.svg)
}}

\item{scale}{\itemize{
\item scaling of Value; as-is (default) vs 'log'
}}

\item{...}{}
}
\description{
Plots 2D semantic space of ReviGO results for gene ontology analysis
data format:
(df format TermID,	Name,	Value,	LogSize,	Frequency,	Uniqueness,	Dispensability,	PC_0,	PC_1,	Representative)
}
