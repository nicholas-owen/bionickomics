#' Get genes for a specific KEGG pathway
#'
#' returns geneList a list of all human and zebrafish gene IDs related to pathway KEGG with keyword used.
#' To get a complete list of pathways to subsequently search, use the keyword 'ALL' .
#'
#' @param pathwayquery search term keyword (human)



get_KEGG_genes <- function(pathwayquery) {
  load_pack("KEGGREST")
  load_pack("EnrichmentBrowser")
  load_pack("gage")
  load_pack("biomaRt")

  hsa <-
    getGenesets(
      org = "hsa",
      db = "kegg",
      cache = TRUE,
      return.type = "list"
    )

  if (pathwayquery == "ALL"){
    message("Returning list object of all KEGG pathways.")
    return(names(hsa))
  }

  keggPathway <-
    hsa[grep(pathwayquery, names(hsa),  ignore.case = TRUE)]

  message("Obtaining KEGG pathway information using search term: ",
          pathwayquery)

  message("Pathway found - ", names(hsa)[grep(pathwayquery, names(hsa),  ignore.case = TRUE)])
  ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  keggIDs <-
    getBM(
      attributes = c('ensembl_gene_id', 'hgnc_symbol', "entrezgene_id") ,
      filters = 'entrezgene_id',
      values = keggPathway,
      mart = ensembl
    )
  message("Grabbing zebrafish orthologues as a bonus!")
  zfIDs <-
    getBM(
      attributes = c(
        'ensembl_gene_id',
        "drerio_homolog_ensembl_gene",
        "drerio_homolog_associated_gene_name"
      ) ,
      filters = 'ensembl_gene_id',
      values = keggIDs$ensembl_gene_id,
      mart = ensembl
    )


  keggIDs <- merge (
    x = keggIDs,
    y = zfIDs,
    by.x = "ensembl_gene_id",
    by.y = "ensembl_gene_id"
  )

  keggIDs<-(keggIDs[order(keggIDs$ensembl_gene_id, decreasing = FALSE), ]   )

  message("Here you go.")
  #keggPathway
  #write.table(keggIDs, file=paste0(names(keggPathway), ".csv"), quote = FALSE, sep=",", row.names = FALSE)

  return(keggIDs)
}

