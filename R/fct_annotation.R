

#' Connection to Ensembl via biomaRt to get GO terms
#'
#' @param input_select_ensembl chosen version of Ensembl
#' @param input.data filtered input data
#'
#' @return dataframe with GO annotation 
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom dplyr mutate_all

getGObiomaRt <- function(input_select_ensembl, 
                         input.data){
    
    genesA <- unlist(strsplit(input.data$geneA, ","))
    genesB <- unlist(strsplit(input.data$geneB, ","))
    genes <- unique(c(genesA,genesB))
    
    ensembl <- biomaRt::useEnsembl(biomart = 'genes', 
                          dataset = 'hsapiens_gene_ensembl',
                          version = input_select_ensembl)
    GO.biomart <-  biomaRt::getBM(mart = ensembl,
                                  attributes = c("hgnc_symbol",
                                                 "go_id", 
                                                 "name_1006", 
                                                 "namespace_1003", 
                                                 "go_linkage_type"), 
                                  filters = "hgnc_symbol",
                                  values = genes)
    GO.biomart <- GO.biomart %>% 
        dplyr::mutate_all(~ifelse(. %in% c("N/A", "null", ""), NA, .)) %>% 
        na.omit()
    
    # Rename columns 
    colnames(GO.biomart) <- plyr::mapvalues(colnames(GO.biomart), 
                                            from = c("hgnc_symbol",
                                                     "name_1006", 
                                                     "namespace_1003") , 
                                            to = c("gene_symbol",  
                                                   "go_term", "domain"))
    return(GO.biomart)
}

#' Get pathway databases from graphite package
#'
#' @param species for now h_sapiens
#' @param database name of the database to download
#'
#' @return dataframe containing the downloaded database to pre-process
#' @importFrom graphite pathways convertIdentifiers
getGraphiteDB <- function(species, database){
    db <- graphite::pathways(species = species, database = database)
    db <- graphite::convertIdentifiers(db, to = "SYMBOL")
    
    return(db)
}