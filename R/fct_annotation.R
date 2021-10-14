

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
    out <- tryCatch({
    GO.biomart <-  biomaRt::getBM(mart = ensembl,
                                  attributes = c("hgnc_symbol",
                                                 "go_id", 
                                                 "name_1006", 
                                                 "namespace_1003", 
                                                 "go_linkage_type"), 
                                  filters = "hgnc_symbol",
                                  values = genes,
                                  useCache = TRUE)
    },
    error = function(cond){
        message(cond)
    },
    warning = function(cond){
        message(cond)
    })
    
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

