#' Read custom input file and re-structure it with InterCellar format
#'
#' @param tab custom input table
#' @param separator character that separates two elements of an interaction pair
#'
#' @return preprocessed table
#' @importFrom plyr mapvalues

read.customInput <- function(tab, separator){
    # Re-format int_pair
    tab$int_pair <- sub(separator, " & ", tab$int_pair)
    
    # Check if gene columns are there, if not, generate them from int_pair
    if(!all(c("geneA", "geneB") %in% colnames(tab))){
        tab$geneA <- unlist(sapply(strsplit(
            tab$int_pair, " & "), function(x) x[1]))
        tab$geneB <- unlist(sapply(strsplit(
            tab$int_pair, " & "), function(x) x[2]))
    }
    # Add column for interaction type
    tab$int.type <- ifelse(tab$clustA == tab$clustB, "autocrine", "paracrine")
    # rename value column to score
    colnames(tab) <- plyr::mapvalues(colnames(tab), from = "value", to = "score")
    # order columns
    ord.cols <- c("int_pair", "geneA", "geneB", "typeA", "typeB", "clustA", 
                  "clustB", "score", "p_value", "int.type")
    tab <- tab[, na.omit(match(ord.cols, colnames(tab)))]
    
    return(tab)

}