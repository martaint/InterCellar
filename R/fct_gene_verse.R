# 
#' Get table for gene-verse
#'
#' @param input.data 
#'
#' @return
#' @importFrom dplyr select %>% distinct
#' @importFrom utils read.csv
#' @importFrom tibble add_column
#' @examples
getGeneTable <- function(input.data){
    gene_tab <- input.data %>%
        dplyr::select(int_pair, geneA, geneB, 
                      typeA, typeB, annotation_strategy) %>%
        distinct()
    # get protein info from cpdb_v2 file
    gene_input <- read.csv(app_sys("app", "extdata", "cpdb_gene_input.csv"))
    
    ### Adding for geneA
    gene_input_sub <- gene_input[match(gene_tab$geneA, gene_input$hgnc_symbol), 
                                 c("uniprot", "ensembl")]
    gene_tab <- add_column(gene_tab, 
                           "uniprotA" = uniprotLink(as.character(gene_input_sub$uniprot)), 
                           .after = "geneA")
    gene_tab <- add_column(gene_tab, 
                           "ensemblA" = ensemblLink(as.character(gene_input_sub$ensembl)), 
                           .after = "uniprotA")
    # considering complexes
    ind_compl <- grep(",", gene_tab$geneA)
    for(i in ind_compl){
        genes <- unlist(strsplit(gene_tab$geneA[i], ","))
        prots <- uniprotLink(as.character(
            gene_input[match(genes, gene_input$hgnc_symbol), "uniprot"]))
        ensembls <- ensemblLink(as.character(
            gene_input[match(genes, gene_input$hgnc_symbol), "ensembl"]))
        gene_tab$uniprotA[i] <- paste(prots, collapse = ",")
        gene_tab$ensemblA[i] <- paste(ensembls, collapse = ",")
    }
    
    ### Adding for geneB
    gene_input_sub <- gene_input[match(gene_tab$geneB, gene_input$hgnc_symbol), 
                                 c("uniprot", "ensembl")]
    gene_tab <- add_column(gene_tab, "uniprotB" = uniprotLink(as.character(gene_input_sub$uniprot)), 
                           .after = "geneB")
    gene_tab <- add_column(gene_tab, "ensemblB" = ensemblLink(as.character(gene_input_sub$ensembl)), 
                           .after = "uniprotB")
    # considering complexes
    ind_compl <- grep(",", gene_tab$geneB)
    for(i in ind_compl){
        genes <- unlist(strsplit(gene_tab$geneB[i], ","))
        prots <- uniprotLink(as.character(
            gene_input[match(genes, gene_input$hgnc_symbol), "uniprot"]))
        ensembls <- ensemblLink(as.character(
            gene_input[match(genes, gene_input$hgnc_symbol), "ensembl"]))
        gene_tab$uniprotB[i] <- paste(prots, collapse = ",")
        gene_tab$ensemblB[i] <- paste(ensembls, collapse = ",")
    }
    
    
    return(gene_tab)
}


#' Get html link to uniprot
#'
#' @param uniprot 
#'
#' @return
#'
#' @examples
uniprotLink <- function(uniprot){
    paste0('<a href="https://www.uniprot.org/uniprot/', uniprot, 
           '" target="_blank">', uniprot, '</a>')
}

#' Get html link to ensembl
#'
#' @param ensembl 
#'
#' @return
#'
#' @examples
ensemblLink <- function(ensembl){
    paste0('<a href="https://www.ensembl.org/Homo_sapiens/geneview?gene=', 
           ensembl, '" target="_blank">', ensembl, '</a>')
}


# 
#' Get number of unique ligands and receptors
#'
#' @param input.data 
#' @param type 
#'
#' @return
#'
#' @examples
#' @importFrom dplyr distinct filter
getNumLR <- function(input.data, type){
    LR.df <- rbind(input.data[, c("geneA", "typeA")], 
                   input.data[, c("geneB", "typeB")],
                   use.names = FALSE)
    LR.df <- distinct(LR.df)
    if(type == "L"){
        return(nrow(filter(LR.df, typeA == "L")))
    } else if(type == "R"){
        return(nrow(filter(LR.df, typeA == "R")))
    }
}






#' Functions to plot DotPlots 
#'
#' @param selected_tab 
#' @param clust.order 
#' @param low_color 
#' @param high_color 
#'
#' @return
#' @importFrom tidyr unite
#' @import ggplot2 

getDotPlot_selInt <- function(selected_tab, clust.order, 
                              low_color = "aquamarine", 
                              high_color = "#131780"){
    selected_tab$groups_x <- factor(selected_tab$clustA, levels = clust.order)
    selected_tab <- selected_tab %>%
        unite(col="cluster_pair", clustA:clustB, sep = "::", remove = TRUE)
    p <- ggplot(selected_tab, aes(x = int_pair, y = cluster_pair)) +
        geom_point(aes(color = score, size=3.5)) + 
        theme_minimal() +
        scale_color_gradient(low = low_color, high = high_color) + 
        facet_grid(groups_x ~ ., scales = "free", space = "free") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
              text = element_text(size=20),
              strip.text = element_blank()) +
        guides(size = FALSE) + 
        labs(x = "Int-pairs", y = "Cluster-pairs")
    
    return(list(data_dot = selected_tab, p = p))
}
