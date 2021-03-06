
#' Get table for gene-verse
#'
#' @param input.data preprocessed input data
#'
#' @return gene table with unique intpairs (no connection to clusters)
#' 
#' @importFrom dplyr select %>% distinct
#' @importFrom utils read.csv
#' @importFrom tibble add_column
#' @importFrom scales rescale
#' @importFrom biomaRt getBM useMart
#' @export
#' @examples 
#' data(input.data)
#' gene_table <- getGeneTable(input.data)
getGeneTable <- function(input.data){
    if("annotation_strategy" %in% colnames(input.data)){
    gene_tab <- input.data %>%
        dplyr::select(int_pair, geneA, geneB, 
                      typeA, typeB, annotation_strategy) %>%
        distinct()
    } else {
        gene_tab <- input.data %>%
            dplyr::select(int_pair, geneA, geneB, 
                          typeA, typeB) %>%
            distinct()
    }
    if("annotation_strategy" %in% colnames(input.data)){
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
    } else {
        ensembl <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
        # query biomaRt for geneA
        bm <- biomaRt::getBM(attributes = c('hgnc_symbol', 
                                   'uniprotswissprot', 
                                   'ensembl_gene_id'), 
              filters = 'hgnc_symbol', 
              values = gene_tab$geneA, 
              mart = ensembl)
        # remove rows with empty values
        bm <- bm %>%
            filter(uniprotswissprot != "" & ensembl_gene_id != "")
        ### Adding for geneA
        bm_sub <- bm[match(gene_tab$geneA, bm$hgnc_symbol),]
        gene_tab <- add_column(gene_tab, 
                               "uniprotA" = uniprotLink(as.character(bm_sub$uniprotswissprot)), 
                               .after = "geneA")
        gene_tab <- add_column(gene_tab, 
                               "ensemblA" = ensemblLink(as.character(bm_sub$ensembl_gene_id)), 
                               .after = "uniprotA")
        # query biomaRt for geneA
        bm <- biomaRt::getBM(attributes = c('hgnc_symbol', 
                                   'uniprotswissprot', 
                                   'ensembl_gene_id'), 
                    filters = 'hgnc_symbol', 
                    values = gene_tab$geneB, 
                    mart = ensembl)
        # remove rows with empty values
        bm <- bm %>%
            filter(uniprotswissprot != "" & ensembl_gene_id != "")
        ### Adding for geneB
        bm_sub <- bm[match(gene_tab$geneB, bm$hgnc_symbol),]
        
        gene_tab <- add_column(gene_tab, "uniprotB" = uniprotLink(as.character(bm_sub$uniprotswissprot)), 
                               .after = "geneB")
        gene_tab <- add_column(gene_tab, "ensemblB" = ensemblLink(as.character(bm_sub$ensembl_gene_id)), 
                               .after = "uniprotB")
    }
    
    
    # computing uniqueness score
    n_tot_paths <- length(getClusterNames(input.data))^2
    score <- c()
    for(i in seq_len(nrow(gene_tab))){
        ip <- gene_tab$int_pair[i]
        score[i] <- 1- (nrow(filter(input.data, int_pair == ip))/n_tot_paths)
    }
    
    gene_tab <- add_column(gene_tab, 
                           "uniqueness_score" = round(
                               scales::rescale(score, to = c(0,1)), digits = 2),
                           .after = "int_pair")
    

    return(gene_tab)
}


#' Get html link to uniprot
#'
#' @param uniprot symbol
#'
#' @return html link to website
#'

uniprotLink <- function(uniprot){
    paste0('<a href="https://www.uniprot.org/uniprot/', uniprot, 
           '" target="_blank">', uniprot, '</a>')
}

#' Get html link to ensembl
#'
#' @param ensembl symbol
#'
#' @return html link to website
#'

ensemblLink <- function(ensembl){
    paste0('<a href="https://www.ensembl.org/Homo_sapiens/geneview?gene=', 
           ensembl, '" target="_blank">', ensembl, '</a>')
}


# 
#' Get number of unique ligands and receptors
#'
#' @param gene.table gene table of unique int-pairs
#' @param type either L or R
#'
#' @return number of L or R genes
#'

#' @importFrom dplyr distinct filter
getNumLR <- function(gene.table, type){
    LR.df <- data.frame(gene = c(gene.table$geneA, gene.table$geneB),
                        type = c(gene.table$typeA, gene.table$typeB))
    LR.df <- distinct(LR.df)
    if(type == "L"){
        return(nrow(filter(LR.df, type == "L")))
    } else if(type == "R"){
        return(nrow(filter(LR.df, type == "R")))
    }
}






#' Functions to plot DotPlots 
#'
#' @param selected_tab table of selected rows from gene tableeeeeee
#' @param clust.order how to order clusters
#' @param low_color of dotplot
#' @param high_color of dotplot
#'
#' @return list with modified selected data and ggplot2 dotplot
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

#' Plot dotplot containing only unique int-pair/cluster pairs with many conditions
#'
#' @param data_dotplot table with selected int_pairs for multiple conditions
#'
#' @return ggplot object

getUniqueDotplot <- function(data_dotplot){
    distinct_pairs_clust <- data_dotplot %>%
        group_by(int_pair, cluster_pair) %>%
        mutate(n = n()) %>%
        filter(n == 1)
    
    
    g <- ggplot(distinct_pairs_clust, aes(x = int_pair, y = cluster_pair)) +
        geom_point(aes(color = condition, size=4)) + 
        theme_minimal() +
        scale_color_manual(values = c("#00BA38","#619CFF","#F8766D")) + 
        facet_grid(groups_x ~ ., scales = "free", space = "free") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
              text = element_text(size=20),
              strip.text.y = element_blank(),
              strip.text.x = element_text(angle = 0)) +
        guides(size = FALSE) + 
        labs(x = "Int-pairs", y = "Cluster-pairs")
    return(g)
}