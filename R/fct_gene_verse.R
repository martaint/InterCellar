
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
    } else if("pathway_cellchat" %in% colnames(input.data)){
        gene_tab <- input.data %>%
            dplyr::select(int_pair, geneA, geneB, 
                          typeA, typeB, pathway_cellchat,
                          annotation_cellchat, evidence_cellchat) %>%
            distinct()
    } else {
        gene_tab <- input.data %>%
            dplyr::select(int_pair, geneA, geneB, 
                          typeA, typeB) %>%
            distinct()
    }
    if("annotation_strategy" %in% colnames(input.data)){ # CPDB
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
    } else if("pathway_cellchat" %in% colnames(input.data)){ # cellchat
        ### prepare db with all possible genes to query biomart
        all_genes <- character()
        # separating simple genes from complexes
        # simple A
        ind_compl <- grepl(",", gene_tab$geneA)
        all_genes <- c(all_genes, gene_tab$geneA[!ind_compl])
        # simple B
        ind_compl <- grepl(",", gene_tab$geneB)
        all_genes <- c(all_genes, gene_tab$geneB[!ind_compl])
        # complex A
        ind_compl <- grep(",", gene_tab$geneA)
        for(i in ind_compl){
            genes <- unlist(strsplit(gene_tab$geneA[i], ","))
            all_genes <- c(all_genes, genes)
        }
        # complex B
        ind_compl <- grep(",", gene_tab$geneB)
        for(i in ind_compl){
            genes <- unlist(strsplit(gene_tab$geneB[i], ","))
            all_genes <- c(all_genes, genes)
        }
        
        # unique all genes
        all_genes <- unique(all_genes)
        
        ensembl <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
        # query biomaRt for all genes
        bm <- biomaRt::getBM(attributes = c('hgnc_symbol', 
                                            'uniprotswissprot', 
                                            'ensembl_gene_id'), 
                             filters = 'hgnc_symbol', 
                             values = all_genes, 
                             mart = ensembl,
                             useCache = TRUE)
        # remove rows with empty values
        bm <- bm %>%
            filter(uniprotswissprot != "" & ensembl_gene_id != "") %>%
            group_by(hgnc_symbol)
        
        ### Adding for geneA
        bm_sub <- bm[match(gene_tab$geneA, bm$hgnc_symbol),]
        gene_tab <- add_column(gene_tab, 
                               "uniprotA" = uniprotLink(as.character(bm_sub$uniprotswissprot)), 
                               .after = "geneA")
        gene_tab <- add_column(gene_tab, 
                               "ensemblA" = ensemblLink(as.character(bm_sub$ensembl_gene_id)), 
                               .after = "uniprotA")
        # considering complexes for geneA
        ind_compl <- grep(",", gene_tab$geneA)
        for(i in ind_compl){
            genes <- unlist(strsplit(gene_tab$geneA[i], ","))
            prots <- uniprotLink(as.character(unlist(
                bm[match(genes, bm$hgnc_symbol), "uniprotswissprot"])))
            ensembls <- ensemblLink(as.character(unlist(
                bm[match(genes, bm$hgnc_symbol), "ensembl_gene_id"])))
            gene_tab$uniprotA[i] <- paste(prots, collapse = ",")
            gene_tab$ensemblA[i] <- paste(ensembls, collapse = ",")
        }
        
        ### Adding for geneB
        bm_sub <- bm[match(gene_tab$geneB, bm$hgnc_symbol),]
        
        gene_tab <- add_column(gene_tab, "uniprotB" = uniprotLink(as.character(bm_sub$uniprotswissprot)), 
                               .after = "geneB")
        gene_tab <- add_column(gene_tab, "ensemblB" = ensemblLink(as.character(bm_sub$ensembl_gene_id)), 
                               .after = "uniprotB")
        # considering complexes for geneB
        ind_compl <- grep(",", gene_tab$geneB)
        for(i in ind_compl){
            genes <- unlist(strsplit(gene_tab$geneB[i], ","))
            prots <- uniprotLink(as.character(unlist(
                bm[match(genes, bm$hgnc_symbol), "uniprotswissprot"])))
            ensembls <- ensemblLink(as.character(unlist(
                bm[match(genes, bm$hgnc_symbol), "ensembl_gene_id"])))
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
        # query biomaRt for geneB
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
#' @param selected_tab selected rows of filt.data by selection from gene table
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
        guides(size = "none") + 
        labs(x = "Int-pairs", y = "Cluster-pairs")
    
    return(list(data_dot = selected_tab, p = p))
}

