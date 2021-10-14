#' Read output from CellPhoneDB v2. 
#' @description  Output is a folder containing 4 .txt files 
#' - deconvoluted.txt: containing list of single genes and their mean 
#' expression in each cluster (not considered);
#' - means.txt: containing list of interacting pairs with info regarding L/R, 
#' annotation strategy and mean value of all pairs over cluster couples.
#' - pvalues.txt: same as means, but containing pvalue of each pair, for each 
#' cluster couple.
#' - significant_means.txt: only means of those pairs that have pvalue < 0.05. 
#' Has one more column:rank. 
#' If the statistical analysis is not run, the folder would contain only 
#' deconvoluted and means
#'
#' @param folder folder containing output
#'
#' @return input.data which is the pre-processed object with annotated L-R pairs
#' @importFrom data.table data.table
#' @importFrom tidyr gather
#' @importFrom utils read.csv read.table
#' @importFrom dplyr %>%


read.CPDBv2 <- function(folder){
    files <- list.files(folder)
    if("significant_means.txt" %in% files){
        # CellPhoneDB run with statistical analysis
        s_mean <- read.table(
            file = file.path(folder, "significant_means.txt"), 
            header = TRUE, sep = "\t", 
            check.names = FALSE, stringsAsFactors = FALSE)
        pvalues <- read.table(
            file = file.path(folder, "pvalues.txt"), 
            header = TRUE, sep = "\t", 
            check.names = FALSE, stringsAsFactors = FALSE)
        # remove erroneously duplicated rows
        s_mean <- s_mean[!duplicated(s_mean$interacting_pair), ]
        pvalues <- pvalues[!duplicated(pvalues$interacting_pair), ]
        
        # get metadata for interaction pairs
        metadata <- s_mean[, which(colnames(s_mean) == "id_cp_interaction"):
                               which(colnames(s_mean) == "rank")]
        # get table of interactions vs cluster couples
        table.means <- cbind(interacting_pair = s_mean[, "interacting_pair"], 
                             s_mean[, (which(colnames(s_mean) == "rank")+1)
                                    :ncol(s_mean)])
        table.pvalues <- cbind(interacting_pair = pvalues[, "interacting_pair"],
                               pvalues[, (which(colnames(pvalues) == 
                                                    "is_integrin")+1):
                                           ncol(pvalues)])
        
        # convert tables to long format with tidyr
        table.means.long <- tidyr::gather(table.means, cluster_pair, mean_value,
                                          -interacting_pair, na.rm = TRUE)
        table.pvalues.long <- tidyr::gather(table.pvalues, cluster_pair, `p_value`,
                                            -interacting_pair, na.rm = TRUE)
        # table pvalues contains pvalues for all pairs, not only the significant ones!
        # input data built using only significant means
        table.means.long <- merge(table.means.long, metadata, 
                                  by="interacting_pair", all.x=TRUE)
        table.long <- merge(table.means.long, table.pvalues.long, 
                            by=c("interacting_pair", "cluster_pair"), 
                            all.x=TRUE)
        # Check if CPDB table contains ensembl ids instead of gene names
        gene_input <- read.csv(app_sys("app", "extdata", "cpdb_gene_input.csv"))
        if(any(grepl("ENSG00", table.long$gene_a)) | 
           any(grepl("ENSG00", table.long$gene_b))){
            table.long$gene_a[grep("ENSG00", table.long$gene_a)] <- as.character(
                gene_input[match(grep("ENSG00", table.long$gene_a, value = TRUE), 
                                 gene_input$ensembl), "hgnc_symbol"])
            table.long$gene_b[grep("ENSG00", table.long$gene_b)] <- as.character(
                gene_input[match(grep("ENSG00", table.long$gene_b, value = TRUE), 
                                 gene_input$ensembl), "hgnc_symbol"])
            
        }
        
        int_pair_a <- ifelse(grepl("complex", table.long$partner_a), 
                             sub("complex:", "", table.long$partner_a), 
                             table.long$gene_a)
        int_pair_b <- ifelse(grepl("complex", table.long$partner_b), 
                             sub("complex:", "", table.long$partner_b), 
                             table.long$gene_b)
        
        input.data <- data.table(int_pair = paste(int_pair_a, int_pair_b, 
                                                  sep = " & "),
                                 geneA=table.long$gene_a, 
                                 geneB=table.long$gene_b, 
                                 typeA=ifelse(table.long$receptor_a == "True", 
                                              "R", "L"), 
                                 typeB=ifelse(table.long$receptor_b == "True", 
                                              "R", "L"), 
                                 clustA=unlist(sapply(strsplit(
                                     table.long$cluster_pair, "\\|"), 
                                     function(x) x[1])), 
                                 clustB=unlist(sapply(strsplit(
                                     table.long$cluster_pair, "\\|"), 
                                     function(x) x[2])), 
                                 score=table.long$mean_value, 
                                 `p_value`=table.long$p_value,
                                 annotation_strategy=table.long$annotation_strategy,
                                 stringsAsFactors = FALSE)
        
        
        input.data$int.type <- ifelse(input.data$clustA == input.data$clustB, 
                                      "autocrine", "paracrine")
        
        
        
    } else{
        # CellPhoneDB run without statistical analysis: considering means
        means <- read.table(file = file.path(folder, "means.txt"), 
                            header = TRUE, sep = "\t", check.names = FALSE, 
                            stringsAsFactors = FALSE)
        # remove erroneously duplicated rows
        means <- means[!duplicated(means$interacting_pair), ]
        
        # get metadata for interaction pairs
        metadata <- means[, which(colnames(means) == "id_cp_interaction"):
                              which(colnames(means) == "is_integrin")]
        # get table of interactions vs cluster couples
        table.means <- cbind(interacting_pair = means[, "interacting_pair"], 
                             means[, (which(colnames(means) == "is_integrin")+1)
                                   :ncol(means)])
        
        # convert tables to long format with tidyr
        table.means.long <- gather(table.means, cluster_pair, mean_value, 
                                   -interacting_pair, na.rm = TRUE)
        table.long <- merge(table.means.long, metadata, by="interacting_pair", 
                            all.x=TRUE)
        # remove interactions that have mean = 0
        table.long <- table.long[table.long$mean_value > 0,]
        
        # Check if CPDB table contains ensembl ids instead of gene names
        gene_input <- read.csv(app_sys("app", "extdata", "cpdb_gene_input.csv"))
        if(any(grepl("ENSG00", table.long$gene_a)) | 
           any(grepl("ENSG00", table.long$gene_b))){
            table.long$gene_a[grep("ENSG00", table.long$gene_a)] <- as.character(
                gene_input[match(grep("ENSG00", table.long$gene_a, value = TRUE), 
                                 gene_input$ensembl), "hgnc_symbol"])
            table.long$gene_b[grep("ENSG00", table.long$gene_b)] <- as.character(
                gene_input[match(grep("ENSG00", table.long$gene_b, value = TRUE), 
                                 gene_input$ensembl), "hgnc_symbol"])
            
        }
        
        int_pair_a <- ifelse(grepl("complex", table.long$partner_a), 
                             sub("complex:", "", table.long$partner_a), 
                             table.long$gene_a)
        int_pair_b <- ifelse(grepl("complex", table.long$partner_b), 
                             sub("complex:", "", table.long$partner_b), 
                             table.long$gene_b)
        
        input.data <- data.table(int_pair = paste(int_pair_a, int_pair_b, 
                                                  sep = " & "),
                                 geneA=table.long$gene_a, 
                                 geneB=table.long$gene_b, 
                                 typeA=ifelse(table.long$receptor_a == "True",
                                              "R", "L"), 
                                 typeB=ifelse(table.long$receptor_b == "True",
                                              "R", "L"), 
                                 clustA=unlist(sapply(strsplit(
                                     table.long$cluster_pair, "\\|"), 
                                     function(x) x[1])), 
                                 clustB=unlist(sapply(strsplit(
                                     table.long$cluster_pair, "\\|"), 
                                     function(x) x[2])), 
                                 score=table.long$mean_value, 
                                 annotation_strategy=table.long$annotation_strategy,
                                 stringsAsFactors = FALSE)
        
        
        input.data$int.type <- ifelse(input.data$clustA == input.data$clustB, 
                                      "autocrine", "paracrine")
    }
    
    # get gene info on complexes 
    complex_input <- read.csv(app_sys("app", "extdata", "cpdb_complex_input.csv"))
    
    
    complex_input$hgnc_symbol_1 <- as.character(gene_input[match(
        complex_input$uniprot_1, gene_input$uniprot), "hgnc_symbol"])
    complex_input$hgnc_symbol_2 <- as.character(gene_input[match(
        complex_input$uniprot_2, gene_input$uniprot), "hgnc_symbol"])
    complex_input$hgnc_symbol_3 <- as.character(gene_input[match(
        complex_input$uniprot_3, gene_input$uniprot), "hgnc_symbol"])
    complex_input$hgnc_symbol_4 <- as.character(gene_input[match(
        complex_input$uniprot_4, gene_input$uniprot), "hgnc_symbol"])
    
    
    ind.compl.a <- which(input.data$geneA == "")
    if(length(ind.compl.a)>0){
        compl.a <- unlist(sapply(strsplit(input.data$int_pair[ind.compl.a], "&"),
                                 function(x) trimws(x[1])))
        genes.compl.a <- complex_input[match(compl.a, complex_input$complex_name),
                                       c("hgnc_symbol_1", "hgnc_symbol_2", 
                                         "hgnc_symbol_3", "hgnc_symbol_4")]
        input.data[ind.compl.a, "geneA"] <- apply(genes.compl.a, 1, function(x) 
            paste0(x[!is.na(x)], collapse = ","))
    }
    
    ind.compl.b <- which(input.data$geneB == "")
    if(length(ind.compl.b)>0){
        compl.b <- unlist(sapply(strsplit(input.data$int_pair[ind.compl.b], "&"), 
                                 function(x) trimws(x[2])))
        genes.compl.b <- complex_input[match(compl.b, complex_input$complex_name), 
                                       c("hgnc_symbol_1", "hgnc_symbol_2", 
                                         "hgnc_symbol_3", "hgnc_symbol_4")]
        input.data[ind.compl.b, "geneB"] <- apply(genes.compl.b, 1, function(x) 
            paste0(x[!is.na(x)], collapse = ","))
    }
    
    ## Re-annotate specific genes to R / L
    input.data <- checkLL_RR(input.data)
    ## Update input.data with ordered L-R interactions
    input.data <- updateInputLR(input.data)
    
    return(input.data)
}

#' Manually change the annotation of L-L and R-R pairs
#'
#' @param input.data preprocessed table
#'
#' @return input.data
#' @export
#' @examples
#' data(input.data)
#' checked.input.data <- checkLL_RR(input.data)
#' 
checkLL_RR <- function(input.data){

    
    # to re-define as L-R
    
    # integrins
    integrins <- unique(input.data$int_pair[intersect(grep("complex", input.data$int_pair),
                                                      grep("ITG", input.data$geneB))])
    
    l_r <- c("ALOX5 & ALOX5AP", "CD6 & ALCAM", integrins)
    input.data[input.data$int_pair %in% l_r, "typeB"] <- "R"
    
    
    ## CADM -> cell adhesion molecules are transmembrane -> R
    input.data[grep("CADM", input.data$geneA), "typeA"] <- "R"
    input.data[grep("CADM", input.data$geneB), "typeB"] <- "R"
    input.data[grep("CEACAM", input.data$geneA), "typeA"] <- "R"
    input.data[grep("CEACAM", input.data$geneB), "typeB"] <- "R"
    input.data[grep("ESAM", input.data$geneA), "typeA"] <- "R"
    input.data[grep("ESAM", input.data$geneB), "typeB"] <- "R"
    input.data[grep("PTPRZ1", input.data$geneA), "typeA"] <- "R"
    input.data[grep("PTPRZ1", input.data$geneB), "typeB"] <- "R"
    input.data[grep("TNFRSF6B", input.data$geneA), "typeA"] <- "R"
    input.data[grep("TNFRSF6B", input.data$geneB), "typeB"] <- "R"
    input.data[grep("TNFRSF11B", input.data$geneA), "typeA"] <- "R"
    input.data[grep("TNFRSF11B", input.data$geneB), "typeB"] <- "R"
    input.data[grep("CXCR", input.data$geneA), "typeA"] <- "R"
    input.data[grep("CXCR", input.data$geneB), "typeB"] <- "R"
    input.data[grep("CXCL", input.data$geneA), "typeA"] <- "L"
    input.data[grep("CXCL", input.data$geneB), "typeB"] <- "L"
    
    
    
    
    
    return(input.data)
}




#' Read output from SingleCellSignalR
#' @description SCSR description: the output folder is a collection of txt files, 
#' one for each 
#' clusters pair considered. The "paracrine" option looks for ligands expressed 
#' in cluster A and their associated receptors according to LRdb that are
#'  expressed in any other cluster but A. These interactions are labelled 
#'  "paracrine". The interactions that involve a ligand and a receptor, both 
#'  differentially expressed in their respective cell clusters according to 
#'  the **edgeR** analysis performed by the **cluster_analysis()** function, 
#'  are labelled "specific". The "autocrine" option searches for ligands 
#'  expressed in cell cluster A and their associated receptors also expressed 
#'  in A. These interactions are labelled "autocrine". Additionally, it searches
#'  for those associated receptors in the other cell clusters (not A) to cover 
#'  the part of the signaling that is "autocrine" and "paracrine" 
#'  simultaneously. These interactions are labelled "autocrine/paracrine".
#'  This file is a 4-column table: ligands, receptors, interaction types 
#'  ("paracrine", "autocrine", "autocrine/paracrine" and "specific"),
#'  and the associated LRscore. 
#'  InterCellar: rename autocrine|paracrine to paracrine
#'
#' @param folder containing output from SingleCellSignalR, named cell-signaling
#'
#' @return input.data: preprocessed object with annotated L-R pairs
#' @importFrom utils read.csv read.table
#' 
read.SCsignalR <- function(folder){
    files <- list.files(folder)
    input.data <- data.table(int_pair=character(), geneA=character(), 
                             geneB=character(), typeA=character(), 
                             typeB=character(), clustA=character(), 
                             clustB=character(), int.type=character(), 
                             score=double(), scSignalR_specific=character(), 
                             stringsAsFactors = FALSE)
    
    
    for(f in files){
        # read table
        table.tmp <- read.table(file = file.path(folder, f), header = TRUE, 
                                check.names = FALSE, stringsAsFactors = FALSE)
        input.tmp <- data.table(int_pair = paste(table.tmp[,1], table.tmp[,2], 
                                                 sep = " & "), 
                                geneA = table.tmp[,1],
                                geneB = table.tmp[,2],
                                typeA = "L",
                                typeB = "R",
                                int.type = table.tmp$interaction.type,
                                score = round(table.tmp$LRscore, 3),
                                scSignalR_specific = NA)
        
        input.tmp$clustA <- colnames(table.tmp)[1]
        input.tmp$clustB <- sub(pattern = ".1", replacement = "", 
                                colnames(table.tmp)[2])
        # rename autocrine|paracrine interactions to paracrine
        input.tmp$int.type[grepl("^autocrine\\|paracrine$", input.tmp$int.type)] <-
            "paracrine"
        # fix column scSignalR_specific
        ind.specific <- grep("specific", input.tmp$int.type)
        input.tmp$scSignalR_specific[ind.specific] <-
            "specific"
        # rename specific interaction to either auto or paracrine
        for(ind in ind.specific){
            input.tmp$int.type[ind] <- ifelse(input.tmp$clustA[ind] == 
                                                  input.tmp$clustB[ind], 
                                              "autocrine",
                                              "paracrine")
        }
        
        
        input.data <- merge(input.data, input.tmp, all = TRUE)
    }
    ## Update input.data with ordered L-R interactions
    input.data <- updateInputLR(input.data)
    
    return(input.data)
}

#' Read dataframe of cell-cell communication from CellChat (ligand/receptor)
#'
#' @param file_tab dataframe from cellchat
#' @importFrom stats complete.cases
#' @return input.data formatted for InterCellar

read.cellchat <- function(file_tab){
    # remove empty rows from dataframe
    file_tab <- file_tab[stats::complete.cases(file_tab),]
    input.data <- data.table(int_pair=character(), geneA=character(), 
                             geneB=character(), typeA=character(), 
                             typeB=character(), clustA=character(), 
                             clustB=character(), int.type=character(), 
                             score=double(), p_value=double(),
                             pathway_cellchat=character(), 
                             annotation_cellchat=character(),
                             evidence_cellchat=character(),
                             stringsAsFactors = FALSE)
    
    input.data <- data.table(int_pair = gsub(" - ", " & ", file_tab$interaction_name_2),
                             geneA = character(),
                             geneB = character(),
                             typeA = "L", typeB = "R",
                             clustA = file_tab$source,
                             clustB = file_tab$target,
                             int.type = ifelse(file_tab$source == file_tab$target, 
                                               "autocrine", "paracrine"),
                             score = round(file_tab$prob, 3),
                             p_value = file_tab$pval,
                             pathway_cellchat = file_tab$pathway_name,
                             annotation_cellchat = file_tab$annotation,
                             evidence_cellchat = file_tab$evidence)
    
    input.data$geneA <- unlist(sapply(strsplit(input.data$int_pair, "&"), 
                          function(x) trimws(x[1])))
    input.data$geneA <- gsub("\\+", ",", input.data$geneA)
    input.data$geneA <- gsub("\\(|\\)", "", input.data$geneA)
    input.data$geneB <- unlist(sapply(strsplit(input.data$int_pair, "&"), 
                                      function(x) trimws(x[2])))
    input.data$geneB <- gsub("\\+", ",", input.data$geneB)
    input.data$geneB <- gsub("\\(|\\)", "", input.data$geneB)

    
    return(input.data)
}

#' Read ICELLNET dataframe
#'
#' @param tab dataframe with int-pairs in "X" column, other columns as cell types
#' @param input_icellnet_CC central cell name
#' @param input_icellnet_dir direction of interaction either out or in
#'
#' @return pre-processed input data
#' @importFrom tidyr pivot_longer
#' @importFrom data.table data.table

read.icellnet <- function(tab, input_icellnet_CC, input_icellnet_dir){
    
    if(input_icellnet_dir == "out"){ # L on CC, R on other cells
        data.long <- tab %>%
            pivot_longer(!X,
                         names_to = "clustB",
                         values_to = "score"
            ) %>%
            mutate(clustA = input_icellnet_CC) %>%
            filter(!(is.na(score) | score == 0 ))
    } else {
        data.long <- tab %>%
            pivot_longer(!X,
                         names_to = "clustA",
                         values_to = "score"
            ) %>%
            mutate(clustB = input_icellnet_CC) %>%
            filter(!(is.na(score) | score == 0 ))
    }
    
    
    input.data <- data.table(int_pair = gsub(pattern = " / ", " & ", data.long$X),
                             geneA = unlist(sapply(strsplit(
                                 data.long$X, " / "), 
                                 function(x) x[1])),
                             geneB = unlist(sapply(strsplit(
                                 data.long$X, " / "), 
                                 function(x) x[2])),
                             typeA = "L", 
                             typeB = "R", 
                             clustA= data.long$clustA, 
                             clustB= data.long$clustB, 
                             score= round(data.long$score, 3),
                             stringsAsFactors = FALSE)
    
    input.data$geneA <- gsub(" \\+ ", ",", input.data$geneA)                                 
    input.data$geneB <- gsub(" \\+ ", ",", input.data$geneB)
    
    input.data$int.type <- ifelse(input.data$clustA == input.data$clustB, 
                                    "autocrine", "paracrine")
    
    return(input.data)
}