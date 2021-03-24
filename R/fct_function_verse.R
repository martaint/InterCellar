
#' Perform GO annotation of input data
#'
#' @param input_select_ensembl ensembl version selected by user
#' @param input_go_evidence_exclude evidence codes to exclude by user
#' @param input_go_sources_checkbox GO sources to use by user
#' @param input.data preprocessed input data
#'
#' @return GO_annotation
#' @importFrom dplyr distinct transmute filter select group_by
#' @importFrom data.table data.table
annotateGO <- function(input_select_ensembl, 
                       input_go_evidence_exclude, 
                       input_go_sources_checkbox,
                       input.data){
    # Get annotation from Ensembl (biomaRt)
    GO.biomart <- getGObiomaRt(input_select_ensembl, input.data)
    
    # Filter evidence code
    if(!is.null(input_go_evidence_exclude)){
        GO.biomart <- GO.biomart %>%
            filter(!(go_linkage_type %in% input_go_evidence_exclude))
    }
    # Filter ontologies
    GO.biomart <- GO.biomart %>%
        filter(domain %in% input_go_sources_checkbox)
    
    # get unique pairs of interactions and associated genes
    unique.interactions <- distinct(input.data, int_pair, .keep_all = TRUE)
    # separate simple interactions from complexes
    complex.index <- grepl(",", unique.interactions$geneA) | 
        grepl(",", unique.interactions$geneB)
    # only done if there are complex interactions
    if(length(complex.index) > 0){
        simple.interactions <- unique.interactions[!complex.index, 
                                                   c("int_pair", "geneA", "geneB")]
        complex.interactions <- unique.interactions[complex.index, 
                                                    c("int_pair", "geneA", "geneB")] %>%
            transmute(int_pair = int_pair, 
                      geneA.1 = vapply(strsplit(geneA, ","), function(x) x[1], 
                                       character(1)), 
                      geneA.2 = vapply(strsplit(geneA, ","), function(x) x[2],
                                       character(1)),
                      geneA.3 = vapply(strsplit(geneA, ","), function(x) x[3],
                                       character(1)),
                      geneA.4 = vapply(strsplit(geneA, ","), function(x) x[4],
                                       character(1)),
                      geneB.1 = vapply(strsplit(geneB, ","), function(x) x[1],
                                       character(1)), 
                      geneB.2 = vapply(strsplit(geneB, ","), function(x) x[2],
                                       character(1)),
                      geneB.3 = vapply(strsplit(geneB, ","), function(x) x[3],
                                       character(1)),
                      geneB.4 = vapply(strsplit(geneB, ","), function(x) x[4],
                                       character(1))
            )
    } else {
        simple.interactions <- unique.interactions[, c("int_pair", "geneA", "geneB")] 
    }
    
    GO_annotation <- data.table(int_pair=character(), 
                                id=character(), 
                                term=character(), 
                                source=character(), 
                                stringsAsFactors = FALSE)
    
    # annotate simple
    for(i in 1:nrow(simple.interactions)){
        go.ann <- GO.biomart %>%
            filter(gene_symbol == simple.interactions$geneA[i] | 
                       gene_symbol == simple.interactions$geneB[i]) %>%
            select(-go_linkage_type) %>% distinct(.keep_all = TRUE) %>%
            group_by(go_id) %>% filter(n()==2) %>%
            select(-gene_symbol) %>% distinct(.keep_all = TRUE)
        GO_annotation_tmp <- data.table(
            int_pair = simple.interactions$int_pair[i], 
            id = go.ann$go_id, 
            term = go.ann$go_term, 
            source = go.ann$domain)
        GO_annotation <- rbind(GO_annotation, GO_annotation_tmp)
    }
    
    # annotate complex
    if(length(complex.interactions) > 0){
        n_genes <- complex.interactions %>%
            select(geneA.1:geneB.4)
        n_genes$n_gene <- rowSums(!is.na(n_genes), na.rm = TRUE)
        for(i in 1:nrow(complex.interactions)){
            go.ann <- GO.biomart %>%
                filter(gene_symbol == complex.interactions$geneA.1[i] | 
                           gene_symbol == complex.interactions$geneA.2[i] |
                           gene_symbol == complex.interactions$geneA.3[i] |
                           gene_symbol == complex.interactions$geneA.4[i] |
                           gene_symbol == complex.interactions$geneB.1[i] |
                           gene_symbol == complex.interactions$geneB.2[i] |
                           gene_symbol == complex.interactions$geneB.3[i] |
                           gene_symbol == complex.interactions$geneB.4[i] ) %>%
                select(-go_linkage_type) %>% distinct(.keep_all = TRUE) %>%
                group_by(go_id) %>% filter(n()== n_genes$n_gene[i]) %>%
                select(-gene_symbol) %>% distinct(.keep_all = TRUE)
            GO_annotation_tmp <- data.table(
                int_pair = complex.interactions$int_pair[i], 
                id = go.ann$go_id, 
                term = go.ann$go_term, 
                source = go.ann$domain)
            GO_annotation <- rbind(GO_annotation, GO_annotation_tmp)
        }
    }
    
    # remove NA rows
    GO_annotation <- GO_annotation[!is.na(GO_annotation$term),]
    # Rename columns
    colnames(GO_annotation) <- c("int_pair", "GO_id", "functional_term", 
                                 "source")
    # Rename sources
    GO_annotation$source <- plyr::mapvalues(GO_annotation$source, 
                                            from = c("biological_process", 
                                                     "cellular_component", 
                                                     "molecular_function"),
                                            to = c("GO:BP", "GO:CC", "GO:MF"))
    return(GO_annotation)
    
}


#' Annotate pathways for input data
#'
#' @param selected.db pathways sources to use
#' @param input.data filtered input data
#'
#' @return pathways_annotation
#' @importFrom dplyr distinct transmute

annotatePathways <- function(selected.db, input.data){
    
    # get unique pairs of interactions and associated genes
    unique.interactions <- distinct(input.data, int_pair, .keep_all = TRUE)
    # separate simple interactions from complexes
    complex.index <- grepl(",", unique.interactions$geneA) | 
        grepl(",", unique.interactions$geneB)
    # only done if there are complex interactions
    if(length(complex.index) > 0){
        simple.interactions <- unique.interactions[!complex.index, 
                                                   c("int_pair", "geneA", "geneB")] %>%
            transmute(int_pair = int_pair, 
                      geneA = paste0("SYMBOL:", geneA), 
                      geneB = paste0("SYMBOL:", geneB))
        complex.interactions <- unique.interactions[complex.index, 
                                                    c("int_pair", "geneA", "geneB")] %>%
            transmute(int_pair = int_pair, 
                      geneA.1 = paste0("SYMBOL:", unlist(vapply(strsplit(geneA, ","), function(x) x[1],
                                                                character(1)))), 
                      geneA.2 = paste0("SYMBOL:", unlist(vapply(strsplit(geneA, ","), function(x) x[2],
                                                                character(1)))),
                      geneA.3 = paste0("SYMBOL:", unlist(vapply(strsplit(geneA, ","), function(x) x[3],
                                                                character(1)))),
                      geneA.4 = paste0("SYMBOL:", unlist(vapply(strsplit(geneA, ","), function(x) x[4],
                                                                character(1)))),
                      geneB.1 = paste0("SYMBOL:", unlist(vapply(strsplit(geneB, ","), function(x) x[1],
                                                                character(1)))), 
                      geneB.2 = paste0("SYMBOL:", unlist(vapply(strsplit(geneB, ","), function(x) x[2],
                                                                character(1)))),
                      geneB.3 = paste0("SYMBOL:", unlist(vapply(strsplit(geneB, ","), function(x) x[3],
                                                                character(1)))),
                      geneB.4 = paste0("SYMBOL:", unlist(vapply(strsplit(geneB, ","), function(x) x[4],
                                                                character(1)))),
            )
    } else {
        simple.interactions <- unique.interactions[, c("int_pair", "geneA", "geneB")] %>%
            transmute(int_pair = int_pair, 
                      geneA = paste0("SYMBOL:", geneA), 
                      geneB = paste0("SYMBOL:", geneB))
    }
    
    
    pathways_annotation <- data.table(int_pair=character(), 
                                      term= character(), 
                                      source=character(), 
                                      stringsAsFactors = FALSE)
    for(db.name in selected.db){
        
        switch(db.name,
               biocarta = {
                   db.symbol <- hs_biocarta
               },
               kegg = {
                   db.symbol <- hs_kegg
               },
               nci = {
                   db.symbol <- hs_nci
               },
               panther = {
                   db.symbol <- hs_panther
               },
               pharmgkb = {
                   db.symbol <- hs_pharmgkb
               },
               reactome = {
                   db.symbol <- hs_reactome
               }
               )
        
        for(p in 1:length(db.symbol)){
            # annotate simple
            int.included <- simple.interactions$int_pair[
                simple.interactions$geneA %in% db.symbol[[p]] & 
                    simple.interactions$geneB %in% db.symbol[[p]]]
            if(length(int.included)>0){
                pathways_tmp <- data.table(int_pair = int.included, 
                                           term = names(db.symbol[p]), 
                                           source = db.name, 
                                           stringsAsFactors = FALSE)
                pathways_annotation <- rbind(pathways_annotation, pathways_tmp)
            }
            # annotate complex
            if(length(complex.interactions) > 0){
                int.included <- complex.interactions$int_pair[complex.interactions$geneA.1 %in% c("SYMBOL:NA", db.symbol[[p]])
                                                              & complex.interactions$geneA.2 %in% c("SYMBOL:NA", db.symbol[[p]])
                                                              & complex.interactions$geneA.3 %in% c("SYMBOL:NA", db.symbol[[p]])
                                                              & complex.interactions$geneA.4 %in% c("SYMBOL:NA", db.symbol[[p]])
                                                              & complex.interactions$geneB.1 %in% c("SYMBOL:NA", db.symbol[[p]])
                                                              & complex.interactions$geneB.2 %in% c("SYMBOL:NA", db.symbol[[p]])
                                                              & complex.interactions$geneB.3 %in% c("SYMBOL:NA", db.symbol[[p]])
                                                              & complex.interactions$geneB.4 %in% c("SYMBOL:NA", db.symbol[[p]])]
                if(length(int.included)>0){
                    pathways_tmp <- data.table(int_pair = int.included, 
                                               term = names(db.symbol[p]), 
                                               source = db.name, 
                                               stringsAsFactors = FALSE)
                    pathways_annotation <- rbind(pathways_annotation, 
                                                 pathways_tmp)
                }
            }
        }
        
    }
    # Rename columns
    colnames(pathways_annotation) <- c("int_pair", "functional_term", "source")
    return(pathways_annotation)
    
}


#' Combine GO annotation and pathways in a unique object
#'
#' @param GO_annotation data
#' @param pathways_annotation  data
#'
#' @return combined annotation dataframe
#' @importFrom tibble add_column
#' @importFrom dplyr arrange group_by filter

combineAnnotations <- function(GO_annotation, pathways_annotation){
    # Add empty for GO_id in pathways
    pathways_annotation <- add_column(pathways_annotation, GO_id = "", 
                                      .after = "int_pair")
    combined <- rbind(GO_annotation, pathways_annotation)
    combined <- combined %>%
        arrange(int_pair, functional_term)
    
    dup <- combined %>%
        group_by(int_pair) %>%
        filter(duplicated(tolower(functional_term)))
    combined <- combined %>%
        group_by(int_pair) %>%
        filter(!duplicated(tolower(functional_term)))
    for(r in 1:nrow(dup)){
        dup_row <- paste0(dup$int_pair[r], tolower(dup$functional_term[r]))
        comb_rows <- paste0(combined$int_pair, tolower(combined$functional_term))
        ind_dup <- match(dup_row, comb_rows)
        combined[ind_dup, "source"] <- paste(combined[ind_dup, "source"], 
                                             dup$source[r], sep = ",")
    }
    
    return(combined)
}

#' Calculate number of terms of a database
#'
#' @param annotation data from either pathways, GO or combined
#'
#' @return number of terms by dataset
#' @importFrom dplyr group_by n summarise

getNtermsBYdb <- function(annotation){
    nTermsBYdataset <- annotation %>%
        group_by(source) %>%
        summarise(n_terms = n())
    return(nTermsBYdataset)
}


#' Build binary matrix with int-pairs in rows, functions in cols
#'
#' @param functions_df annotated df (GO/path/combined)
#'
#' @return binary matrix
#' @importFrom data.table dcast as.data.table

buildPairsbyFunctionMatrix <- function(functions_df){
    functions_df$value <- 1
    pairs_func_matrix <- as.matrix(dcast(as.data.table(functions_df), 
                                         int_pair ~ functional_term, 
                                         fun = function(x) min(sum(x), 1), 
                                         value.var = "value"))
    gene_pairs <- pairs_func_matrix[, "int_pair"]
    pairs_func_matrix <- pairs_func_matrix[, -1]
    rownames(pairs_func_matrix) <- gene_pairs
    storage.mode(pairs_func_matrix) <- "numeric"
    return(pairs_func_matrix)
}

#' Get GO link
#'
#' @param go_id string
#'
#' @return html link to website
#'
goLink <- function(go_id){
    paste0('<a href="http://amigo.geneontology.org/amigo/term/', go_id, 
           '" target="_blank">', go_id, '</a>')
}


#' Get table with ranked functional terms
#'
#' @param data.fun.annot annotated df (GO/path/combined)
#' @param gene.table of unique intpairs
#' 
#' @return table with ranking
#' @importFrom dplyr group_by summarise arrange

getRankedTerms <- function(data.fun.annot, gene.table){
    data.fun.annot$uniq_score <- gene.table$uniqueness_score[
        match(data.fun.annot$int_pair, gene.table$int_pair)]
    rank.terms <- data.fun.annot %>%
        group_by(tolower(functional_term)) %>%
        summarise(n_occurrence = n(),
                  # Add average uniqueness score
                  avg_uniqueness = round(mean(uniq_score), digits = 2),
                  int_pair_list = paste(int_pair, collapse = ","), 
                  source = paste(source, collapse = ",")) %>%
        arrange(plyr::desc(avg_uniqueness))
    colnames(rank.terms)[1] <- "functional_term"
    rank.terms$source <- sapply(rank.terms$source, 
                                   function(x) paste(unique(
                                       unlist(strsplit(x, split=","))), 
                                       collapse = ","))
    
    
    
    return(rank.terms)
}

#' Get Sunburst plot of selected functional terms
#'
#' @param sel.data dataframe of selected functions
#' @param func_selected the selected functional term
#' @param int_p_fun dataframe with int pairs annotated to this function
#' @param cluster.colors for plotting
#'
#' @return plotly figure
#' @importFrom dplyr group_by summarise n
#' @importFrom plotly plot_ly
getSunburst <- function(sel.data, func_selected, int_p_fun, cluster.colors){
    trace1 <- sel.data %>%
        group_by(clustA, int_pair) %>%
        dplyr::summarise(n = n())
    
    trace1 <- trace1 %>%
        group_by(clustA) %>%
        dplyr::summarise(n_tot = sum(n))
    
    trace2 <- sel.data %>%
        group_by(clustA, clustB) %>%
        dplyr::summarise(n = n(), int_pair = paste(int_pair, collapse = ","))
    
    sunburst.df <- data.frame(
        parents = c(rep(func_selected, nrow(trace1)), 
                    paste0("tr1_",trace2$clustA)), 
        ids = c(paste0("tr1_", trace1$clustA), 
                paste0("tr2_",trace2$clustA, trace2$clustB)),
        labels = c(trace1$clustA, trace2$clustB), 
        values = c(trace1$n_tot, trace2$n),
        text = c(paste(trace1$n_tot, nrow(sel.data), sep="/"), 
                 trace2$int_pair),
        stringsAsFactors = FALSE)
    
    
    
    # insert break lines in functional term only
    ind.func_sel <- which(sunburst.df$parents == func_selected)
    sunburst.df$parents[ind.func_sel] <- gsub(" ", "<br>", 
                                              sunburst.df$parents[ind.func_sel])
    sunburst.df$text <- gsub(",", "<br>", sunburst.df$text)
    sunburst.df$colors <- cluster.colors[sunburst.df$labels]
    fig <- plot_ly(sunburst.df, ids = ~ids, labels = ~labels, 
            parents = ~parents, values = ~values, 
            hovertext = ~text, type = 'sunburst', 
            hoverinfo = "text", marker = list(colors= ~colors))
    return(fig)
}
