#' Subset pairs-function matrix by selected flow
#'
#' @param pairs_func_matrix binary
#' @param flow_df subset of input data by flow
#'
#' @return subset of binary mat
#' @export
subsetFuncMatBYFlow <- function(pairs_func_matrix, flow_df){
    sub.mat <- pairs_func_matrix[rownames(pairs_func_matrix) %in% 
                                     unique(flow_df$int_pair),]
    #remove empty columns
    sub.mat <- sub.mat[, colSums(sub.mat) != 0]
    return(sub.mat)
}


#' Get dendrogram of int pair modules
#'
#' @param pairs_func_matrix binary matrix pairs x functions
#' @param seed for reproducibility
#'
#' @return list with dendrogram, hclust and umap
#' @importFrom umap umap
#' @importFrom stats hclust dist
#' @export

dendroIntPairModules <- function(pairs_func_matrix, seed = NULL){
    if(!is.null(seed)){
        set.seed(seed)
    }
    intPairs_umap <- umap(pairs_func_matrix, 
                          n_neighbors = ifelse(nrow(pairs_func_matrix) > 15,
                                               15,
                                               nrow(pairs_func_matrix)-1), 
                          n_components = 2,
                          metric = "cosine", input= "data", min_dist = 0.001)
    umap.embed <- data.frame(UMAP_1 = intPairs_umap$layout[,1], 
                             UMAP_2 = intPairs_umap$layout[,2],
                             int_pair = dimnames(intPairs_umap$layout)[[1]])
    
    ## Hierarchical clust
    d <- dist(umap.embed[, c("UMAP_1", "UMAP_2")], method="euclidean")
    h_clust <- hclust(d, method = "ward.D2")
    
    return(list(d = d, h_clust = h_clust, umap = umap.embed))
}


#' Determine the elbow point on a curve (from package akmedoids)
#' @description Given a list of x, y coordinates on a curve, function determines the elbow point of the curve.
#' 
#' @param x vector of x coordinates of points on the curve
#' @param y vector of y coordinates of points on the curve
#' 
#' @details highlight the maximum curvature to identify the elbow point (credit: 'github.com/agentlans')

#' @return an x, y coordinates of the elbow point.
#' @importFrom stats approx approxfun optimize predict smooth.spline
#' @importFrom signal sgolayfilt

elbowPoint <- function(x, y) {
    
    # check for non-numeric or infinite values in the inputs
    is.invalid <- function(x) {
        any((!is.numeric(x)) | is.infinite(x))
    }
    if (is.invalid(x) || is.invalid(y)) {
        stop("x and y must be finite and numeric. Missing values are not allowed.")
    }
    if (length(x) != length(y)) {
        stop("x and y must be of equal length.")
    }
    
    # generate value of curve at equally-spaced points
    new.x <- seq(from=min(x), to=max(x), length.out=length(x))
    # Smooths out noise using a spline
    sp <- smooth.spline(x, y)
    new.y <- predict(sp, new.x)$y
    
    # Finds largest odd number below given number
    largest.odd.num.lte <- function(x) {
        x.int <- floor(x)
        if (x.int %% 2 == 0) {
            x.int - 1
        } else {
            x.int
        }
    }
    
    # Use Savitzky-Golay filter to get derivatives
    smoothen <- function(y, p=p, filt.length=NULL, ...) {
        # Time scaling factor so that the derivatives are on same scale as original data
        ts <- (max(new.x) - min(new.x)) / length(new.x)
        p <- 3 # Degree of polynomial to estimate curve
        # Set filter length to be fraction of length of data
        # (must be an odd number)
        if (is.null(filt.length)) {
            filt.length <- min(largest.odd.num.lte(length(new.x)), 7)
        }
        if (filt.length <= p) {
            stop("Need more points to find cutoff.")
        }
        signal::sgolayfilt(y, p=p, n=filt.length, ts=ts, ...)
    }
    
    # Calculate first and second derivatives
    first.deriv <- smoothen(new.y, m=1)
    second.deriv <- smoothen(new.y, m=2)
    
    # Check the signs of the 2 derivatives to see whether to flip the curve
    # (Pick sign of the most extreme observation)
    pick.sign <- function(x) {
        most.extreme <- which(abs(x) == max(abs(x), na.rm=TRUE))[1]
        sign(x[most.extreme])
    }
    first.deriv.sign <- pick.sign(first.deriv)
    second.deriv.sign <- pick.sign(second.deriv)
    
    # The signs for which to flip the x and y axes
    x.sign <- 1
    y.sign <- 1
    if ((first.deriv.sign == -1) && (second.deriv.sign == -1)) {
        x.sign <- -1
    } else if ((first.deriv.sign == -1) && (second.deriv.sign == 1)) {
        y.sign <- -1
    } else if ((first.deriv.sign == 1) && (second.deriv.sign == 1)) {
        x.sign <- -1
        y.sign <- -1
    }
    # If curve needs flipping, then run same routine on flipped curve then
    # flip the results back
    if ((x.sign == -1) || (y.sign == -1)) {
        results <- elbowPoint(x.sign * x, y.sign * y)
        return(list(x = x.sign * results$x, y = y.sign * results$y))
    }
    
    # Find cutoff point for x
    cutoff.x <- NA
    # Find x where curvature is maximum
    curvature <- abs(second.deriv) / (1 + first.deriv^2)^(3/2)
    
    if (max(curvature) < min(curvature) | max(curvature) < max(curvature)) {
        cutoff.x = NA
    } else {
        # Interpolation function
        f <- approxfun(new.x, curvature, rule=1)
        # Minimize |f(new.x) - max(curvature)| over range of new.x
        cutoff.x = optimize(function(new.x) abs(f(new.x) - max(curvature)), range(new.x))$minimum
    }
    
    if (is.na(cutoff.x)) {
        warning("Cutoff point is beyond range. Returning NA.")
        list(x=NA, y=NA)
    } else {
        # Return cutoff point on curve
        approx(new.x, new.y, cutoff.x)
    }
}

#' Get UMAP for IP modules
#'
#' @param intPairs.dendro list output of dendrogram
#' @param gpModules_assign named vector of module assignment
#' @param gene.table unique intpairs table
#' @param ipm_colors for intpair modules
#' @param input_ipM_UMAPcolors user choice for coloring umap
#'
#' @return plotly umap
#' @importFrom plotly plot_ly layout config
getUMAPipModules <- function(intPairs.dendro, 
                             gpModules_assign, 
                             gene.table,
                             ipm_colors, 
                             input_ipM_UMAPcolors){
    umap.embed <- intPairs.dendro$umap
    umap.embed$hclust <- as.factor(gpModules_assign[
        match(umap.embed$int_pair, names(gpModules_assign))])
    umap.embed$uniq_score <- gene.table$uniqueness_score[
        match(umap.embed$int_pair, gene.table$int_pair)]
    
    # Choose colors
    if(input_ipM_UMAPcolors == "ipM_col"){
        colors <- ipm_colors
        color_var <- "hclust"
    } else if(input_ipM_UMAPcolors == "US_col"){
        colors <- scales::viridis_pal()(10)
        color_var <- "uniq_score"
    }
    
    ax <- list(zeroline=FALSE)
    fig <- plot_ly(data = umap.embed, 
                   x= ~UMAP_1, y= ~UMAP_2, 
                   type='scatter', mode='markers', 
                   color = umap.embed[, color_var],
                   text = ~as.character(int_pair), 
                   hoverinfo='text', colors = colors)
    fig <- fig %>% layout(xaxis = ax, yaxis = ax, 
                          title="<b>UMAP of Int-pairs</b>")
    fig <- fig %>% config(modeBarButtonsToRemove = c(
        'sendDataToCloud', 'autoScale2d', 'resetScale2d',
        'hoverClosestCartesian', 'hoverCompareCartesian',
        'zoom2d','pan2d','select2d','lasso2d'))
    return(fig)
}


#' Plot circle plot
#'
#' @param data subset of input data by flow / intpair module
#' @param cluster_colors global
#' @param int_flow string specifying the flow 
#' @param link.color of intpair module
#'
#' @return circle plot
#' 
#' @importFrom circlize circos.par chordDiagram circos.trackPlotRegion 
#' get.cell.meta.data circos.text highlight.sector circos.clear uh CELL_META

circlePlot <- function(data, cluster_colors, int_flow, link.color){
    
    cell_types <- unique(c(data$clustA, data$clustB))
    partnerA <- unlist(sapply(strsplit(data$int_pair, " & "), function(x) x[1]))
    partnerB <- unlist(sapply(strsplit(data$int_pair, " & "), function(x) x[2]))
    genes <- c(structure(partnerA, names = data$clustA), 
               structure(partnerB, names = data$clustB))
    genes <- genes[!duplicated(paste(names(genes), genes))]
    genes <- genes[order(names(genes))]
    
    
    
    if(length(cell_types)!=1){
        gap.degree <- do.call("c", lapply(table(names(genes)), 
                                          function(i) c(rep(1, i-1), 8)))
    }else{
        gap.degree <- do.call("c", lapply(table(names(genes)), 
                                          function(i) c(rep(1, i))))
    }
    
    # parameters
    if(int_flow == "undirected"){
        directional <- 0
        direction.type <- "diffHeight"
    } else{
        directional <- 1
        direction.type <- c("diffHeight", "arrows")
    }
    
    track.height.genes <- ifelse(max(nchar(c(partnerA, partnerB))) >= 10, 
                                 0.25, 
                                 0.2)
    cex.genes <- 0.9

    
    
    df <- data.frame(from = paste(data$clustA,partnerA), 
                     to = paste(data$clustB,partnerB), 
                     stringsAsFactors = FALSE)
    
    
    circos.par(gap.degree = gap.degree)
    
    chordDiagram(df, order=paste(names(genes),genes),
                 grid.col = link.color, 
                 transparency = 0.6, 
                 directional = directional, 
                 direction.type = direction.type,
                 link.arr.type = "big.arrow", 
                 annotationTrack = "grid", 
                 preAllocateTracks = list(
                     list(track.height = uh(1.2,'mm')), 
                     list(track.height = track.height.genes)),  
                 annotationTrackHeight = c(0.01,0.01))
    
    
    
    
    circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
        sector.index = genes[get.cell.meta.data("sector.numeric.index")]
        circos.text(CELL_META$xcenter, 
                    CELL_META$cell.ylim[1], 
                    sector.index, 
                    col = "black", 
                    cex = cex.genes, 
                    adj = c(0, 0.5), 
                    facing = 'clockwise', 
                    niceFacing = TRUE)
    }, bg.border = NA)
    
    
    for(c in unique(names(genes))) {
        gene = as.character(genes[names(genes) == c])
        highlight.sector(sector.index = paste(c,gene), 
                         track.index = 1, 
                         col = ifelse(length(cluster_colors)==1,
                                      cluster_colors,
                                      cluster_colors[c]), 
                         text = c, 
                         text.vjust = '0.4cm', 
                         niceFacing = TRUE, 
                         lwd=1,
                         facing = "bending.inside")
    }
    
    circos.clear()
    
    
    
}


#' Subfunction to calculate significant functions by permutation test
#'
#' @param mat binary matrix of intpairs by functions
#' @param gpModules_assign assignment of intpairs to modules
#'
#' @return matrix with hits

getHitsf <- function(mat, gpModules_assign){
    hits <- matrix(0, nrow = nrow(mat), ncol = length(unique(gpModules_assign)))
    rownames(hits) <- rownames(mat)
    colnames(hits) <- unique(gpModules_assign)
    for(gi in unique(gpModules_assign)){
        sub.mat <- mat[, names(gpModules_assign)[gpModules_assign == gi]]
        hits[, gi] <- rowSums(sub.mat)/ncol(sub.mat)
    }
    return(hits)
}


#' Calculate significant function per intpair module
#'
#' @param subGenePairs_func_mat subset of binary mat
#' @param gpModules_assign assignment of intpairs to modules
#' @param rank.terms table of ranked functions
#' @param input_maxPval threshold of significance
#'
#' @return table with significant functions
#' @importFrom tidyr gather


getSignificantFunctions <- function(subGenePairs_func_mat, 
                                    gpModules_assign,
                                    rank.terms,
                                    input_maxPval){
    permMat <- t(subGenePairs_func_mat)
 
    hits_true <- getHitsf(permMat, gpModules_assign)
    hits_perm <- list()
    
    for(np in 1:999){
        # shuffle cols of original matrix
        shufMat <- permMat[,sample(colnames(permMat), ncol(permMat), 
                                   replace = FALSE)]
        colnames(shufMat) <- colnames(permMat)
        hits_perm[[np]] <- getHitsf(shufMat, gpModules_assign)
    }
    
    # calculate empirical pvalue
    pvalue <- matrix(0, nrow = nrow(permMat), 
                     ncol = length(unique(gpModules_assign)))
    rownames(pvalue) <- rownames(permMat)
    colnames(pvalue) <- unique(gpModules_assign)
    for(gM in 1:ncol(hits_true)){
        for(fM in 1:nrow(hits_true)){
            hits_gm_fm <- unlist(lapply(hits_perm, function(x) x[fM, gM]))
            pvalue[fM,gM] <- (1 + sum(hits_gm_fm >= hits_true[fM,gM]))/1000
        }
    }
    
    pvalue_df <- cbind(pvalue, functionalTerm = rownames(pvalue))
    pvalue_df <- gather(as.data.frame(pvalue_df), 
                        key = "int_pairModule", 
                        value = "pvalue", 
                        unique(gpModules_assign), 
                        factor_key = FALSE)
    
    signFun <- pvalue_df[pvalue_df$pvalue <= input_maxPval,]
    
    ## Adding int_pairs from selected Module to each functional term
    
    for(r in 1:nrow(signFun)){
        int_pairs_all <- rownames(subGenePairs_func_mat)[
            subGenePairs_func_mat[, signFun$functionalTerm[r]] == 1]
        signFun[r, "int_pair_list"] <- paste(
            intersect(int_pairs_all, names(gpModules_assign)[
                gpModules_assign == signFun$int_pairModule[r]]), collapse = ",")
    }
 
    genes_all <- rownames(subGenePairs_func_mat)[subGenePairs_func_mat[, signFun$fTerm[1]] == 1]
    paste(intersect(genes_all, names(gpModules_assign)[gpModules_assign == 1]), collapse = ",")
    signFun$source <- rank.terms$source[
        match(tolower(signFun$functionalTerm), 
              tolower(rank.terms$functional_term))]
    signFun$avg_uniqueness <- rank.terms$avg_uniqueness[
        match(tolower(signFun$functionalTerm), 
              tolower(rank.terms$functional_term))]
    return(signFun)
}



#' Construct table with occurrence values for significant terms
#'
#' @param signFun table of significant terms
#' @param ipMselected number of intpair module selected by user
#' @param gpModules_assign assignment of intpairs to modules
#' @param subGenePairs_func_mat subset of binary matrix
#'
#' @return table of occurrence

getOccurrenceTab4wordcloud <- function(signFun, 
                                       ipMselected, 
                                       gpModules_assign, 
                                       subGenePairs_func_mat){
    functionSelected <- as.character(unlist(
        signFun[signFun$`int_pairModule` == as.numeric(ipMselected), "functionalTerm"]))
    # create table of occurrence for these terms based on int-pair module selected
    gpSelected <- names(gpModules_assign)[
        gpModules_assign == as.numeric(ipMselected)]
    t_occurrence <- subGenePairs_func_mat[gpSelected, functionSelected]
    t_occurrence <- data.frame(term = functionSelected, 
                               occ = colSums(t_occurrence))

    return(t_occurrence)
}

#' Plot wordcloud of significant terms
#'
#' @param occurTab table with occurrences of terms
#' @param input_minFreq_wordcloud threshold on min occurrence by user
#'
#' @return wordloud2 html object
#' @importFrom wordcloud2  wordcloud2
plotWordCloud <- function(occurTab, input_minFreq_wordcloud){
    colnames(occurTab) <- c("word", "freq")
    occurTab <- occurTab %>%
        filter(freq >= input_minFreq_wordcloud)
    wordcloud2(data = occurTab, 
               shuffle = FALSE, 
               shape = "diamond", 
               color = "random-dark",
               fontWeight = "normal",
               size = .2)
    
}

#' Function to fix the html generated by wordcloud2
#'
#' @param inputFile html file 
#' @param outputFile html file 
#'
#' @return fixed html

wordcloud2FixHTML <- function(inputFile, outputFile){
    a = readLines(inputFile)
    output = paste(a, collapse = "\n")
    output = gsub(">\n\n</div>","></div>",output)
    writeLines(output, outputFile)
    invisible(NULL)
}
