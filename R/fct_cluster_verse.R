
#' Get clusters names from initial input data
#'
#' @param input.data preprocessed input data
#'
#' @return named list of clusters

getClusterNames <- function(input.data){
    cl <- unique(c(input.data$clustA, input.data$clustB))
    # check if cluster names are numbers
    if(all(!grepl("\\D", cl))){
        cl <- cl[order(as.numeric(cl))]
    } else {
        cl <- cl[order(cl)]
    }
    cluster.list <- as.list(cl)
    names(cluster.list) <- unlist(cluster.list)
    return(cluster.list)
}



#' Creating edges dataframe for network of clusters
#'
#' @param input.data  preprocessed input data
#'
#' @return edges dataframe
#' @importFrom scales rescale

getClusterNetwork <- function(input.data){
    edges.df <- input.data %>%
        dplyr::select(clustA, clustB) %>%
        dplyr::group_by(clustA, clustB) %>%
        dplyr::summarise(nInteractions = dplyr::n())
    edges.df$width <- rescale(edges.df$nInteractions, 
                              to = c(0,5), 
                              from = c(0, max(edges.df$nInteractions)))
    return(edges.df)
}

#' Get Clusters size
#'
#' @param cl cluster name
#' @param edges.df dataframe with edges for network
#'
#' @return sum of interactions for that cluster

getClusterSize <- function(cl, edges.df){
    edges_sub <- edges.df %>% 
        filter(clustA == cl | clustB == cl) 
    
    return(sum(edges_sub$nInteractions))
}

#' Create Network of clusters
#'
#' @param data.filt.cluster filtered input data (by clusters)
#'
#' @return list containing nodes and edges for network
#' @importFrom scales hue_pal

createNetwork <- function(data.filt.cluster){
    # Get cluster names
    cluster.list <- getClusterNames(data.filt.cluster)
    edges.filt <- getClusterNetwork(data.filt.cluster)
    cluster.size <- lapply(cluster.list, function(x) 
        getClusterSize(x, edges.filt))
    # Control shape of nodes depending on name type
    if(all(!grepl("\\D", names(cluster.list)))){
        shape <- "circle" #numbers, text inside
    } else {
        shape <- "dot" # letters, text outside
    }
    nodes <- data.frame(id = names(cluster.list), shape = shape,
                        label = names(cluster.list),
                        value = as.numeric(as.vector(cluster.size)),
                        title = paste0("<p>", as.numeric(as.vector(cluster.size)), "</p>"),
                        color = hue_pal(c = 80, l = 80)(length(names(cluster.list))))
    edges <- data.frame(from = edges.filt$clustA, to = edges.filt$clustB,
                        width = edges.filt$width, label = edges.filt$nInteractions, 
                        length = 40*length(names(cluster.list)),
                        arrows.to.type = "arrow")
    return(list(nodes = nodes, edges = edges))
}


#' Create Barplot cluster-verse
#'
#' @param data.filt.bar filtered object (checkbox auto/para)
#' @param input_cluster_selected_checkbox checkbox input
#'
#' @return plotly barplot
#' @importFrom plotly plot_ly add_trace layout config
#' @importFrom dplyr filter
#' @importFrom scales hue_pal

createBarPlot_CV <- function(data.filt.bar, input_cluster_selected_checkbox){
    tot.int <- data.frame(clusters = input_cluster_selected_checkbox)
    tot.int$n_paracrine <- lapply(input_cluster_selected_checkbox, 
                                  function(x)  nrow(data.filt.bar %>%
                                                filter(clustA == x | clustB == x) %>%
                                                filter(int.type == "paracrine")))
    tot.int$n_autocrine <- lapply(input_cluster_selected_checkbox, 
                                  function(x)  nrow(data.filt.bar %>%
                                                filter(clustA == x & clustB == x)))
    # make clusters factor to have ordered x axis
    tot.int$clusters <- factor(tot.int$clusters, 
                               levels = input_cluster_selected_checkbox)
    
    cluster.colors <- hue_pal(c = 80, l = 80)(length(input_cluster_selected_checkbox))
    fig <- plot_ly(tot.int, x = ~clusters, y = ~n_paracrine, type = "bar", 
                   name = "paracrine", 
                   marker = list(line = list(color = cluster.colors, width = 3),
                                 color = "#C0C0C0"))
    fig <- fig %>% add_trace(y = ~n_autocrine, name= "autocrine", 
                             marker = list(color = "#606060") )
    fig <- fig %>% layout(title = "Total number of interactions per cluster",
                          xaxis = list(title = "Clusters"),
                          yaxis = list(title = "# Interactions"),
                          barmode = "stack")
    fig <- fig %>% config(modeBarButtonsToRemove = c(
                              'sendDataToCloud', 'autoScale2d', 'resetScale2d', 
                              'hoverClosestCartesian', 'hoverCompareCartesian',
                              'zoom2d','pan2d','select2d','lasso2d'
                          ))
    return(fig)
}



#' Get subset of interactions corresponding to a certain viewpoint and flow
#'
#' @param vp viewpoint cluster
#' @param input.data preprocessed/filtered input data
#' @param flow one among directed_out, directed_in or undirected
#'
#' @return subset of data
#' @importFrom dplyr filter
#' @export
getIntFlow <- function(vp, input.data, flow){
    RRint <- input.data %>%
        filter(typeA == "R" & typeB == "R")
    LLint <- input.data %>%
        filter(typeA == "L" & typeB == "L")
    LRint <- input.data %>%
        filter(typeA == "L" & typeB == "R")
    
    directed.int <- LRint
    undirected.int <- rbind(RRint, LLint)
    
    if(flow == "directed_out"){
        # 1. Directed, outgoing (L on vp)
        data1 <- directed.int %>%
            filter(clustA == vp)
        return(data1)
    }
    else if(flow == "directed_in"){
        # 2. Directed, incoming (R on vp)
        data2 <- directed.int %>%
            filter(clustB == vp)
        return(data2)
    }
    else if(flow == "undirected"){
        # 3. Undirected, autocrine and paracrine
        data3 <- undirected.int %>%
            filter(clustA == vp | clustB == vp)
        return(data3)
    }
 
}