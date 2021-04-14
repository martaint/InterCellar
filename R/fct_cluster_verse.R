
#' Get clusters names from initial input data
#'
#' @param input.data preprocessed input data
#'
#' @return named list of clusters
#' @export
#' @examples 
#' data(input.data)
#' cluster_list <- getClusterNames(input.data)

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



#' Get dataframe for plotting barplot (all clusters)
#'
#' @param data.filt.bar filtered object (checkbox auto/para)
#' @param input_cluster_selected_checkbox checkbox input
#'
#' @return dataframe with number of interactions per cluster auto/para
#' @importFrom dplyr filter
getBarplotDF <- function(data.filt.bar, input_cluster_selected_checkbox){
    barplotDF <- data.frame(clusters = input_cluster_selected_checkbox)
    barplotDF$n_paracrine <- unlist(lapply(input_cluster_selected_checkbox, 
                                  function(x)  nrow(data.filt.bar %>%
                                                        filter(clustA == x | clustB == x) %>%
                                                        filter(int.type == "paracrine"))))
    barplotDF$n_autocrine <- unlist(lapply(input_cluster_selected_checkbox, 
                                  function(x)  nrow(data.filt.bar %>%
                                                        filter(clustA == x & clustB == x))))
    # make clusters factor to have ordered x axis
    barplotDF$clusters <- factor(barplotDF$clusters,
                                 levels = input_cluster_selected_checkbox)
    return(barplotDF)
}

#' Create Barplot cluster-verse
#'
#' @param barplotDF dataframe with N interactions per cluster (auto/para)
#' @param input_cluster_selected_checkbox checkbox input
#'
#' @return plotly barplot
#' @importFrom plotly plot_ly add_trace layout config
#' @importFrom dplyr filter
#' @importFrom scales hue_pal

createBarPlot_CV <- function(barplotDF, input_cluster_selected_checkbox){
    
    
    cluster.colors <- hue_pal(c = 80, l = 80)(length(input_cluster_selected_checkbox))
    fig <- plot_ly(barplotDF, x = ~clusters, y = ~n_paracrine, type = "bar", 
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


#' Create ggplot barplot to be saved in tiff
#'
#' @param barplotDF dataframe with N interactions per cluster (auto/para)
#' @param input_cluster_selected_checkbox checkbox input 
#'
#' @return ggplot barplot
#' @importFrom tidyr pivot_longer
#' @importFrom scales hue_pal

createBarPlot1_ggplot <- function(barplotDF, input_cluster_selected_checkbox){
    bar_ggplot_df <- tidyr::pivot_longer(barplotDF, 
                                         cols = c("n_paracrine","n_autocrine"), 
                                         names_to = "type", 
                                         names_prefix = "n_",
                                         values_to = "n_int")
    cluster.colors <- scales::hue_pal(c = 80, l = 80)(length(input_cluster_selected_checkbox))
    
    g <- ggplot(bar_ggplot_df, aes(x = clusters, 
                                   y = n_int, 
                                   fill = type, 
                                   color = clusters)) +
        geom_bar(position="stack", stat="identity") +
        theme_minimal() +
        scale_fill_manual(values = c("#606060", "#C0C0C0")) + 
        scale_color_manual(values = cluster.colors) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
              text = element_text(size=20),
              strip.text = element_blank()) +
        guides(color = FALSE) + 
        labs(x = "Clusters", y = "# Interactions", 
             title = "Total number of interactions per cluster") 
    return(g)
}

#' Get dataframe for barplot (by cluster)
#'
#' @param filt.data input data filtered in cluster-verse
#' @param input_cluster_selected_checkbox selected clusters to keep
#' @param input_clust_barplot2 selected cluster to plot
#'
#' @return dataframe with num int per cluster
#' @importFrom dplyr filter
getBarplotDF2 <- function(filt.data, input_cluster_selected_checkbox,
                          input_clust_barplot2){
    # Get paracrine and autocrine data
    para <- filt.data %>%
        filter(clustA == input_clust_barplot2 | clustB == input_clust_barplot2) %>%
        filter(!(clustA == input_clust_barplot2 & clustB == input_clust_barplot2))
    auto <- filt.data %>%
        filter(clustA == input_clust_barplot2 & clustB == input_clust_barplot2)
    
    bar.data <- data.frame(Clusters = input_cluster_selected_checkbox, 
                           Num_int = 0)
    for(c in input_cluster_selected_checkbox){
        if(c == input_clust_barplot2){
            bar.data$Num_int[bar.data$Clusters == c] <- nrow(auto)
        } else {
            bar.data$Num_int[bar.data$Clusters == c] <- sum(para$clustA == c) +
                sum(para$clustB == c)
        }
        
    }
    
    bar.data$Clusters <- factor(bar.data$Clusters, 
                                levels = input_cluster_selected_checkbox)
    return(bar.data)
}

#' Create barplot of number of interaction for selected cluster
#'
#' @param barplotDF2 dataframe with barplot data
#' @param input_cluster_selected_checkbox selected clusters to keep
#' @param input_clust_barplot2 selected cluster to plot
#'
#' @return plotly fig
#' @importFrom plotly plot_ly add_annotations layout config

createBarPlot2_CV <- function(barplotDF2, input_cluster_selected_checkbox,
                              input_clust_barplot2){
    
    cluster.colors <- hue_pal(c = 80, l = 80)(length(input_cluster_selected_checkbox))
    fig <- plot_ly(barplotDF2, 
                   x = ~Clusters, y = ~Num_int, type = "bar",
                   marker = list(color = cluster.colors),
                   text = ~Num_int,
                   textposition = "outside") %>%
        layout(title = paste0("Number of Interactions for Cluster ", 
                              input_clust_barplot2),
                          xaxis = list(title = "Clusters"),
                          yaxis = list(title = "# Interactions"))
    fig <- fig %>% config(modeBarButtonsToRemove = c(
        'sendDataToCloud', 'autoScale2d', 'resetScale2d', 
        'hoverClosestCartesian', 'hoverCompareCartesian',
        'zoom2d','pan2d','select2d','lasso2d'
    ))
    return(fig)
    
    
}


#' Create ggplot barplot of Nint per cluster selected
#'
#' @param barplotDF2 dataframe with barplot data
#' @param input_cluster_selected_checkbox selected clusters to keep
#' @param input_clust_barplot2 selected cluster to plot
#'
#' @return ggplot barplot
#' @importFrom ggplot2 geom_bar geom_text theme_minimal scale_fill_manual
#' theme guides labs
createBarPlot2_ggplot <- function(barplotDF2, input_cluster_selected_checkbox,
                                  input_clust_barplot2){
    
    cluster.colors <- scales::hue_pal(c = 80, l = 80)(length(input_cluster_selected_checkbox))
    
    g <- ggplot(barplotDF2, aes(x = Clusters, 
                                y = Num_int, 
                                fill = Clusters)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = Num_int), vjust = -0.3, size = 5)+
        theme_minimal() +
        scale_fill_manual(values = cluster.colors) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
              text = element_text(size=20),
              strip.text = element_blank()) +
        guides(fill = FALSE) + 
        labs(x = "Clusters", y = "# Interactions", 
             title = paste0("Number of Interactions for Cluster ", 
                            input_clust_barplot2)) 
    return(g)
    
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
#' @examples 
#' data(input.data)
#' caf_out <- getIntFlow(vp = "CAF", input.data, flow = "directed_out")
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

#' Get back-to-back barplot for 2 conditions comparison
#'
#' @param tab_c1 table from csv file (barplot#1) containing data for condition 1
#' @param tab_c2 table from csv file (barplot#1)containing data for condition 2
#' @param lab_c1 label for condition 1
#' @param lab_c2 label for condition 2
#'
#' @return ggplot object
#' @importFrom tidyr pivot_longer
getBack2BackBarplot <- function(tab_c1, tab_c2, 
                                lab_c1,
                                lab_c2){
    difference_df <- data.frame(clusters = unique(c(tab_c1$clusters,
                                                    tab_c2$clusters)))
    difference_df$tot_c1 <- tab_c1$n_paracrine + tab_c1$n_autocrine
    difference_df$tot_c2 <- tab_c2$n_paracrine + tab_c2$n_autocrine
    difference_df$diff_c1_c2 <- difference_df$tot_c1 - difference_df$tot_c2
    
    # From wide to long
    tab_c1 <- tidyr::pivot_longer(tab_c1,
                                  cols = c("n_paracrine","n_autocrine"), 
                                  names_to = "type", 
                                  names_prefix = "n_",
                                  values_to = "n_int")
    tab_c2 <- tidyr::pivot_longer(tab_c2, 
                                  cols = c("n_paracrine","n_autocrine"), 
                                  names_to = "type", 
                                  names_prefix = "n_",
                                  values_to = "n_int")
    # Add column for condition
    tab_c1$condition <- lab_c1
    tab_c2$condition <- lab_c2
    
    # Bind dfs
    barplot_df <- rbind(tab_c1, tab_c2)
    
    cluster.colors <- scales::hue_pal(c = 80, l = 80)(
        length(unique(barplot_df$clusters)))
    
    g <- ggplot() +
        geom_bar(data = subset(barplot_df, condition == lab_c1), 
                 aes(y = n_int, x = clusters, 
                     fill = type, 
                     color = clusters),
                 position="stack", stat="identity") +
        geom_bar(data = subset(barplot_df, condition == lab_c2), 
                 aes(y = -n_int, x = clusters, 
                     fill = type, 
                     color = clusters),
                 position="stack", stat="identity") +
        geom_point(data = difference_df, aes(x = clusters, y = diff_c1_c2)) +
        geom_line(data = difference_df, aes(x = clusters, y = diff_c1_c2, group = 1), 
                  stat = "identity") +
        theme_minimal() +
        scale_fill_manual(values = c("#606060", "#C0C0C0")) + 
        scale_color_manual(values = cluster.colors) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
              text = element_text(size=20),
              strip.text = element_blank()) +
        guides(color = FALSE) + 
        labs(x = "Clusters", y = "# Interactions", 
             title = "Total number of interactions per cluster",
             subtitle = paste0(lab_c1, " (positive y-axis) vs ", lab_c2,
                               "  (negative y-axis)"))
    
    return(g)
    
}



#' Get radar plot of relative numbers of interactions for a certain cell type
#'
#' @param tab_c1 table from csv file (barplot#2) containing data for condition 1
#' @param tab_c2 table from csv file (barplot#2) containing data for condition 2
#' @param lab_c1 label for condition 1
#' @param lab_c2 label for condition 2
#' @param cell_name label of cell type of interest
#'
#' @return plot
#' @importFrom fmsb radarchart
#' @importFrom data.table transpose
getRadarPlot <- function(tab_c1, tab_c2, lab_c1,
                         lab_c2, cell_name){
    df <- merge(tab_c1, tab_c2, by = "Clusters", all = TRUE)
    colnames(df) <- c("Clusters", "nint_c1", "nint_c2")
    df[is.na(df)] <- 0
    
    cluster_names <- df$Clusters
    # add max and min
    max_nint <- max(df[, -1])
    df <- add_column(df, max_nint, .after = "Clusters")
    df <- add_column(df, "min_nint" = 0, .after = "max_nint")
    
    radar_df <- data.table::transpose(df[, -1])
    
    rownames(radar_df) <- c("max", "min", lab_c1, lab_c2)
    colnames(radar_df) <- cluster_names
    
    color <- c("#438ECC", "#E97778")
    
    fmsb::radarchart(
        radar_df, axistype = 1,
        # Customize the polygon
        pcol = color, 
        pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
        # Customize the grid
        cglcol = "grey", cglty = 1, cglwd = 0.8,
        # Customize the axis
        axislabcol = "grey30", 
        # Variable labels
        vlcex = 1.2, vlabels = colnames(radar_df),
        caxislabels = round(seq(from = 0, to = radar_df["max",1], length.out = 5)), 
        title = cell_name
    )
    legend(
        x = "bottomleft", legend = rownames(radar_df[-c(1,2),]), horiz = FALSE,
        bty = "n", pch = 20 , col = color,
        text.col = "black", cex = 1, pt.cex = 1.5
    )
    
}
