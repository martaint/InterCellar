#' int_pair_modules UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList selectInput plotOutput downloadButton uiOutput
#' numericInput verbatimTextOutput
#' @importFrom plotly plotlyOutput
#' @importFrom DT DTOutput
#' @importFrom shinycssloaders withSpinner
#' @importFrom wordcloud2 wordcloud2Output
mod_int_pair_modules_ui <- function(id){
  ns <- NS(id)
  n <- 1:10
  names(n) <- n
  tagList(
    fluidRow(
      box(title = h3("Int-Pair Modules analysis"),
          width = 12,
          
          hr(),
          column(width = 3,
                 selectInput(ns("ipM_vp"), 
                             label = "Select Viewpoint:", 
                             choices = list("Cluster1", "Cluster2"),
                             multiple = FALSE
                 )),
          column(width = 3,
                 selectInput(ns("ipM_flow"), 
                             label = "Select Flow:", 
                             choices = list("Directed, outgoing" = "directed_out",
                                            "Directed, incoming" = "directed_in",
                                            "Undirected" = "undirected"),
                             multiple = FALSE
                 ))
      ), #box
      box(width = 12,
          status = "info",
          title = "Int-Pair Modules: definition",
          solidHeader = TRUE,
          
          column(width=5,
                 selectInput(ns("ipM_Nmodules"), 
                             label = "Choose Number of Int-Pair Modules:", 
                             choices = as.list(n),
                             multiple = FALSE
                 ),
                 selectInput(ns("ipM_UMAPcolors"), 
                             label = "Color UMAP by:", 
                             choices = list("Int-Pair Module" = "ipM_col",
                                            "Uniqueness Score" = "US_col"),
                             multiple = FALSE
                 ),
                 plotOutput(ns("ipM_elbow")) %>% withSpinner(),
                 plotOutput(ns("ipM_silhouette")) %>% withSpinner(),
          ),
          column(width=7,
                 downloadButton(ns("download_dendro_IPM"), "Dendrogram"),
                 plotOutput(ns("ipM_dendro")) %>% withSpinner(),
                 br(),
                 downloadButton(ns("download_umap_IPM"), "UMAP"),
                 plotlyOutput(ns("ipM_umap")) %>% withSpinner()
          ),
      ), #box 
      box(width = 3,
          status = "info",
          title = "Int-Pair Modules: visualization",
          solidHeader = TRUE,
          
          
          uiOutput(ns("chooseIPModuleUI")),
          downloadButton(ns("download_table_IPM"), "Download Table"),
          downloadButton(ns("download_circle_IPM"), "Download Circle Plot")
      ),
      
      tabBox(
        id = 'chooseIPM_tabbox',
        width = 9,
        
        tabPanel(h4("Circle Plot"),
                 uiOutput(ns("IPM_circle_ui"))
        ),
        tabPanel(h4("Table"),
                 DT::DTOutput(ns("IPM_table")) 
        )
      ),
      
      box(width = 3,
          status = "danger",
          title = "Significant Functional Terms",
          solidHeader = TRUE,
          
          uiOutput(ns("chooseIPModuleUI_signF")),
          numericInput(ns("maxPval"),
                       label = "Maximum significant pValue",
                       value = 0.05,
                       min = 0, max = 1, step = 0.01),
          numericInput(ns("minFreq_wordcloud"), 
                       label = "Min occurrence terms in Word Cloud:",
                       value = 5,
                       min = 1,
                       max = 50
          ),
          
          downloadButton(ns("download_table_signF"), "Download Table"),
          downloadButton(ns("download_wordcloud"), "Download Word Cloud")
      ),
      tabBox(
        id = 'signFunc_tabbox',
        width = 9,
        
        tabPanel(h4("Table"),
                 DT::DTOutput(ns("signF_table")) 
        ),
        tabPanel(h4("Word Cloud"),
                 wordcloud2Output(ns("signF_cloud"))
        )
      )
    )
 
  )
}
    
#' int_pair_modules Server Functions
#' @importFrom ggplot2 geom_vline ggtitle xlab
#' @importFrom shiny updateSelectInput renderPrint renderPlot
#' @importFrom factoextra fviz_nbclust hcut
#' @importFrom dendextend cutree
#' @importFrom plotly renderPlotly plot_ly layout config
#' @importFrom DT renderDT 
#' @importFrom scales hue_pal
#' @importFrom colorspace rainbow_hcl
#' @importFrom utils write.csv
#' @importFrom wordcloud2 renderWordcloud2
#' @importFrom dplyr filter arrange
#' @importFrom grDevices tiff dev.off
#' @importFrom stats as.dendrogram na.omit
#' @noRd 
mod_int_pair_modules_server <- function(id, 
                                        input_sidebarmenu,
                                        filt.data, 
                                        genePairs_func_mat, 
                                        gene.table, 
                                        rank.terms){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    n <- 1:10
    names(n) <- n
    rv <- reactiveValues(flag_nModules = 0,
                         subGenePairs_func_mat = NULL)
    
    ##--- Alert module if functional annotation has not been run
    observeEvent(input_sidebarmenu(), {
      if(input_sidebarmenu() == "ipModules"){
        if(is.null(genePairs_func_mat())){
        shinyalert(text = "Please perform the functional annotation in 
                 Function-verse before proceeding with the analysis!",
                   type = "warning",
                   showCancelButton = FALSE)
        }
      }
      
      
    })
    
    
    
    
    ##------- Int-Pair Modules analysis
    observeEvent(filt.data(), {
      
      clusts <- getClusterNames(filt.data())
      updateSelectInput(session, "ipM_vp",
                        choices = clusts)
    })
    
    #  Subset data to Viewpoint and flow
    data.vp.flow <- reactive({
      req(filt.data())
      getIntFlow(vp = input$ipM_vp, 
                 input.data = filt.data(), 
                 flow = input$ipM_flow)
    })
    
    
    
    # Check subset data  
    observeEvent(data.vp.flow(), {
      req(genePairs_func_mat())
      uniq.int.pair <- isolate(unique(data.vp.flow()$int_pair))
      # Case okay: enough int-pairs to create modules
      if(length(uniq.int.pair) > 10){
        rv$subGenePairs_func_mat <- subsetFuncMatBYFlow(genePairs_func_mat(), 
                                                        data.vp.flow())
        rv$flag_nModules <- 0
      } # Case < 10: only one module
      else if(length(uniq.int.pair) <= 10 & length(uniq.int.pair) > 1){
        shinyalert(text = "Less than 10 unique int-pairs are enriched for your 
                   choice of viewpoint and flow: defining only 1 module.",
                   type = "warning",
                   showCancelButton = FALSE)
        rv$subGenePairs_func_mat <- subsetFuncMatBYFlow(genePairs_func_mat(), 
                                                        data.vp.flow())
        rv$flag_nModules <- 1
        updateSelectInput(session, "ipM_Nmodules",
                          selected = 1,
                          choices = list("1" = 1))
      } # Case 1: only one unique int-pair
      else if(length(uniq.int.pair) == 1){
        shinyalert(text = "Only one int-pair is enriched to the 
                   viewpoint and flow chosen. Definition of Modules 
                     is not possible.",
                   type = "warning",
                   showCancelButton = FALSE)
        rv$subGenePairs_func_mat <- NULL
      } # Case 0: no data in subset
      else if(nrow(data.vp.flow()) == 0){
        shinyalert(text = "Sorry, there in no int-pair corresponding to the 
                   viewpoint and flow chosen. Definition of Modules 
                     is not possible.",
                   type = "warning",
                   showCancelButton = FALSE)
        rv$subGenePairs_func_mat <- NULL
      }
    })
    
    # Build dendrogram int-pair modules
    intPairs.dendro <- reactive({
      req(rv$subGenePairs_func_mat)
      dendroIntPairModules(rv$subGenePairs_func_mat, seed = 123)
      })
    # Predict best number of clusters by elbow method
    elbow_plot <- reactive({
      req(intPairs.dendro())
      if(rv$flag_nModules == 1){
        NULL
      } else{
        factoextra::fviz_nbclust(intPairs.dendro()$umap[, c("UMAP_1", "UMAP_2")], 
                                 factoextra::hcut, method = "wss",
                                 k.max = ifelse(nrow(intPairs.dendro()$umap) > 10, 
                                                10,
                                                nrow(intPairs.dendro()$umap) - 1))
      }
       
    })
    elbow_x <- reactive({
      req(elbow_plot())
      round(elbowPoint(as.numeric(elbow_plot()$data$clusters), 
                       elbow_plot()$data$y)$x)
    })
    
    # Plot elbow plot
    output$ipM_elbow <- renderPlot({
      req(elbow_x())
      isolate(elbow_plot()) + 
        geom_vline(xintercept = elbow_x(), linetype = 2, colour = "steelblue") +
        ggtitle("Optimal number of Modules: Elbow plot") +
        xlab("Number of Modules")
    })
  
    # Predict best number of clusters by average silhouette
    output$ipM_silhouette <- renderPlot({
      req(intPairs.dendro())
      if(rv$flag_nModules == 1){
        NULL
      } else{
        factoextra::fviz_nbclust(intPairs.dendro()$umap[, c("UMAP_1", "UMAP_2")],
                                 factoextra::hcut, method = "silhouette",
                                 k.max = ifelse(nrow(intPairs.dendro()$umap) > 10, 
                                                10,
                                                nrow(intPairs.dendro()$umap) - 1)) +
          ggtitle("Optimal number of Modules: average silhouette") +
          xlab("Number of Modules")
      }
    })
    
    # update selectInput with predicted number based on elbow method
    observeEvent(elbow_x(), {
      updateSelectInput(session, "ipM_Nmodules",
                        selected = elbow_x(),
                        choices = as.list(n))
    })
    
    gpModules_assign <- reactive({
      req(intPairs.dendro())
      dendextend::cutree(intPairs.dendro()$h_clust, 
             k = input$ipM_Nmodules, 
             order_clusters_as_data = FALSE)
      })
    
    ## Plot dendrogram of int-pair modules
    dendro <- eventReactive(c(input$ipM_Nmodules, rv$flag_nModules), {
      if(rv$flag_nModules == 1){
        NULL
      } else {
        req(gpModules_assign())
        d <- isolate(as.dendrogram(intPairs.dendro()$h_clust))
        d <- dendextend::color_branches(d, 
                                        k=input$ipM_Nmodules, 
                                        groupLabels = TRUE) %>%
          dendextend::set("labels", 
                          rep("", times = 
                                isolate(length(intPairs.dendro()$h_clust$labels))))
        d
      }
      
    })
    
    output$ipM_dendro <- renderPlot({
      req(dendro())
      plot(dendro(), 
           horiz = TRUE, 
           main="Dendrogram of Int-pairs", 
           family = "sans", 
           cex.main = 1.3)
      
    })
    
    output$download_dendro_IPM <- downloadHandler(
      filename = function() {
        paste0("IpModules_dendro_", 
               input$ipM_vp, "_", 
               input$ipM_flow, "_", 
               input$ipM_Nmodules, ".tiff")
      },
      content = function(file) {
        tiff(file)
        plot(dendro(), 
             horiz = TRUE, 
             main="Dendrogram of Int-pairs", 
             family.main = "sans", 
             cex.main = 1.3)
        dev.off()
      }
    )
    
    ####--- UMAP ---####
    
    
    ipm_colors <- reactive({
      colorspace::rainbow_hcl(as.numeric(input$ipM_Nmodules))
    })
    
    ipM.umap <- eventReactive(c(ipm_colors(), rv$flag_nModules), {
      if(rv$flag_nModules == 1){
        NULL
      } else{
        req(gpModules_assign())
        getUMAPipModules(isolate(intPairs.dendro()), 
                         gpModules_assign(), 
                         isolate(gene.table()),
                         ipm_colors(),
                         input$ipM_UMAPcolors)
      }
      
    })
    
    output$ipM_umap <- renderPlotly({
      ipM.umap()
    })
    
    output$download_umap_IPM <- downloadHandler(
      filename = function() {
        paste0("IpModules_umap_", 
               input$ipM_vp, "_", 
               input$ipM_flow, "_", 
               input$ipM_Nmodules, ".html")},
      content = function(file) {
        fig <- getUMAPipModules(intPairs.dendro(), 
                                gpModules_assign(), 
                                gene.table(),
                                ipm_colors(),
                                input$ipM_UMAPcolors)
        htmlwidgets::saveWidget(fig, file = file, selfcontained = TRUE)
      }
    )
    
    output$chooseIPModuleUI <- renderUI({
      selectInput(ns("chooseIPModule"),
                  label = "Choose Int-Pair Module:",
                  choices = as.list(unique(gpModules_assign())),
                  multiple = FALSE
      )
    })
    
    selected.data <- reactive({
      req(gpModules_assign())
      data.vp.flow() %>%
        filter(int_pair %in% names(gpModules_assign())[
          gpModules_assign() == as.numeric(input$chooseIPModule)])
    })
    
    ## Plot table selected int-pair module
    output$IPM_table <- DT::renderDT({
      req(selected.data())
      selected.data()
    }, options = list(scrollX= TRUE, 
                      scrollCollapse = TRUE, 
                      processing = FALSE,
                      pageLength = 5), escape = FALSE)
    
    
    ####--- Circle plot ---####
    
    output$IPM_circle_ui <- renderUI({
      req(selected.data())
      plotOutput(ns("IPM_circle"), 
                 height = ifelse(nrow(selected.data()) <= 200, 
                                 600, 
                                 round(3*nrow(selected.data())))) %>% 
        withSpinner()
    })
    
    output$IPM_circle <- renderPlot({
      req(selected.data())
      cluster.list <- getClusterNames(isolate({filt.data()}))
      # assign a color to each cluster 
      cluster.colors <- hue_pal(c = 80, l = 80)(length(names(cluster.list)))
      names(cluster.colors) <- names(cluster.list)
      circlePlot(selected.data(), 
                 cluster_colors = cluster.colors, 
                 int_flow = input$ipM_flow,
                 link.color = ipm_colors()[as.numeric(input$chooseIPModule)])
      
    })
    
    output$download_table_IPM <- downloadHandler(
      filename = function() {
        paste0("IpModule", 
               input$chooseIPModule, "_", 
               input$ipM_vp, "_", 
               input$ipM_flow, "_table.csv")
      },
      content = function(file) {
        write.csv(selected.data(), file, quote = TRUE, row.names = FALSE)
      }
    )
    
    output$download_circle_IPM <- downloadHandler(
      filename = function() {
        paste0("IpModule", 
               input$chooseIPModule, "_", 
               input$ipM_vp, "_", 
               input$ipM_flow, "_circleplot.tiff")
      },
      content = function(file) {
        cluster.list <- getClusterNames(isolate({filt.data()}))
        # assign a color to each cluster 
        cluster.colors <- hue_pal(c = 80, l = 80)(length(names(cluster.list)))
        names(cluster.colors) <- names(cluster.list)
        tiff(file, height = 600, width = 700)
        circlePlot(selected.data(), 
                   cluster_colors = cluster.colors, 
                   int_flow = input$ipM_flow,
                   link.color = ipm_colors()[as.numeric(input$chooseIPModule)])
        dev.off()
      }
    )
    
    output$chooseIPModuleUI_signF <- renderUI({
      selectInput(ns("chooseIPModule_signF"),
                  label = "Choose Int-Pair Module:",
                  choices = as.list(unique(gpModules_assign())),
                  multiple = FALSE
      )
    })
    
    
    
    
    # Permutation test to get significant functional terms for all int-pair modules
    significantFunc <- reactive({
      req(rv$subGenePairs_func_mat)
      getSignificantFunctions(rv$subGenePairs_func_mat, 
                              gpModules_assign(),
                              rank.terms(),
                              input$maxPval)
    })
    
   
    
    # get terms significant for the selected int-pair module
    sign_table <- reactive({
      req(significantFunc(), input$chooseIPModule_signF)
      significantFunc() %>%
        filter(int_pairModule == as.integer(input$chooseIPModule_signF)) %>%
        arrange(`pvalue`)
    })
    
    output$signF_table <- DT::renderDT({
      sign_table()
    }, options = list(scrollX= TRUE, 
                      scrollCollapse = TRUE, 
                      processing = FALSE,
                      pageLength = 5), escape = FALSE)
    
    output$download_table_signF <- downloadHandler(
      filename = function() {
        paste0("SignFun_for_IPM", 
               input$chooseIPModule_signF,"_",
               input$ipM_vp, "_", input$ipM_flow, "_table.csv")
      },
      content = function(file) {
        write.csv(sign_table(), file, quote = TRUE, row.names = FALSE)
      }
    )
    
    
  


    occurTab <- reactive({
      req(significantFunc())
      getOccurrenceTab4wordcloud(significantFunc(),
                                 input$chooseIPModule_signF,
                                 gpModules_assign(),
                                 rv$subGenePairs_func_mat)
    })
    
    output$signF_cloud <- renderWordcloud2({
      req(occurTab())
      plotWordCloud(occurTab(), input$minFreq_wordcloud)
    })
    
    output$download_wordcloud <- downloadHandler(
      filename = function() {
        paste0("Wordcloud_IPMod", 
               input$chooseIPModule_signF, "_signFunctions.html")
      },
      content = function(file) {
        graph <- plotWordCloud(occurTab(), input$minFreq_wordcloud)
        htmlwidgets::saveWidget(graph, file = file, selfcontained = TRUE)
        wordcloud2FixHTML(file, file)
      }
    )
    
    
        

        
        

    
    
 
  })
}
    

