#' int_pair_modules UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList selectInput plotOutput downloadButton uiOutput
#' numericInput verbatimTextOutput repeatable
#' @importFrom plotly plotlyOutput
#' @importFrom DT DTOutput
#' @importFrom shinycssloaders withSpinner
#' @importFrom wordcloud2 wordcloud2Output
mod_int_pair_modules_ui <- function(id){
  ns <- NS(id)
  n <- seq_len(10)
  names(n) <- n
  tagList(
    fluidRow(
      box(title = h3("Int-Pair Modules analysis"),
          width = 12,
          
          hr(),
          column(width = 3,
                 uiOutput(ns("ipM_vp_ui"))
                 ),
          column(width = 3,
                 uiOutput(ns("ipM_flow_ui"))
                 ),
          column(width = 3,
                 # br(),
                 # actionButton(ns("go"), 
                 #             label = "Go!", 
                 #             class = "btn-info")
                 )
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
                 # selectInput(ns("ipM_UMAPcolors"), 
                 #             label = "Color UMAP by:", 
                 #             choices = list("Int-Pair Module" = "ipM_col",
                 #                            "Uniqueness Score" = "US_col"),
                 #             multiple = FALSE
                 # ),
                 uiOutput(ns("ipM_elbow_ui")),
                 uiOutput(ns("ipM_silhouette_ui"))
                 
          ),
          column(width=7,
                 actionButton(ns("download_dendro_IPM"), "Dendrogram", icon = icon("download")),
                 uiOutput(ns("ipM_dendro_ui")),
                 br(),
                 actionButton(ns("download_umap_IPM"), "UMAP", icon = icon("download")),
                 uiOutput(ns("ipM_umap_ui"))
          ),
      ), #box 
      fluidRow(
        column(width = 12,
          box(width = 3,
              status = "warning",
              title = "Int-Pair Modules: visualization",
              solidHeader = TRUE,
              
              
              uiOutput(ns("chooseIPModuleUI")),
              selectInput(ns("link_color"), label = "Color links by:",
                          choices = list("Int-Pair Module" = "ipm",
                                         "Scaled Int Score" = "score"),
                          selected = "ipm",
                          multiple = FALSE),
              actionButton(ns("download_table_IPM"), "Table (csv)", icon = icon("download")),
              actionButton(ns("download_circle_IPM_tiff"), 
                             "Circle Plot (tiff)", icon = icon("download")),
              actionButton(ns("download_circle_IPM_pdf"), 
                             "Circle Plot (pdf)", icon = icon("download"))
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
          )
        )
      ),
      fluidRow(
        column(width = 12,
          box(width = 3,
              status = "danger",
              title = "Significant Functional Terms",
              solidHeader = TRUE,
              
              uiOutput(ns("chooseIPModuleUI_signF")),
              numericInput(ns("maxPval"),
                           label = "Maximum significant p value",
                           value = 0.05,
                           min = 0, max = 1, step = 0.01),
              
              actionButton(ns("download_table_signF"), "Table (csv)", icon = icon("download"))
              
          ),
          tabBox(
            id = 'signFunc_tabbox',
            width = 9,
            
            tabPanel(h4("Table"),
                     DT::DTOutput(ns("signF_table")) 
            )
          )
        )
        
      ),
      
      
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
#' @importFrom grDevices tiff pdf dev.off
#' @importFrom stats as.dendrogram na.omit
#' @noRd 
mod_int_pair_modules_server <- function(id, 
                                        seed,
                                        input_sidebarmenu,
                                        filt.data, 
                                        genePairs_func_mat, 
                                        rank.terms,
                                        ipM_vp_selected,
                                        ipM_flow_selected,
                                        out_folder){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    n <- seq_len(10)
    names(n) <- n
    rv <- reactiveValues(okay_flag = FALSE,
                         data.vp.flow = NULL,
                         flag_nModules = 0,
                         subGenePairs_func_mat = NULL,
                         elbow_x = NULL,
                         intPairs.dendro = NULL,
                         ipM.umap = NULL,
                         ipM_vp_selected = NULL,
                         ipM_flow_selected = NULL)
    
    ##--- Alert module if functional annotation has not been run
    observeEvent(input_sidebarmenu(), {
      if(input_sidebarmenu() == "ipModules"){
        out <- tryCatch({
          req(genePairs_func_mat())
          rv$okay_flag <- TRUE
        },
        error = function(cond){
          shinyalert(text = "Please perform the functional annotation in 
                 Function-verse before proceeding with the analysis!",
                     type = "warning",
                     showCancelButton = FALSE)
        },
        warning = function(cond){
          shinyalert(text = "Please perform the functional annotation in 
                 Function-verse before proceeding with the analysis!",
                     type = "warning",
                     showCancelButton = FALSE)
        })
        
      }
      
      
    })
    
    
    
    
    ##------- Int-Pair Modules analysis
    observeEvent({
      req(rv$okay_flag)
      filt.data()}, {
      
      clusts <- getClusterNames(filt.data())
      
      
      output$ipM_vp_ui <- renderUI({
        selectInput(session$ns("ipM_vp"), 
                    label = "Select Viewpoint:", 
                    choices = clusts,
                    selected = ipM_vp_selected(),
                    multiple = FALSE
        )
        
      })
      output$ipM_flow_ui <- renderUI({
      selectInput(session$ns("ipM_flow"), 
                  label = "Select Flow:", 
                  choices = list("Directed, outgoing" = "directed_out",
                                 "Directed, incoming" = "directed_in",
                                 "Undirected" = "undirected"),
                  selected = ipM_flow_selected(),
                  multiple = FALSE
      )
      })
      
    })
    
    observeEvent(input$ipM_vp, {
      rv$ipM_vp_selected <- input$ipM_vp
    })
    
    observeEvent(input$ipM_flow, {
      rv$ipM_flow_selected <- input$ipM_flow
    })
    
    
    
    #  Subset data to Viewpoint and flow
    
    observeEvent({
      req(rv$okay_flag, input$ipM_vp)
          c(input$ipM_vp, input$ipM_flow)}, {
      req(filt.data())
            
      freezeReactiveValue(input, "ipM_Nmodules")
      rv$data.vp.flow <- getIntFlow(vp = input$ipM_vp,
                                    input.data = filt.data(), 
                                    flow = input$ipM_flow)
      # Check subset data  
      req(genePairs_func_mat())
      uniq.int.pair <- unique(rv$data.vp.flow$int_pair)
      # Case okay: enough int-pairs to create modules
      if(length(uniq.int.pair) > 10){
        rv$subGenePairs_func_mat <- subsetFuncMatBYFlow(genePairs_func_mat(), 
                                                        rv$data.vp.flow)
        rv$flag_Nmodules <- 0
      } # Case < 10: only one module
      else if(length(uniq.int.pair) <= 10 & length(uniq.int.pair) > 1){
        shinyalert(text = "Less than 10 unique int-pairs are enriched for your 
                   choice of viewpoint and flow: defining only 1 module.",
                   type = "warning",
                   showCancelButton = FALSE)
        rv$subGenePairs_func_mat <- subsetFuncMatBYFlow(genePairs_func_mat(), 
                                                        rv$data.vp.flow)
        rv$elbow_x <- NULL
        rv$flag_Nmodules <- 1
        updateSelectInput(session, "ipM_Nmodules",
                          selected = 1,
                          choices = list("1" = 1))
      } # Case 1: only one unique int-pair or Case 0: no data in subset
      else if(length(uniq.int.pair) == 1 | nrow(rv$data.vp.flow) == 0){
        rv$subGenePairs_func_mat <- NULL
        rv$elbow_x <- NULL
        rv$flag_Nmodules <- -1
        output$ipM_silhouette <- renderPlot({ NULL })
        output$ipM_elbow <- renderPlot({ NULL })
        output$ipM_dendro <- renderPlot({ NULL })
        output$ipM_umap <- renderPlotly({ NULL })
        output$IPM_circle_ui <- renderUI({ NULL })
        
        # Show message
        if(length(uniq.int.pair) == 1){
          shinyalert(text = "Only one int-pair is enriched to the 
                   viewpoint and flow chosen. Definition of Modules 
                     is not possible.",
                     type = "warning",
                     showCancelButton = FALSE)
        } else {
          shinyalert(text = "Sorry, there in no int-pair corresponding to the 
                   viewpoint and flow chosen. Definition of Modules 
                     is not possible.",
                     type = "warning",
                     showCancelButton = FALSE)
        }
      }
        
        
      
      
      if(is.null(rv$subGenePairs_func_mat)){
        rv$intPairs.dendro <- NULL
        
      } else {
        # Generate UIs for plotting
        output$ipM_elbow_ui <- renderUI({
          plotOutput(ns("ipM_elbow")) %>% withSpinner()
        })
        
        output$ipM_silhouette_ui <- renderUI({
          plotOutput(ns("ipM_silhouette")) %>% withSpinner()
        })
        
        output$ipM_dendro_ui <- renderUI({
          plotOutput(ns("ipM_dendro")) %>% withSpinner()
        })
        
        output$ipM_umap_ui <- renderUI({
          plotlyOutput(ns("ipM_umap")) %>% withSpinner()
        })
        
        
        
        # Build dendrogram int-pair modules
        if(!is.null(seed())){
          rep_dendroIntPairModules <- shiny::repeatable(dendroIntPairModules, 
                                                        seed = seed())
          rv$intPairs.dendro <- rep_dendroIntPairModules(rv$subGenePairs_func_mat)
        } else {
          rep_dendroIntPairModules <- shiny::repeatable(dendroIntPairModules)
          rv$intPairs.dendro <- rep_dendroIntPairModules(rv$subGenePairs_func_mat)
        }
        
        # Predict best number of clusters by elbow method
        if(rv$flag_Nmodules == 1){
          output$ipM_silhouette <- renderPlot({ NULL })
          output$ipM_elbow <- renderPlot({ NULL })
        } else {
          elbow_plot <- factoextra::fviz_nbclust(
            rv$intPairs.dendro$umap[, c("UMAP_1", "UMAP_2")],
            factoextra::hcut, method = "wss",
            k.max = ifelse(nrow(rv$intPairs.dendro$umap) > 10, 
                           10,
                           nrow(rv$intPairs.dendro$umap) - 1))
          # Get number
          rv$elbow_x <- round(elbowPoint(as.numeric(elbow_plot$data$clusters), 
                                         elbow_plot$data$y)$x)
          # Plot elbow plot
          output$ipM_elbow <- renderPlot({
            elbow_plot + 
              geom_vline(xintercept = rv$elbow_x, linetype = 2, 
                         colour = "steelblue") +
              ggtitle("Optimal number of Modules: Elbow plot") +
              xlab("Number of Modules")
          })
          # Predict best number of clusters by average silhouette
          output$ipM_silhouette <- renderPlot({
            factoextra::fviz_nbclust(
              rv$intPairs.dendro$umap[, c("UMAP_1", "UMAP_2")],
              factoextra::hcut, method = "silhouette",
              k.max = ifelse(nrow(rv$intPairs.dendro$umap) > 10, 
                             10,
                             nrow(rv$intPairs.dendro$umap) - 1)) +
              ggtitle("Optimal number of Modules: average silhouette") +
              xlab("Number of Modules")
            
          })
        }
   
      }
  
    })
    

    # update selectInput with predicted number based on elbow method
    observeEvent(rv$elbow_x, {
      updateSelectInput(session, "ipM_Nmodules",
                        selected = rv$elbow_x,
                        choices = as.list(n))
    })
    
    gpModules_assign <- reactive({
      req(rv$intPairs.dendro)
      dendextend::cutree(rv$intPairs.dendro$h_clust,
                         k = ifelse(rv$flag_Nmodules == 1, 1, input$ipM_Nmodules),
                         order_clusters_as_data = FALSE)
    })
    
    dendro <- reactive({
      req(rv$intPairs.dendro)
      ## dendrogram of int-pair modules
      d <- as.dendrogram(rv$intPairs.dendro$h_clust)
      dendro <- dendextend::color_branches(d,
                                           k = ifelse(rv$flag_Nmodules == 1, 1, input$ipM_Nmodules),
                                           groupLabels = TRUE) %>%
        dendextend::set("labels",
                        rep("", times =
                              length(rv$intPairs.dendro$h_clust$labels)))
      return(dendro)
    })
    
    
    output$ipM_dendro <- renderPlot({
      plot(dendro(),
           horiz = TRUE,
           main = "Dendrogram of Int-pairs",
           #family = "sans",
           cex.main = 1.3)
      
    })
    
    # Download Dendro (tiff)
    observeEvent(input$download_dendro_IPM, {
      dir.create(file.path(out_folder(), "ipModules"), showWarnings = FALSE)
      file <- file.path(out_folder(), "ipModules", 
                        paste(input$ipM_vp, input$ipM_flow, input$ipM_Nmodules, 
                              "dendrogram.tiff", sep = "_"))
      tiff(file)
      plot(dendro(),
           horiz = TRUE,
           main = "Dendrogram of Int-pairs",
           #family = "sans",
           cex.main = 1.3)
      dev.off()
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
    })
    
    
    
    
    ####--- UMAP ---####
    ipM.umap <- reactive({
      req(rv$intPairs.dendro, gpModules_assign())
      
      fig <- getUMAPipModules(rv$intPairs.dendro,
                       gpModules_assign(),
                       ipm_colors = colorspace::rainbow_hcl(as.numeric(input$ipM_Nmodules)))
      return(fig)
    })
      
    output$ipM_umap <- renderPlotly({
      req(ipM.umap())
      ipM.umap()
    })
    
    # Download UMAP
    observeEvent(input$download_umap_IPM, {
      dir.create(file.path(out_folder(), "ipModules"), showWarnings = FALSE)
      file <- file.path(out_folder(), "ipModules", 
                        paste(input$ipM_vp, input$ipM_flow, input$ipM_Nmodules, 
                              "UMAP.html", sep = "_"))
      htmlwidgets::saveWidget(ipM.umap(), file = file, selfcontained = TRUE)
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
    })
      
   
    
    ####--- visualization
    output$chooseIPModuleUI <- renderUI({
      req(gpModules_assign())
      selectInput(ns("chooseIPModule"),
                  label = "Choose Int-Pair Module:",
                  choices = as.list(unique(gpModules_assign())),
                  multiple = FALSE
      )
    })
    
    selected.data <- reactive({
      req(input$chooseIPModule, gpModules_assign())
      rv$data.vp.flow %>%
        filter(int_pair %in% names(gpModules_assign())[
          gpModules_assign() == as.numeric(input$chooseIPModule)])
    })
    
    ## Plot table selected int-pair module
    output$IPM_table <- DT::renderDT({
      req(selected.data())
      
      d <- selected.data()
      d$clustA <- as.factor(d$clustA)
      d$clustB <- as.factor(d$clustB)
      d
    }, filter = list(position = 'top', clear = FALSE),
    options = list(scrollX= TRUE,
                   scrollCollapse = TRUE,
                   processing = FALSE,
                   pageLength = 5), escape = FALSE)
    
    
      
      
    ####--- Circle plot ---####
    
    output$IPM_circle_ui <- renderUI({
      req(selected.data())
      if(nrow(selected.data()) <= 200){
        plotOutput(ns("IPM_circle"),
                   height = ifelse(nrow(selected.data()) <= 200,
                                   600,
                                   round(3*nrow(selected.data())))) %>%
          withSpinner()
      } else {
        textOutput(ns("IPM_circle_error"))
      }
      
    })
    
    output$IPM_circle_error <- renderText("Error: circle plot cannot show
                                          more than 200 interactions. Try
                                          choosing a higher number of int-pair
                                          Modules!")
    
    output$IPM_circle <- renderPlot({
      req(selected.data())
      cluster.list <- getClusterNames(isolate({filt.data()}))
      # assign a color to each cluster
      cluster.colors <- hue_pal(c = 80, l = 80)(length(names(cluster.list)))
      names(cluster.colors) <- names(cluster.list)
      # Colors for modules
      ipm_colors <- colorspace::rainbow_hcl(as.numeric(input$ipM_Nmodules))
      circlePlot(selected.data(),
                 cluster_colors = cluster.colors,
                 ipm_color = ipm_colors[as.numeric(input$chooseIPModule)],
                 int_flow = isolate({input$ipM_flow}),
                 link.color = input$link_color)
      
    })
      
    # Download Table
    observeEvent(input$download_table_IPM, {
      dir.create(file.path(out_folder(), "ipModules"), showWarnings = FALSE)
      file <- file.path(out_folder(), "ipModules", 
                        paste("Module", input$chooseIPModule, input$ipM_vp, input$ipM_flow, input$ipM_Nmodules, 
                              "table.csv", sep = "_"))
      write.csv(selected.data(), file, quote = TRUE, row.names = FALSE)
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
    })
    
    
    # Download Circle tiff
    observeEvent(input$download_circle_IPM_tiff, {
      dir.create(file.path(out_folder(), "ipModules"), showWarnings = FALSE)
      file <- file.path(out_folder(), "ipModules", 
                        paste("Module", input$chooseIPModule, input$ipM_vp, input$ipM_flow, input$ipM_Nmodules, 
                              "circleplot.tiff", sep = "_"))
      cluster.list <- getClusterNames(isolate({filt.data()}))
      # assign a color to each cluster
      cluster.colors <- hue_pal(c = 80, l = 80)(length(names(cluster.list)))
      names(cluster.colors) <- names(cluster.list)
      # Colors for modules
      ipm_colors <- colorspace::rainbow_hcl(as.numeric(input$ipM_Nmodules))
      tiff(file, height = 600, width = 700)
      circlePlot(selected.data(),
                 cluster_colors = cluster.colors,
                 ipm_color = ipm_colors[as.numeric(input$chooseIPModule)],
                 int_flow = isolate({input$ipM_flow}),
                 link.color = input$link_color)
      dev.off()
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
    })
    
    # Download Circle pdf
    observeEvent(input$download_circle_IPM_pdf, {
      dir.create(file.path(out_folder(), "ipModules"), showWarnings = FALSE)
      file <- file.path(out_folder(), "ipModules", 
                        paste("Module", input$chooseIPModule, input$ipM_vp, input$ipM_flow, input$ipM_Nmodules, 
                              "circleplot.pdf", sep = "_"))
      cluster.list <- getClusterNames(isolate({filt.data()}))
      # assign a color to each cluster
      cluster.colors <- hue_pal(c = 80, l = 80)(length(names(cluster.list)))
      names(cluster.colors) <- names(cluster.list)
      # Colors for modules
      ipm_colors <- colorspace::rainbow_hcl(as.numeric(input$ipM_Nmodules))
      pdf(file, height = 8, width = 7)
      circlePlot(selected.data(),
                 cluster_colors = cluster.colors,
                 ipm_color = ipm_colors[as.numeric(input$chooseIPModule)],
                 int_flow = isolate({input$ipM_flow}),
                 link.color = input$link_color)
      dev.off()
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
    })
      
      
      
    output$chooseIPModuleUI_signF <- renderUI({
      req(gpModules_assign())
      selectInput(ns("chooseIPModule_signF"),
                  label = "Choose Int-Pair Module:",
                  choices = as.list(unique(gpModules_assign())),
                  multiple = FALSE
      )
    })
    
    
    
    
    # Permutation test to get significant functional terms for all int-pair modules
    significantFunc <- reactive({
      if(!is.null(gpModules_assign()) & rv$flag_Nmodules == 0){
        req(rv$subGenePairs_func_mat)
        getSignificantFunctions(rv$subGenePairs_func_mat,
                                gpModules_assign(),
                                rank.terms(),
                                input$maxPval)
      } else {
        NULL
      }
      
    })
      
      
      
    # get terms significant for the selected int-pair module
    sign_table <- reactive({
      req(significantFunc(), input$chooseIPModule_signF)
      if(!is.null(significantFunc()) & nrow(significantFunc()) > 0){
        significantFunc() %>%
          filter(int_pairModule == as.integer(input$chooseIPModule_signF)) %>%
          arrange(`p_value`)
      } else {
        data.table()
      }
      
    })
    
    
    
    output$signF_table <- DT::renderDT({
      sign_table()
    }, filter = list(position = 'top', clear = FALSE),
    options = list(scrollX= TRUE,
                   scrollCollapse = TRUE,
                   processing = FALSE,
                   pageLength = 5), escape = FALSE)
    
    # Download table
    observeEvent(input$download_table_signF, {
      dir.create(file.path(out_folder(), "ipModules"), showWarnings = FALSE)
      file <- file.path(out_folder(), "ipModules", 
                        paste("SignFun_for_IPM", input$chooseIPModule_signF, input$ipM_vp, input$ipM_flow, input$ipM_Nmodules, 
                              "table.csv", sep = "_"))
      write.csv(sign_table(), file, quote = TRUE, row.names = FALSE)
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
    })
      
      
      
      
     
      
      
      
   

   
    

    return(rv)

    


 
  })
}
    

