#' cluster_verse UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList fluidRow column checkboxGroupInput
#' sliderInput numericInput sidebarLayout sidebarPanel downloadButton 
#' selectInput radioButtons uiOutput fileInput textInput
#' @importFrom shinydashboard valueBoxOutput
#' @importFrom DT DTOutput
#' @importFrom visNetwork visNetworkOutput
#' @importFrom plotly plotlyOutput
#' @importFrom shinycssloaders withSpinner
mod_cluster_verse_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(width = 8,
             h2("Cluster-verse")),
      valueBoxOutput(ns("tot_inter"))
    ),
    fluidRow(
      box(title = "Filters",
          width = 12,
          status = "primary",
          solidHeader = TRUE,
          collapsible = TRUE,
          fluidRow(
            column(width = 4,
                   uiOutput(ns("clust_checkbox_ui"))
            ),
            column(width = 4,
                   uiOutput(ns("minScore_slider_ui"))
            ),
            column(width = 4,
                   uiOutput(ns("maxPval_slider_ui"))
            )
          )
          
  
      )
    ),
    fluidRow(
      tabBox(
        id = ns("cluster_verse_tabbox"),
        width = 12,
        height = "auto",
        
        tabPanel(h4("Network"),
                 sidebarLayout(
                   sidebarPanel(width = 3,
                                uiOutput(ns("clust_net_select_ui"))
                                ,
                                hr(),
                                radioButtons(ns("num_or_weight_radio"),
                                             label = "Show",
                                             choices = list("Number of interactions" = "n_int",
                                                            "Weighted number of interactions (by score)" = "weighted"),
                                             ),
                                radioButtons(ns("edge_weight"),
                                             label = "Scale edges weight",
                                             choices = list("small", "medium", "large"),
                                ),
                                hr(),
                                checkboxGroupInput(ns("autocrine_checkbox_net"), 
                                                   label = "Interaction Type",
                                                   choices = list("Autocrine", 
                                                                  "Paracrine"), 
                                                   inline = TRUE,
                                                   selected = c("Autocrine", 
                                                                "Paracrine")),
                                hr(),
                                actionButton(ns("download_network"), 
                                               "Network (html)",
                                             icon = icon("download"))
                                
                   ),
                   mainPanel(width = 9,
                             visNetworkOutput(ns("cluster.net"), 
                                              height = "550px") 
                   )
                 )
                 
                 
        ),
        tabPanel(h4("Barplot"),
                 h4("Overall interactions"),
                 sidebarLayout(
                   sidebarPanel(width = 3,
                                radioButtons(ns("num_or_weight_bar1"),
                                             label = "Show",
                                             choices = list("Number of interactions" = "n_int",
                                                            "Weighted number of interactions (by score)" = "weighted"),
                                ),
                                checkboxGroupInput(ns("autocrine_checkbox_bar"), 
                                                   label = "Interaction Type",
                                                   choices = list("Autocrine", 
                                                                  "Paracrine"), 
                                                   inline = TRUE,
                                                   selected = c("Autocrine", 
                                                                "Paracrine")),
                                hr(),
                                actionButton(ns("download_barClust_html"), 
                                               "Barplot (html)",
                                             icon = icon("download")),
                                actionButton(ns("download_barClust_tiff"), 
                                               "Barplot (tiff)",
                                             icon = icon("download")),
                                actionButton(ns("download_barClust_table"), 
                                               "Barplot (csv)",
                                             icon = icon("download"))
                                
                                
                   ),
                   mainPanel(width = 9,
                             plotlyOutput(ns("cluster.bar")) %>% withSpinner()
                   )
                 ),
                 br(),
                 br(),
                 h4("Overall interactions from Viewpoint cluster"),
                 sidebarLayout(
                   sidebarPanel(width = 3,
                                selectInput(ns("clust_barplot2"), 
                                            label = "Select Cluster",
                                            choices = list("clust1"),
                                            selected = NULL
                                            ),
                                hr(),
                                actionButton(ns("download_barClust2_html"), 
                                               "Barplot (html)",
                                             icon = icon("download")),
                                actionButton(ns("download_barClust2_tiff"), 
                                               "Barplot (tiff)",
                                             icon = icon("download")),
                                actionButton(ns("download_barClust2_table"), 
                                               "Barplot (csv)",
                                             icon = icon("download"))
                                
                                
                   ),
                   mainPanel(width = 9,
                             plotlyOutput(ns("cluster.bar2")) %>% withSpinner()
                   )
                 )
                 
        ),
        
        tabPanel(h4("Table"),
                 sidebarLayout(
                   sidebarPanel(width = 4,
                                selectInput(ns("vp_table"), 
                                            label = "Select Viewpoint:",
                                            choices = list("cl1", "cl2"), 
                                            multiple = FALSE),
                                radioButtons(ns("flow_table"), 
                                             label = "Select Flow:",
                                             choices = list("Directed, outgoing" = "directed_out",
                                                            "Directed, incoming" = "directed_in",
                                                            "Undirected" = "undirected")),
                                hr(),
                                actionButton(ns("download_table_cluster"),
                                               "Table",
                                             icon = icon("download"))
                   ),
                   mainPanel(width = 8, 
                             DT::DTOutput(ns("cluster_table"), 
                                                 height = "240px") %>% 
                               withSpinner()
                   )
                 )
        )
      ))
 
  )
}
    
#' cluster_verse Server Functions
#' @importFrom shiny renderUI numericInput updateCheckboxGroupInput
#' updateSliderInput reactiveValues downloadHandler
#' @importFrom shinydashboard renderValueBox valueBox
#' @importFrom visNetwork renderVisNetwork visNetwork visNodes visIgraphLayout
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly renderPlotly
#' @importFrom utils write.csv
#' @importFrom DT renderDT
#' @noRd 
mod_cluster_verse_server <- function(id, 
                                     input_sidebarmenu, 
                                     input.data, 
                                     cluster.list, 
                                     checkbox_selected, 
                                     minScore, 
                                     maxPval, 
                                     out_folder){
  moduleServer( id, function(input, output, session){
    
    rv <- reactiveValues(filt.data = NULL,
                         checkbox_selected_out = NULL,
                         minScore_out = NULL,
                         maxPval_out = NULL,
                         okay_flag = FALSE)
    
    observeEvent(input_sidebarmenu(), {
      if(input_sidebarmenu() == "cluster-verse"){
       out <- tryCatch({
         req(input.data())
         rv$okay_flag <- TRUE
       },
       error = function(cond){
         message("Error! Please upload you data")
       },
       warning = function(cond){
         message("war")
       })
      }
    })


    
      
    observeEvent(rv$okay_flag, {
      output$clust_checkbox_ui <- renderUI({
        req(rv$okay_flag)
        tagList(
          checkboxGroupInput(session$ns("cluster_selected_checkbox"),
                             label = h4("Cluster Selection"),
                             choices = cluster.list(),
                             selected = checkbox_selected(),
                             inline = TRUE)
        )
        
      })
      
      


      # Generate score slider

      output$minScore_slider_ui <- renderUI({
        req(rv$okay_flag, input.data())
        sliderInput(session$ns("minScore_slider"),
                    label = h4("Minimum Interaction Score"),
                    value = minScore(),
                    min = min(input.data()$score),
                    max = max(input.data()$score),
                    step = 0.01)
      })

      if(rv$okay_flag){
        if("p_value" %in% colnames(input.data())){

          output$maxPval_slider_ui <- renderUI({
            req(rv$okay_flag)
            numericInput(session$ns("maxPval_slider"),
                         label = h4("Maximum Interaction p value"),
                         value = maxPval(),
                         min = 0, max = 1, step = 0.01)
          })
        }
      }
      
      
      
      
    })
      

  observeEvent(rv$filt.data, {
    output$tot_inter <- renderValueBox({
      req(rv$filt.data)
      valueBox(nrow(rv$filt.data), h4("Num total interactions"),
               icon = icon("list"), color= "blue")
    })
  })

    



    observeEvent(c(input$cluster_selected_checkbox,
                   input$minScore_slider), {
                     req(input.data())
                     # checkboxgroup clusters
                     rv$filt.data <-  input.data() %>%
                       filter(clustA %in% input$cluster_selected_checkbox &
                                clustB %in% input$cluster_selected_checkbox) %>%
                       # slider on minimum score
                       filter(score >= input$minScore_slider)

                     rv$checkbox_selected_out <- input$cluster_selected_checkbox
                     rv$minScore_out <- input$minScore_slider

    })

    # filter on max p value separated cause not always present
    observeEvent(input$maxPval_slider, {
      req(rv$filt.data)
      rv$filt.data <- rv$filt.data[rv$filt.data$p_value <= input$maxPval_slider,]
      rv$maxPval_out <- input$maxPval_slider

    })




    # visual filter types for network
    observeEvent(input$cluster_selected_checkbox, {
      clusters.selected <- as.list(input$cluster_selected_checkbox)
      names(clusters.selected) <- input$cluster_selected_checkbox
      output$clust_net_select_ui <- renderUI({
        selectInput(session$ns("clust_net_select"),
                    label = "Select Viewpoint",
                    choices = c(list("-" = "all"), clusters.selected),
                    multiple = FALSE)
      })

    })

    data.filt.net <- reactive({
      req(rv$filt.data, input$clust_net_select)
      d <- rv$filt.data %>%
        filter(int.type %in% tolower(input$autocrine_checkbox_net))
      if(input$clust_net_select == "all"){
        d
      } else {
        d_oneclust <- d %>%
          filter(clustA == input$clust_net_select | clustB == input$clust_net_select)
      }
    })

    net <- reactive({
      req(data.filt.net())
      createNetwork(data.filt.net(), input$num_or_weight_radio,
                                     input$edge_weight)})

    # Plot network
    output$cluster.net <- renderVisNetwork({
      validate(
        need(!is.null(input$autocrine_checkbox_net), 'Check at least one interaction type!')
      )

      req(data.filt.net())
      if(any("circle" %in% net()$nodes$shape)){
        # cluster names are numbers -> no background
        visNetwork(net()$nodes, net()$edges, width = "100%") %>%
          visNodes(font = list(size = 18),
                   scaling = list(min = 10, max = 40)) %>%
          visIgraphLayout(smooth = TRUE)
      } else {
        visNetwork(net()$nodes, net()$edges, width = "100%") %>%
          visNodes(font = list(size = 18, background = "#ffffff"),
                   scaling = list(min = 10, max = 40)) %>%
          visIgraphLayout(smooth = TRUE)
      }

    })

    # download network
    observeEvent(input$download_network, {
      dir.create(file.path(out_folder(), "cluster_verse"), showWarnings = FALSE)
      file <- file.path(out_folder(), "cluster_verse", 
                        paste(input$clust_net_select, input$num_or_weight_radio,
                              input$edge_weight,  "network.html", sep = "_"))
      network <- visNetwork(net()$nodes, net()$edges, width = "100%") %>%
        visNodes(font = list(size = 18, background = "#ffffff"),
                 scaling = list(min = 10, max = 40)) %>%
        visIgraphLayout(smooth = TRUE)
      htmlwidgets::saveWidget(network, file = file, selfcontained = TRUE)
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
    })
    


    ####---Barplot---####

    # visual filter on interaction type for barplot
    data.filt.bar <- reactive({
      d <- rv$filt.data %>%
        filter(int.type %in% tolower(input$autocrine_checkbox_bar))
    })

    # Get barplot dataframe
    barplotDF <- reactive({
      getBarplotDF(data.filt.bar(), input$cluster_selected_checkbox,
                   input$num_or_weight_bar1)
    })

    # Plot barplot
    output$cluster.bar <- renderPlotly({
      createBarPlot_CV(barplotDF(), input$cluster_selected_checkbox,
                       input$num_or_weight_bar1)
      })

    # Download Barplot (html)
    observeEvent(input$download_barClust_html, {
      dir.create(file.path(out_folder(), "cluster_verse"), showWarnings = FALSE)
      file <- file.path(out_folder(), "cluster_verse", 
                        paste("all", input$num_or_weight_bar1,  "barplot.html", sep = "_"))
      fig <- createBarPlot_CV(barplotDF(),
                              input$cluster_selected_checkbox,
                              input$num_or_weight_bar1)
      htmlwidgets::saveWidget(fig, file = file, selfcontained = TRUE)
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
    })
    
    
    # Download Barplot (tiff)
    observeEvent(input$download_barClust_tiff, {
      dir.create(file.path(out_folder(), "cluster_verse"), showWarnings = FALSE)
      file <- file.path(out_folder(), "cluster_verse", 
                        paste("all", input$num_or_weight_bar1,  "barplot.tiff", sep = "_"))
      fig <- createBarPlot1_ggplot(barplotDF(),
                                   input$cluster_selected_checkbox,
                                   input$num_or_weight_bar1)
      tiff(file, width = 700)
      plot(fig)
      dev.off()
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
    })
    
    

    # Download table
    observeEvent(input$download_barClust_table, {
      dir.create(file.path(out_folder(), "cluster_verse"), showWarnings = FALSE)
      file <- file.path(out_folder(), "cluster_verse", 
                        paste("all", input$num_or_weight_bar1,  "barplot.csv", sep = "_"))
      write.csv(barplotDF(), file, quote = TRUE, row.names = FALSE)
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
    })
    
  


    ##--- Barplot per cluster

    observeEvent(input$cluster_selected_checkbox, {
      clusters.selected <- as.list(input$cluster_selected_checkbox)
      names(clusters.selected) <- input$cluster_selected_checkbox
      updateSelectInput(session, "clust_barplot2",
                        choices = clusters.selected)
    })

    # Get barplot dataframe
    barplotDF2 <- reactive({
      getBarplotDF2(rv$filt.data,
                    input$cluster_selected_checkbox,
                    input$clust_barplot2)
    })

    # Plot barplot
    output$cluster.bar2 <- renderPlotly({
      createBarPlot2_CV(barplotDF2(),
                        input$cluster_selected_checkbox,
                        input$clust_barplot2)
    })

    
    
    
    
    
    
    # Download Barplot (html)
    observeEvent(input$download_barClust2_html, {
      dir.create(file.path(out_folder(), "cluster_verse"), showWarnings = FALSE)
      file <- file.path(out_folder(), "cluster_verse", 
                        paste(as.character(input$clust_barplot2),  "barplot.html", sep = "_"))
      fig <- createBarPlot2_CV(barplotDF2(),
                               input$cluster_selected_checkbox,
                               input$clust_barplot2)
      htmlwidgets::saveWidget(fig, file = file, selfcontained = TRUE)
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
    })
    
    # Download Barplot (tiff)
    observeEvent(input$download_barClust2_tiff, {
      dir.create(file.path(out_folder(), "cluster_verse"), showWarnings = FALSE)
      file <- file.path(out_folder(), "cluster_verse", 
                        paste(as.character(input$clust_barplot2),  "barplot.tiff", sep = "_"))
      fig <- createBarPlot2_ggplot(barplotDF2(),
                                   input$cluster_selected_checkbox,
                                   input$clust_barplot2)
      tiff(file, width = 700)
      plot(fig)
      dev.off()
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
    })
    
    
    # Download table
    observeEvent(input$download_barClust2_table, {
      dir.create(file.path(out_folder(), "cluster_verse"), showWarnings = FALSE)
      file <- file.path(out_folder(), "cluster_verse", 
                        paste(as.character(input$clust_barplot2) ,  "barplot.csv", sep = "_"))
      write.csv(barplotDF2(), file, quote = TRUE, row.names = FALSE)
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
    })
    


    


    ####---Table---####
    # Update inputs Table and plot table
    observeEvent(input$cluster_selected_checkbox, {
      clusts <- as.list(input$cluster_selected_checkbox)
      names(clusts) <- input$cluster_selected_checkbox
      updateSelectInput(session, "vp_table",
                        choices = clusts)
    })

    # Table data
    table.data <- reactive({
      req(rv$filt.data)
      # select viewpoint and flow
      getIntFlow(vp = input$vp_table, rv$filt.data, flow = input$flow_table)
    })


    output$cluster_table <- DT::renderDT({
      table.data()
      d <- table.data()
      d$clustA <- as.factor(d$clustA)
      d$clustB <- as.factor(d$clustB)
      d
    }, filter = list(position = 'top', clear = FALSE),
    options = list(scrollX= TRUE,
                        scrollCollapse = TRUE,
                        processing = FALSE))
    
    
    # Download table
    observeEvent(input$download_table_cluster, {
      dir.create(file.path(out_folder(), "cluster_verse"), showWarnings = FALSE)
      file <- file.path(out_folder(), "cluster_verse", 
                        paste(input$vp_table, input$flow_table, "table.csv", sep = "_"))
      write.csv(table.data(), file, quote = TRUE, row.names = FALSE)
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
    })

    return(rv)

  })
}
    

