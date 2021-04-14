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
          
          column(width = 4,
                 checkboxGroupInput(ns("cluster_selected_checkbox"),
                                    label = h4("Cluster Selection"),
                                    choices = list("cluster1", "cluster2"),
                                    selected = NULL,
                                    inline = TRUE)
          ),
          column(width = 4,
                 sliderInput(ns("minScore_slider"),
                             label = h4("Minimum Interaction Score"),
                             value = 0,
                             min = 0, max = 1, step = 0.1)
          ),
          column(width = 4,
                 uiOutput(ns("maxPval_slider_ui"))
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
                                checkboxGroupInput(ns("autocrine_checkbox_net"), 
                                                   label = "Interaction Type",
                                                   choices = list("Autocrine", 
                                                                  "Paracrine"), 
                                                   inline = TRUE,
                                                   selected = c("Autocrine", 
                                                                "Paracrine")),
                                hr(),
                                downloadButton(ns("download_network"), 
                                               "Download Network")
                                
                   ),
                   mainPanel(width = 9,
                             visNetworkOutput(ns("cluster.net"), 
                                              height = "550px") 
                   )
                 )
                 
                 
        ),
        tabPanel(h4("Barplot"),
                 h4("Barplot #1: Total number of interactions per cell type"),
                 sidebarLayout(
                   sidebarPanel(width = 3,
                                checkboxGroupInput(ns("autocrine_checkbox_bar"), 
                                                   label = "Interaction Type",
                                                   choices = list("Autocrine", 
                                                                  "Paracrine"), 
                                                   inline = TRUE,
                                                   selected = c("Autocrine", 
                                                                "Paracrine")),
                                hr(),
                                downloadButton(ns("download_barClust_html"), 
                                               "Download Barplot (html)"),
                                downloadButton(ns("download_barClust_tiff"), 
                                               "Download Barplot (tiff)"),
                                downloadButton(ns("download_barClust_table"), 
                                               "Download Barplot (table)")
                                
                                
                   ),
                   mainPanel(width = 9,
                             plotlyOutput(ns("cluster.bar")) %>% withSpinner()
                   )
                 ),
                 br(),
                 br(),
                 h4("Barplot #2: Relative number of interactions for a certain cell type"),
                 sidebarLayout(
                   sidebarPanel(width = 3,
                                selectInput(ns("clust_barplot2"), 
                                            label = "Select Cluster",
                                            choices = list("clust1"),
                                            selected = NULL
                                            ),
                                hr(),
                                downloadButton(ns("download_barClust2_html"), 
                                               "Download Barplot (html)"),
                                downloadButton(ns("download_barClust2_tiff"), 
                                               "Download Barplot (tiff)"),
                                downloadButton(ns("download_barClust2_table"), 
                                               "Download Barplot (table)")
                                
                                
                   ),
                   mainPanel(width = 9,
                             plotlyOutput(ns("cluster.bar2")) %>% withSpinner()
                   )
                 )
                 
        ),
        tabPanel(h4("1 vs 1"),
                 h4("Here you can compare the total number of interactions, 
                    in two conditions!"),
                 p("For each condition, please download the table output of 
                   Barplot #1 (in the previous tab) and upload it below."),
                 sidebarLayout(
                   sidebarPanel(width = 3,
                                textInput(ns("cond1_bar_lab"), "Condition #1 label"),
                                textInput(ns("cond2_bar_lab"), "Condition #2 label"),
                                hr(),
                                fileInput(ns("bar1_cond1"), 
                                          "Barplot #1, condition #1 (csv)", 
                                          multiple = FALSE, 
                                          accept = ".csv"),
                                fileInput(ns("bar1_cond2"), 
                                          "Barplot #1, condition #2 (csv)", 
                                          multiple = FALSE, 
                                          accept = ".csv"),
                                actionButton(ns("plot_backbar"), "Plot!"),
                                hr(),
                                downloadButton(ns("download_backbar_pdf"), 
                                               "Download Barplot (pdf)"),
                                downloadButton(ns("download_backbar_tiff"), 
                                               "Download Barplot (tiff)")
                                
                                
                   ),
                   mainPanel(width = 9,
                             plotOutput(ns("backbar")) 
                   )
                 ),
                 br(),
                 br(),
                 h4("Here you can compare the relative number of interactions, 
                              in two conditions!"),
                 p("For each condition, please download the table output of 
                             Barplot #2 (in the previous tab) 
                   and upload it below."),
        
                  
                 sidebarLayout(
                   sidebarPanel(width = 3,
                                textInput(ns("cond1_rad_lab"), "Condition #1 label"),
                                textInput(ns("cond2_rad_lab"), "Condition #2 label"),
                                textInput(ns("celltype_lab_rad"), "Cell type label"),
                                hr(),
                                fileInput(ns("bar2_cond1"), 
                                          "Barplot #2, condition #1 (csv)", 
                                          multiple = FALSE, 
                                          accept = ".csv"),
                                fileInput(ns("bar2_cond2"), 
                                          "Barplot #2, condition #2 (csv)", 
                                          multiple = FALSE, 
                                          accept = ".csv"),
                                actionButton(ns("plot_radar"), "Plot!"),
                                hr(),
                                downloadButton(ns("download_radar_pdf"), 
                                               "Download Radar (pdf)"),
                                downloadButton(ns("download_radar_tiff"), 
                                               "Download Radar (tiff)")
                                
                                
                   ),
                   mainPanel(width = 9,
                             plotOutput(ns("radar")) 
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
                                downloadButton(ns("download_table_cluster"),
                                               "Download Table")
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
mod_cluster_verse_server <- function(id, input.data){
  moduleServer( id, function(input, output, session){
    
    rv <- reactiveValues(filt.data = NULL)
    
    observeEvent(input.data(), {
      # Initialize filters 
      # Get cluster names
      cluster.list <- getClusterNames(input.data())
      updateCheckboxGroupInput(session, "cluster_selected_checkbox",
                               choices = cluster.list,
                               selected = names(cluster.list),
                               inline = TRUE)
      updateSliderInput(session, "minScore_slider",
                        value = min(input.data()$score),
                        min = min(input.data()$score),
                        max = max(input.data()$score))
      
      if("p_value" %in% colnames(input.data())){
        output$maxPval_slider_ui <- renderUI(
          numericInput(session$ns("maxPval_slider"),
                       label = h4("Maximum Interaction p value"),
                       value = 0.05,
                       min = 0, max = 1, step = 0.01)
        )
      }
    })
    
    output$tot_inter <- renderValueBox({
      req(rv$filt.data)
      valueBox(nrow(rv$filt.data), h4("Num total interactions"), 
               icon = icon("list"), color= "blue")
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
      
     
    })
    
    # filter on max p value separated cause not always present
    observeEvent(input$maxPval_slider, {
      req(rv$filt.data)
      rv$filt.data <- rv$filt.data[rv$filt.data$p_value <= input$maxPval_slider,]
    })
    
    
    # visual filter on interaction type for network
    
    data.filt.net <- reactive({
      d <- rv$filt.data %>%
        filter(int.type %in% tolower(input$autocrine_checkbox_net))
    })
    

    net <- reactive({
      req(data.filt.net())
      createNetwork(data.filt.net())})

    # Plot network
    output$cluster.net <- renderVisNetwork({
      req(net())
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
    output$download_network <- downloadHandler(
      filename = function() {"Cluster-verse_network.html"},
      content = function(file) {
        network <- visNetwork(net()$nodes, net()$edges, width = "100%") %>%
          visNodes(font = list(size = 18, background = "#ffffff"),
                   scaling = list(min = 10, max = 40)) %>%
          visIgraphLayout(smooth = TRUE)
        htmlwidgets::saveWidget(network, file = file, selfcontained = TRUE)
      }
    )
    
    
    ####---Barplot---#### 

    # visual filter on interaction type for barplot
    data.filt.bar <- reactive({
      d <- rv$filt.data %>%
        filter(int.type %in% tolower(input$autocrine_checkbox_bar))
    })
    
    # Get barplot dataframe
    barplotDF <- reactive({
      getBarplotDF(data.filt.bar(), input$cluster_selected_checkbox)
    })

    # Plot barplot
    output$cluster.bar <- renderPlotly({
      createBarPlot_CV(barplotDF(), input$cluster_selected_checkbox)
      })

    # Download Barplot (html)
    output$download_barClust_html <- downloadHandler(
      filename = function() {"Cluster-verse_barplot.html"},
      content = function(file) {
        fig <- createBarPlot_CV(barplotDF(),
                                input$cluster_selected_checkbox)
        htmlwidgets::saveWidget(fig, file = file, selfcontained = TRUE)
      }
    )
    # Download Barplot (tiff)
    output$download_barClust_tiff <- downloadHandler(
      filename = function() {"Cluster-verse_barplot.tiff"},
      content = function(file) {
        fig <- createBarPlot1_ggplot(barplotDF(), 
                                     input$cluster_selected_checkbox)
        tiff(file, width = 700)
        plot(fig)
        dev.off()
      }
    )
    
    # Download table
    output$download_barClust_table <- downloadHandler(
      filename = function() {"Cluster-verse_barplot.csv"},
      content = function(file) {
        write.csv(barplotDF(), file, quote = TRUE, row.names = FALSE)
      }
    )
    
    
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
    
    # Download Barplot html
    output$download_barClust2_html <- downloadHandler(
      filename = function() {
        paste0("Cluster-verse_barplot_clust", 
               as.character(input$clust_barplot2) ,".html")
        },
      content = function(file) {
        fig <- createBarPlot2_CV(barplotDF2(), 
                                 input$cluster_selected_checkbox,
                                 input$clust_barplot2)
        htmlwidgets::saveWidget(fig, file = file, selfcontained = TRUE)
      }
    )
    
    # Download Barplot (tiff)
    output$download_barClust2_tiff <- downloadHandler(
      filename = function() {
        paste0("Cluster-verse_barplot_clust", 
               as.character(input$clust_barplot2) ,".tiff")
      },
      content = function(file) {
        fig <- createBarPlot2_ggplot(barplotDF2(), 
                                     input$cluster_selected_checkbox,
                                     input$clust_barplot2)
        tiff(file, width = 700)
        plot(fig)
        dev.off()
      }
    )
    
    # Download table
    output$download_barClust2_table <- downloadHandler(
      filename = function() {
        paste0("Cluster-verse_barplot_clust", 
               as.character(input$clust_barplot2) ,".csv")
      },
      content = function(file) {
        write.csv(barplotDF2(), file, quote = TRUE, row.names = FALSE)
      }
    )
    
    
    
    
    ########---------- 1 vs 1 ------#######
    
    ### Back to back barplot
    observeEvent(input$plot_backbar, {
      file_c1b1 <- input$bar1_cond1
      tab_c1b1 <- read.csv(file_c1b1$datapath)

      file_c2b1 <- input$bar1_cond2
      tab_c2b1 <- read.csv(file_c2b1$datapath)
 
      
      
      b2b_barplot <- getBack2BackBarplot(tab_c1 = tab_c1b1, 
                                         tab_c2 = tab_c2b1, 
                                         lab_c1 = input$cond1_bar_lab,
                                         lab_c2 = input$cond2_bar_lab)
                            
   
      output$backbar <- renderPlot({
        b2b_barplot
      })
      
      # Download BackBar (tiff)
      output$download_backbar_tiff <- downloadHandler(
        filename = function() {
          paste0("Cluster-verse_1vs1_back2back_barplot.tiff")
        },
        content = function(file) {
          
          tiff(file, width = 700)
          plot(b2b_barplot)
          dev.off()
        }
      )
      # Download BackBar (pdf)
      output$download_backbar_pdf <- downloadHandler(
        filename = function() {
          paste0("Cluster-verse_1vs1_back2back_barplot.pdf")
        },
        content = function(file) {
          
          pdf(file)
          plot(b2b_barplot)
          dev.off()
        }
      )
      
    })
    
    
    #### Radar plots
    observeEvent(input$plot_radar, {
      file_c1b2 <- input$bar2_cond1
      tab_c1b2 <- read.csv(file_c1b2$datapath)
      
      file_c2b2 <- input$bar2_cond2
      tab_c2b2 <- read.csv(file_c2b2$datapath)
      
      
      
      output$radar <- renderPlot({
        getRadarPlot(tab_c1 = tab_c1b2, 
                     tab_c2 = tab_c2b2, 
                     lab_c1 = input$cond1_rad_lab,
                     lab_c2 = input$cond2_rad_lab, 
                     cell_name = input$celltype_lab_rad)
      })
      
      
      # Download Radar (tiff)
      output$download_radar_tiff <- downloadHandler(
        filename = function() {
          paste0("Cluster-verse_1vs1_radar_", input$celltype_lab_rad , ".tiff")
        },
        content = function(file) {
          
          tiff(file, width = 600)
          getRadarPlot(tab_c1 = tab_c1b2, 
                       tab_c2 = tab_c2b2, 
                       lab_c1 = input$cond1_rad_lab,
                       lab_c2 = input$cond2_rad_lab, 
                       cell_name = input$celltype_lab_rad)
          dev.off()
        }
      )
      
      # Download Radar (pdf)
      output$download_radar_pdf <- downloadHandler(
        filename = function() {
          paste0("Cluster-verse_1vs1_radar_", input$celltype_lab_rad , ".pdf")
        },
        content = function(file) {
          
          pdf(file)
          getRadarPlot(tab_c1 = tab_c1b2, 
                       tab_c2 = tab_c2b2, 
                       lab_c1 = input$cond1_rad_lab,
                       lab_c2 = input$cond2_rad_lab, 
                       cell_name = input$celltype_lab_rad)
          dev.off()
        }
      )
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
    # Download Table
    output$download_table_cluster <- downloadHandler(
      filename = function() {
        paste0("Cluster-verse_Tab_", input$vp_table, "_", 
               input$flow_table, ".csv")
      },
      content = function(file) {
        write.csv(table.data(), file, quote = TRUE, row.names = FALSE)
      }
    )
    
    return(rv)

  })
}
    

