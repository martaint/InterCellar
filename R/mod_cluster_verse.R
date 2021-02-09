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
#' selectInput radioButtons uiOutput
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
                                downloadButton(ns("download_barClust"), 
                                               "Download Barplot")
                                
                                
                   ),
                   mainPanel(width = 9,
                             plotlyOutput(ns("cluster.bar")) %>% withSpinner()
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
      
      if("pvalue" %in% colnames(input.data())){
        output$maxPval_slider_ui <- renderUI(
          numericInput(session$ns("maxPval_slider"),
                       label = h4("Maximum Interaction pValue"),
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
      rv$filt.data <- rv$filt.data[rv$filt.data$pvalue <= input$maxPval_slider,]
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

    # Plot barplot
    output$cluster.bar <- renderPlotly({
      createBarPlot_CV(data.filt.bar(), input$cluster_selected_checkbox)
      })

    # Download Barplot
    output$download_barClust <- downloadHandler(
      filename = function() {"Cluster-verse_barplot.html"},
      content = function(file) {
        fig <- createBarPlot_CV(data.filt.bar(),
                                input$cluster_selected_checkbox)
        htmlwidgets::saveWidget(fig, file = file, selfcontained = TRUE)
      }
    )

    ####---Table---#### updateSelectInput
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
    

