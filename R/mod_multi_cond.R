#' multi_cond UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_multi_cond_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      box(title = h3("Comparison of multiple conditions"),
          width = 12,
          
          hr(),
          column(width = 4,
                 uiOutput(ns("sel_cond1_ui"))
          ),
          column(width = 4,
                 uiOutput(ns("sel_cond2_ui"))
                 ),
          column(width = 3,
                 br(),
                 br(),
                 actionButton(ns("go"), 
                              label = "Compare!", 
                              class = "btn-info"))
      ), #box
      fluidRow(
        column(width = 12,
               box(width = 3,
                   status = "primary",
                   title = "Cluster-verse based",
                   solidHeader = TRUE,
                   
                   h4("Back-to-back Barplot"),
                   downloadButton(ns("download_backbar_pdf"), 
                                  "Download Barplot (pdf)"),
                   downloadButton(ns("download_backbar_tiff"), 
                                  "Download Barplot (tiff)"),
                   hr(),
                   h4("Radar plot"),
                   uiOutput(ns("add_third_radar_ui")),
                   uiOutput(ns("radar_vp_ui")),
                   
                   downloadButton(ns("download_radar_pdf"), 
                                  "Download Radar (pdf)"),
                   downloadButton(ns("download_radar_tiff"), 
                                  "Download Radar (tiff)")
               ),
               
               tabBox(
                 id = 'clust-verse_tabbox',
                 width = 9,
                 
                 tabPanel(h4("Back-to-Back Barplot"),
                          plotOutput(ns("backbar")) ,
                          DT::DTOutput(ns("debug_table1")),
                          DT::DTOutput(ns("debug_table2"))
                 ),
                 tabPanel(h4("Radar Plot"),
                          plotOutput(ns("radar")) 
                 )
               )
        )
      ), #fluidrow
      fluidRow(
        column(width = 12,
               box(width = 3,
                   status = "success",
                   title = "Gene-verse based",
                   solidHeader = TRUE,
                   
                   uiOutput(ns("add_third_gene_ui")),
                   downloadButton(ns("download_dotplot_pdf"), 
                                  "Download Dotplot (pdf)"),
                   downloadButton(ns("download_dotplot_tiff"), 
                                  "Download Dotplot (tiff)")
                   
                   
               ),
               
               tabBox(
                 id = 'gene-verse_tabbox',
                 width = 9,
                 
                 tabPanel(h4("Table"),
                          uiOutput(ns("gene_table_ui"))
                 ),
                 tabPanel(h4("Dot Plot"),
                          uiOutput(ns("dotplot_ui"))
                 )
               )
        )
      )# fluidrow
      
    ) # fluidRow
 
  ) # tagList
}
    
#' multi_cond Server Functions
#'
#' @noRd 
mod_multi_cond_server <- function(id,
                                  db.list,
                                  filt.data.list){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
    # reverse list elements <-> names
    db.names <- as.list(names(db.list()))
    names(db.names) <- unlist(as.character(db.list()))
    
    
    output$sel_cond1_ui <- renderUI({
      selectInput(ns("sel_cond1"), 
                  label = h4("Select condition #1:"),
                  choices = db.list(),
                  multiple = FALSE)
    })
    output$sel_cond2_ui <- renderUI({
      selectInput(ns("sel_cond2"), 
                  label = h4("Select condition #2:"),
                  choices = db.list(),
                  multiple = FALSE)
    })
    
    # Adding third option for radar plot
    output$add_third_radar_ui <- renderUI({
      selectInput(ns("sel_cond3_radar"), 
                  label = "Add condition #3:",
                  choices = db.list(),
                  multiple = FALSE)
    })
    
    # Adding third option for dot plot
    output$add_third_gene_ui <- renderUI({
      selectInput(ns("sel_cond3_dot"), 
                  label = "Add condition #3:",
                  choices = db.list(),
                  multiple = FALSE)
    })
    
    
 
    
    ######------- Back to back barplot
    
    observeEvent(input$go, {
      
      # CCC data condition 1
      data_cond1 <- filt.data.list()[[isolate({input$sel_cond1})]]
      # Get barplot dataframe for condition #1
      barplotDF1 <- reactive({
        getBarplotDF(data_cond1, unlist(getClusterNames(data_cond1)))
      })

      # CCC data condition 2
      data_cond2 <- filt.data.list()[[isolate({input$sel_cond2})]]
      # Get barplot dataframe for condition #2
      barplotDF2 <- reactive({
        getBarplotDF(data_cond2, unlist(getClusterNames(data_cond2)))
      })

      
      
      req(barplotDF1(), barplotDF2())
      
      
      # output$debug_table1 <- DT::renderDT({
      #   barplot_df
      # }, filter = list(position = 'top', clear = FALSE),
      # options = list(scrollX= TRUE,
      #                scrollCollapse = TRUE,
      #                processing = FALSE))
      
      # output$debug_table2 <- DT::renderDT({
      #   barplotDF2()
      # }, filter = list(position = 'top', clear = FALSE),
      # options = list(scrollX= TRUE,
      #                scrollCollapse = TRUE,
      #                processing = FALSE))
      # 
      
      
      
      
      
      b2b_barplot <- getBack2BackBarplot(tab_c1 = barplotDF1(),
                                         tab_c2 = barplotDF2(),
                                         lab_c1 = db.names[[isolate({input$sel_cond1})]],
                                         lab_c2 = db.names[[isolate({input$sel_cond2})]])


      output$backbar <- renderPlot({
        b2b_barplot
      })

      # Download BackBar (tiff)
      output$download_backbar_tiff <- downloadHandler(
        filename = function() {
          paste0("MC_",db.names[[input$sel_cond1]], "VS", db.names[[input$sel_cond2]], "_back2back_barplot.tiff")
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
          paste0("MC_",db.names[[input$sel_cond1]], "VS", db.names[[input$sel_cond2]], "_back2back_barplot.pdf")
        },
        content = function(file) {

          pdf(file)
          plot(b2b_barplot)
          dev.off()
        }
      )

      
      
      # Selection of viewpoint for radar plot
      # get cluster names
      cluster.list <- as.list(intersect(unlist(getClusterNames(data_cond1)), unlist(getClusterNames(data_cond2))))
      output$radar_vp_ui <- renderUI({
        selectInput(ns("sel_vp_radar"),
                    label = "Select Viewpoint cluster:",
                    choices = cluster.list,
                    multiple = FALSE)
      })

      
      ##########------------ Radar plot
      observeEvent(input$sel_vp_radar, {
        
        # CCC data condition 1
        data_cond1 <- filt.data.list()[[isolate({input$sel_cond1})]]
        
        # Get barplot dataframe for condition #1 - from VP
        barplotDF_VP_1 <- reactive({
          getBarplotDF2(data_cond1, 
                        unlist(getClusterNames(data_cond1)),
                        input$sel_vp_radar)
        })
        
        # CCC data condition 2
        data_cond2 <- filt.data.list()[[isolate({input$sel_cond2})]]
        
        # Get barplot dataframe for condition #2 - from VP
        barplotDF_VP_2 <- reactive({
          getBarplotDF2(data_cond2, 
                        unlist(getClusterNames(data_cond2)),
                        input$sel_vp_radar)
        })
        
        req(barplotDF_VP_1(), barplotDF_VP_2())
        
        
        
        output$radar <- renderPlot({
          getRadarPlot(tab_c1 = barplotDF_VP_1(),
                       tab_c2 = barplotDF_VP_2(),
                       lab_c1 = db.names[[isolate({input$sel_cond1})]],
                       lab_c2 = db.names[[isolate({input$sel_cond2})]],
                       cell_name = input$sel_vp_radar)
        })
        
        
        # Download Radar (tiff)
        output$download_radar_tiff <- downloadHandler(
          filename = function() {
            paste0("MC_", db.names[[input$sel_cond1]], "VS", db.names[[input$sel_cond2]],
                   "_radar_", input$sel_vp_radar , ".tiff")
          },
          content = function(file) {

            tiff(file, width = 600)
            getRadarPlot(tab_c1 = barplotDF_VP_1(),
                         tab_c2 = barplotDF_VP_2(),
                         lab_c1 = db.names[[input$sel_cond1]],
                         lab_c2 = db.names[[input$sel_cond2]],
                         cell_name = input$sel_vp_radar)
            dev.off()
          }
        )

        # Download Radar (pdf)
        output$download_radar_pdf <- downloadHandler(
          filename = function() {
            paste0("MC_", db.names[[input$sel_cond1]], "VS", db.names[[input$sel_cond2]],
                   "_radar_", input$sel_vp_radar , ".pdf")
          },
          content = function(file) {

            pdf(file)
            getRadarPlot(tab_c1 = barplotDF_VP_1(),
                         tab_c2 = barplotDF_VP_2(),
                         lab_c1 = db.names[[input$sel_cond1]],
                         lab_c2 = db.names[[input$sel_cond2]],
                         cell_name = input$sel_vp_radar)
            dev.off()
          }
        )
      })
      
      
      
      
      
      
      
      
      
    })
    
    
    
    
    
    
    
    
    
  })
}
    

