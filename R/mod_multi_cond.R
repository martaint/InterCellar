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
          column(width = 3,
                 uiOutput(ns("sel_cond1_ui"))
          ),
          column(width = 3,
                 uiOutput(ns("sel_cond2_ui"))
                 ),
          column(width = 3,
                 uiOutput(ns("sel_cond3_ui"))
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
                          verbatimTextOutput(ns("debug_text2")),
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
                   
                   
                   downloadButton(ns("download_dotplot_pdf"), 
                                  "Download Dotplot (pdf)"),
                   downloadButton(ns("download_dotplot_tiff"), 
                                  "Download Dotplot (tiff)"),
                   br(),
                   downloadButton(ns("download_pie_pdf"), 
                                  "Download Piechart (pdf)"),
                   downloadButton(ns("download_pie_tiff"), 
                                  "Download Piechart (tiff)")
                   
                   
               ),
               
               tabBox(
                 id = 'gene-verse_tabbox',
                 width = 9,
                 height = "auto",
                 # tabPanel(h4("Barplot"),
                 #          plotOutput(ns("info_barplot"))
                 # ),
                 tabPanel(h4("Table"),
                          h4("Select int-pairs/cluster-pairs couplets from the Table to generate a DotPlot!"),
                          column(2,
                                 downloadButton(ns("download_geneTab"), "Download Table"),
                          ),
                          column(2,
                                 actionButton(ns("clear_rows"), "Clear Rows")
                          ),
                          br(),
                          br(),
                          DT::DTOutput(ns("uni_couplets_table")) %>% withSpinner()
                 ),
                 tabPanel(h4("Dot Plot"),
                          uiOutput(ns("dotplot.ui"))
                 ),
                 tabPanel(h4("Pie Chart"),
                          uiOutput(ns("piechart.ui"))
                 )
               )
        )
      ), # fluidrow
      fluidRow(
        column(width = 12,
               box(width = 3,
                   status = "warning",
                   title = "Function-verse based",
                   solidHeader = TRUE,
                   
                   
                   downloadButton(ns("download_sub_annot"),
                                  "Download Dotplot (pdf)"),
                   downloadButton(ns("download_unique_intpa"),
                                  "Download Dotplot (tiff)"),
                   # br(),
                   # downloadButton(ns("download_pie_pdf"), 
                   #                "Download Piechart (pdf)"),
                   # downloadButton(ns("download_pie_tiff"), 
                   #                "Download Piechart (tiff)")
                   
                   
               ),
               
               tabBox(
                 id = 'function-verse_tabbox',
                 width = 9,
                 height = "auto",
                 tabPanel(h4("Condition 1"),
                          column(2,
                                 downloadButton(ns("download_funcTab1"), "Download Table"),
                          ),
                          br(),
                          br(),
                          verbatimTextOutput(ns("debug1")),
                          verbatimTextOutput(ns("debug2")),
                          verbatimTextOutput(ns("debug3")),
                          DT::DTOutput(ns("funcTab1")) %>% withSpinner()
                 ),
                 tabPanel(h4("Condition 2"),
                          column(2,
                                 downloadButton(ns("download_funcTab2"), "Download Table"),
                          ),
                          br(),
                          br(),
                          DT::DTOutput(ns("funcTab2")) %>% withSpinner(),
                          DT::DTOutput(ns("funcTab3")) %>% withSpinner()
                 ),
                 tabPanel(h4("Condition 3"),
                          uiOutput(ns("fun_cond3_ui"))
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
                                  input_sidebarmenu,
                                  db.list,
                                  filt.data.list,
                                  func.annot.mat.list,
                                  ranked.terms.list){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
    rv <- reactiveValues(uni_couplets_tab = NULL)
    
    ##--- Alert module if functional annotation has not been run
    observeEvent(input_sidebarmenu(), {
      if(input_sidebarmenu() == "multiConditions" & is.null(func.annot.mat.list()[[1]])){
          shinyalert(text = "Please perform the functional annotation in 
                 Function-verse before proceeding with the analysis!",
                     type = "warning",
                     showCancelButton = FALSE)
      }
      
      
      
    })
    
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
    
    #Adding third option 
    output$sel_cond3_ui <- renderUI({
      selectInput(ns("sel_cond3"),
                  label = h4("Add condition #3:"),
                  choices = c(list("-" = "none"), db.list()),
                  multiple = FALSE)
      
    })
    
    
 
    
    ######------- Back to back barplot
    
    observeEvent(input$go, {
      
      # CCC data condition 1
      data_cond1 <- filt.data.list()[[isolate({input$sel_cond1})]]
      # Get barplot dataframe for condition #1
      barplotDF1 <- reactive({
        getBarplotDF(data_cond1, unlist(getClusterNames(data_cond1)), "n_int")
      })

      # CCC data condition 2
      data_cond2 <- filt.data.list()[[isolate({input$sel_cond2})]]
      # Get barplot dataframe for condition #2
      barplotDF2 <- reactive({
        getBarplotDF(data_cond2, unlist(getClusterNames(data_cond2)), "n_int")
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
      
      if(!(input$sel_cond3 %in% c("none", isolate({input$sel_cond1}), isolate({input$sel_cond2})))){
        # CCC data condition 3
        data_cond3 <- filt.data.list()[[input$sel_cond3]]
        # get cluster names
        cluster.list <- as.list(intersect(unlist(getClusterNames(data_cond1)), 
                                          intersect(unlist(getClusterNames(data_cond2)),
                                          unlist(getClusterNames(data_cond3)))))
      } else {
        # get cluster names
        cluster.list <- as.list(intersect(unlist(getClusterNames(data_cond1)), unlist(getClusterNames(data_cond2))))
      }
      
      
      
      output$radar_vp_ui <- renderUI({
        selectInput(ns("sel_vp_radar"),
                    label = "Select Viewpoint cluster:",
                    choices = cluster.list,
                    multiple = FALSE)
      })
      
      output$debug_text2 <- renderPrint({
        print(db.names)
      })
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
        
        
        if(!(input$sel_cond3 %in% c("none", isolate({input$sel_cond1}), isolate({input$sel_cond2})))){
          # CCC data condition 3
          data_cond3 <- filt.data.list()[[input$sel_cond3]]
          lab_c3 <- db.names[[input$sel_cond3]]
          barplotDF_VP_3 <- reactive({
            getBarplotDF2(data_cond3, 
                          unlist(getClusterNames(data_cond3)),
                          input$sel_vp_radar)
          })
        } else {
          lab_c3 <- NULL
          barplotDF_VP_3 <- reactive({NULL})
        }
        
        req(barplotDF_VP_1(), barplotDF_VP_2())
        
        
        
        output$radar <- renderPlot({
          getRadarPlot(tab_c1 = barplotDF_VP_1(),
                       tab_c2 = barplotDF_VP_2(),
                       tab_c3 = barplotDF_VP_3(),
                       lab_c1 = db.names[[isolate({input$sel_cond1})]],
                       lab_c2 = db.names[[isolate({input$sel_cond2})]],
                       lab_c3 = lab_c3,
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
                         tab_c3 = barplotDF_VP_3(),
                         lab_c1 = db.names[[isolate({input$sel_cond1})]],
                         lab_c2 = db.names[[isolate({input$sel_cond2})]],
                         lab_c3 = lab_c3,
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
                         tab_c3 = barplotDF_VP_3(),
                         lab_c1 = db.names[[isolate({input$sel_cond1})]],
                         lab_c2 = db.names[[isolate({input$sel_cond2})]],
                         lab_c3 = lab_c3,
                         cell_name = input$sel_vp_radar)
            dev.off()
          }
        )
      })
      
      
      ####---- Gene-verse based table
      # Based on the selected conditions, calculate which int-pairs/cluster-pairs couplets are unique
      # And generate a gene-table that shows only those
      
      
      observeEvent(input$go, {
        # CCC data condition 1
        data_cond1 <- filt.data.list()[[isolate({input$sel_cond1})]]
        # CCC data condition 2
        data_cond2 <- filt.data.list()[[isolate({input$sel_cond2})]]
        
        data_cond3 <- NULL
        lab_c3 <- NULL
        if(!(input$sel_cond3 %in% c("none", isolate({input$sel_cond1}), isolate({input$sel_cond2})))){
          # CCC data condition 3
          data_cond3 <- filt.data.list()[[input$sel_cond3]]
          lab_c3 <- db.names[[input$sel_cond3]]
        } 
        
          
        
        
        # Get table with only unique couplets
        rv$uni_couplets_tab <- getDistinctCouplets(data_cond1 = data_cond1,
                                                data_cond2 = data_cond2,
                                                data_cond3 = data_cond3, 
                                                lab_c1 = db.names[[isolate({input$sel_cond1})]],
                                                lab_c2 = db.names[[isolate({input$sel_cond2})]],
                                                lab_c3 = lab_c3)
        
        output$uni_couplets_table <- DT::renderDT({
          rv$uni_couplets_tab
        }, filter = list(position = 'top', clear = FALSE),
        options = list(scrollX= TRUE, scrollCollapse = TRUE, processing = FALSE),
        escape = FALSE)
        
        # Using a datatable proxy to manipulate the object
        proxy <- DT::dataTableProxy("uni_couplets_table")
        
        # Clear rows button
        observeEvent(input$clear_rows, {
          proxy %>% selectRows(NULL)
        })
        
        # Download table
        output$download_geneTab <- downloadHandler(
          filename = function() {
            paste0("MC_", db.names[[input$sel_cond1]], "VS", db.names[[input$sel_cond2]], "unique_couplets_table.csv")
          },
          content = function(file) {
            write.csv(rv$uni_couplets_tab, file, quote = TRUE, row.names = FALSE)
          }
        )
        
        ### Gene-verse Info barplot with # unique/shared
        
        # output$info_barplot <- renderPlot({
        #   
        # })
        
        
          
      })
      
      ####---- Gene-verse based Dotplot
      
      observeEvent(input$uni_couplets_table_rows_selected, {
        if(length(input$uni_couplets_table_rows_selected) > 0){
          data.dotplot <- reactive({
            rv$uni_couplets_tab[input$uni_couplets_table_rows_selected,]
          })
          cluster.list.dot <- reactive({getClusterA_Names(data.dotplot())})


          # generate UI dotplot
          output$dotplot.ui <- renderUI({
            sidebarLayout(
              sidebarPanel(width = 3,
                           checkboxGroupInput(session$ns("cluster_selected_dotplot"),
                                              label = "Sender clusters:",
                                              choices = cluster.list.dot(),
                                              selected = names(cluster.list.dot()),
                                              inline = FALSE)


              ),
              mainPanel(width = 9,
                        uiOutput(session$ns("unique.dotplot.ui"))

              )
            )

          })

          # React to checkbox
          data.dotplot.filt <- reactive({
            req(data.dotplot())
            data.dotplot() %>%
              filter(clustA %in% input$cluster_selected_dotplot)
          })
          # get dotplot
          unique_dotplot <- reactive({
            req(data.dotplot.filt())
            getUniqueDotplot(data.dotplot.filt(), 
                             clust.order = unique(data.dotplot.filt()$clustA))
          })

          
          # get height size for dotplot
          n_rows_dot <- reactive({
            req(data.dotplot.filt())
            clust_p <- unite(data.dotplot.filt(), col = "clust_p", clustA:clustB)
            n_rows_dot <- length(unique(clust_p$clust_p))
            n_rows_dot
          })





          # generate UI plot
          output$unique.dotplot.ui <- renderUI({
            plotOutput(session$ns("unique.dotplot"),
                       height = max(500, 30*n_rows_dot())) %>% withSpinner()
          })
          # generate plot
          output$unique.dotplot <- renderPlot({
            unique_dotplot()
          })



          # generate download button handler
          output$download_dotplot_tiff <- downloadHandler(
            filename = function() {
              paste0("MC_", db.names[[input$sel_cond1]], "VS", db.names[[input$sel_cond2]], "unique_couplets_dotplot.tiff")
            },
            content = function(file) {
              tiff(file, height = max(500, 30*n_rows_dot()))
              plot(unique_dotplot())
              dev.off()
            }
          )
          # Download dotplot (pdf)
          output$download_dotplot_pdf <- downloadHandler(
            filename = function() {
              paste0("MC_", db.names[[input$sel_cond1]], "VS", db.names[[input$sel_cond2]], "unique_couplets_dotplot.pdf")
            },
            content = function(file) {

              ggsave(filename = file,
                     plot = unique_dotplot(),
                     device = "pdf", width = 12, height = 20, units = "cm", scale = 2)
            }
          )
          
          ####--------------------- Generate Pie Chart
          output$piechart.ui <- renderUI({
            plotOutput(session$ns("unique.piechart"))
          })
          
          output$unique.piechart <- renderPlot({
            req(data.dotplot.filt())
            getPieChart(data.dotplot.filt())
          })
          
          
          # generate download button handler
          output$download_pie_tiff <- downloadHandler(
            filename = function() {
              paste0("MC_", db.names[[input$sel_cond1]], "VS", db.names[[input$sel_cond2]], "unique_couplets_piechart.tiff")
            },
            content = function(file) {
              tiff(file, height = max(500, 30*n_rows_dot()))
              plot(getPieChart(data.dotplot.filt()))
              dev.off()
            }
          )
          # Download dotplot (pdf)
          output$download_pie_pdf <- downloadHandler(
            filename = function() {
              paste0("MC_", db.names[[input$sel_cond1]], "VS", db.names[[input$sel_cond2]], "unique_couplets_piechart.pdf")
            },
            content = function(file) {
              
              ggsave(filename = file,
                     plot = getPieChart(data.dotplot.filt()),
                     device = "pdf", width = 12, height = 20, units = "cm", scale = 2)
            }
          )

        } # end if
        })
      
      
      #### Function-verse based
      observeEvent(input$go, {
        # CCC data condition 1
        data_cond1 <- filt.data.list()[[isolate({input$sel_cond1})]]
        # CCC data condition 2
        data_cond2 <- filt.data.list()[[isolate({input$sel_cond2})]]
        
        
        
        # Annotated int-pairs BY functional terms matrices
        # Cond 1 
        annot_cond1 <- func.annot.mat.list()[[input$sel_cond1]]
        # Cond 2 
        annot_cond2 <- func.annot.mat.list()[[input$sel_cond2]]
        
        # Ranked functional terms
        # Cond1
        ranked_terms_cond1 <- ranked.terms.list()[[input$sel_cond1]]
        # Cond2
        ranked_terms_cond2 <- ranked.terms.list()[[input$sel_cond2]]
        
        lab_c1 <- db.names[[isolate({input$sel_cond1})]]
        lab_c2 <- db.names[[isolate({input$sel_cond2})]]
        
        data_cond3 <- NULL
        lab_c3 <- NULL
        annot_cond3 <- NULL
        ranked_terms_cond3 <- NULL
        if(!(input$sel_cond3 %in% c("none", isolate({input$sel_cond1}), isolate({input$sel_cond2})))){
          # Data condition 3
          data_cond3 <- filt.data.list()[[input$sel_cond3]]
          lab_c3 <- db.names[[input$sel_cond3]]
          annot_cond3 <- func.annot.mat.list()[[input$sel_cond3]]
          ranked_terms_cond3 <- ranked.terms.list()[[input$sel_cond3]]
        } 
        
        # Get table with significant functions annotated to unique int-pairs of each condition
        out <- tryCatch({
          signFunc_table_unique <- getSignif_table(data_cond1 = data_cond1,
                                                   data_cond2 = data_cond2,
                                                   data_cond3 = data_cond3,
                                                   lab_c1 = lab_c1,
                                                   lab_c2 = lab_c2,
                                                   lab_c3 = lab_c3,
                                                   annot_cond1,
                                                   annot_cond2,
                                                   annot_cond3)
        },
        error = function(cond){
          message("error in signFunc_table_unique")
        },
        warning = function(cond){
          message("error in signFunc_table_unique")
        })
        
        
        
        
        
        
        output$download_sub_annot <- downloadHandler(
          filename = function() {
            "sub_annot.rds"
          },
          content = function(file) {
            
            saveRDS(sub_annot, file = file)
          }
        )
        
        output$download_unique_intpa <- downloadHandler(
          filename = function() {
            "uni_intp.rds"
          },
          content = function(file) {
            
            saveRDS(unique_intpairs, file = file)
          }
        )
        
        output$funcTab1 <- DT::renderDT({
          signFunc_table_unique
        }, filter = list(position = 'top', clear = FALSE),
        options = list(scrollX= TRUE, scrollCollapse = TRUE, processing = FALSE),
        escape = FALSE)
        
        
        #signFun <- pvalue_df[pvalue_df$p_value <= input_maxPval,]

        # output$funcTab2 <- DT::renderDT({
        #  signFunc_table_unique
        # }, filter = list(position = 'top', clear = FALSE),
        # options = list(scrollX= TRUE, scrollCollapse = TRUE, processing = FALSE),
        # escape = FALSE)
        
        # output$funcTab3 <- DT::renderDT({
        #   
        # }, filter = list(position = 'top', clear = FALSE),
        # options = list(scrollX= TRUE, scrollCollapse = TRUE, processing = FALSE),
        # escape = FALSE)
      
        
        
        
        
        
        })
      



      
      
      
      
      
   
    
    
    
    
    
    
    
    
    
  })
}
    

