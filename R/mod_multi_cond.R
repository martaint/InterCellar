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
                   actionButton(ns("download_backbar_pdf"), 
                                  "Barplot (pdf)", icon = icon("download")),
                   actionButton(ns("download_backbar_tiff"), 
                                  "Barplot (tiff)", icon = icon("download")),
                   hr(),
                   h4("Radar plot"),
                   
                   uiOutput(ns("radar_vp_ui")),
                   
                   actionButton(ns("download_radar_pdf"), 
                                  "Radar plot (pdf)", icon = icon("download")),
                   actionButton(ns("download_radar_tiff"), 
                                  "Radar plot (tiff)", icon = icon("download"))
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
                   
                   
                   actionButton(ns("download_dotplot_pdf"), 
                                  "Dotplot (pdf)", icon = icon("download")),
                   actionButton(ns("download_dotplot_tiff"), 
                                  "Dotplot (tiff)", icon = icon("download")),
                   br(),
                   actionButton(ns("download_pie_pdf"), 
                                  "Piechart (pdf)", icon = icon("download")),
                   actionButton(ns("download_pie_tiff"), 
                                  "Piechart (tiff)", icon = icon("download"))
                   
                   
               ),
               
               tabBox(
                 id = 'gene-verse_tabbox',
                 width = 9,
                 height = "auto",
                 
                 tabPanel(h4("Table"),
                          h4("Select int-pairs/cluster-pairs couplets from the Table to generate a DotPlot!"),
                          column(2,
                                 actionButton(ns("download_geneTab"), "Table (csv)", icon = icon("download")),
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
                   
                   uiOutput(ns("chooseCond_signF_ui")),
                   numericInput(ns("maxPval"),
                                label = "Maximum significant p value",
                                value = 0.05,
                                min = 0, max = 1, step = 0.01)
                   
                   
                   
               ),
               
               tabBox(
                 id = 'function-verse_tabbox',
                 width = 9,
                 height = "auto",
                 tabPanel(h4("Table"),
                          h4("Select a significant Functional Term from the Table to generate a Sunburst Plot!"),
                          column(2,
                                 actionButton(ns("download_funcTab"), "Table (csv)", icon = icon("download")),
                          ),
                          br(),
                          br(),
                          
                          DT::DTOutput(ns("funcTab")) %>% withSpinner()
                 ),
                 tabPanel(h4("Sunburst Plot"),
                          uiOutput(ns("sunburst.ui"))
                          
                 )
               )
        )
      )# fluidrow
      
    ) # fluidRow
 
  ) # tagList
}
    
#' multi_cond Server Functions
#'
#' @importFrom fmsb radarchart
#' @noRd 
mod_multi_cond_server <- function(id,
                                  input_sidebarmenu,
                                  db.list,
                                  cluster.colors,
                                  filt.data.list,
                                  func.annot.mat.list,
                                  ranked.terms.list,
                                  out_folder,
                                  output_tags){
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
    
    observeEvent({req(db.list())
      db.list()}, {
        # reverse list elements <-> names
        rv$db.names <- as.list(names(db.list()))
        names(rv$db.names) <- unlist(as.character(db.list()))
        
        
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
      
    
      
      rv$b2b_barplot <- getBack2BackBarplot(tab_c1 = barplotDF1(),
                                         tab_c2 = barplotDF2(),
                                         lab_c1 = rv$db.names[[isolate({input$sel_cond1})]],
                                         lab_c2 = rv$db.names[[isolate({input$sel_cond2})]])


      output$backbar <- renderPlot({
        rv$b2b_barplot
      })
      
      
      
      
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
      
      
    })
    
    # Download BackBar (tiff)
    observeEvent(input$download_backbar_tiff, {
      dir.create(file.path(out_folder(), 
                           paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]])), 
                 showWarnings = FALSE)
      file <- file.path(out_folder(), paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]]),
                        "Back2back_Barplot.tiff")
      tiff(file, width = 700)
      plot(rv$b2b_barplot)
      dev.off()
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
    })
    
    # Download BackBar (pdf)
    observeEvent(input$download_backbar_pdf, {
      dir.create(file.path(out_folder(), 
                           paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]])), 
                 showWarnings = FALSE)
      file <- file.path(out_folder(), paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]]),
                        "Back2back_Barplot.pdf")
      pdf(file)
      plot(rv$b2b_barplot)
      dev.off()
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
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
          lab_c3 <- rv$db.names[[input$sel_cond3]]
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
        
        rv$radar_df <- getRadar_df(tab_c1 = barplotDF_VP_1(),
                                   tab_c2 = barplotDF_VP_2(),
                                   tab_c3 = barplotDF_VP_3(),
                                   lab_c1 = rv$db.names[[isolate({input$sel_cond1})]],
                                   lab_c2 = rv$db.names[[isolate({input$sel_cond2})]],
                                   lab_c3 = lab_c3)
        
        output$radar <- renderPlot({
          fmsb::radarchart(
            rv$radar_df, axistype = 1,
            # Customize the polygon
            pcol = c("#438ECC", "#E97778", "#00BA38"), 
            pfcol = scales::alpha(c("#438ECC", "#E97778", "#00BA38"), 0.5), plwd = 2, plty = 1,
            # Customize the grid
            cglcol = "grey", cglty = 1, cglwd = 0.8,
            # Customize the axis
            axislabcol = "grey30", 
            # Variable labels
            vlcex = 1.2, vlabels = colnames(rv$radar_df),
            caxislabels = round(seq(from = 0, to = rv$radar_df["max",1], length.out = 5)), 
            title = input$sel_vp_radar
          )
          legend(
            x = "bottomleft", legend = rownames(rv$radar_df[-c(1,2),]), horiz = FALSE,
            bty = "n", pch = 20 , col = c("#438ECC", "#E97778", "#00BA38"),
            text.col = "black", cex = 1, pt.cex = 1.5
          )
          
        })
        
        
        
        
      })
      
      # Download Radar (tiff)
      observeEvent(input$download_radar_tiff, {
        if(!(input$sel_cond3 %in% c("none", input$sel_cond1, input$sel_cond2))){ # 3 conditions
          dir.create(file.path(out_folder(), 
                               paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]], "VS", output_tags()[[input$sel_cond3]])), 
                     showWarnings = FALSE)
          file <- file.path(out_folder(), paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]], "VS", output_tags()[[input$sel_cond3]]),
                            paste(input$sel_vp_radar, "radarplot.tiff", sep = "_"))
        } else {
          dir.create(file.path(out_folder(), 
                               paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]])), 
                     showWarnings = FALSE)
          file <- file.path(out_folder(), paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]]),
                            paste(input$sel_vp_radar, "radarplot.tiff", sep = "_"))
        }
        
        
        tiff(file, width = 600)
        fmsb::radarchart(
          rv$radar_df, axistype = 1,
          # Customize the polygon
          pcol = c("#438ECC", "#E97778", "#00BA38"), 
          pfcol = scales::alpha(c("#438ECC", "#E97778", "#00BA38"), 0.5), plwd = 2, plty = 1,
          # Customize the grid
          cglcol = "grey", cglty = 1, cglwd = 0.8,
          # Customize the axis
          axislabcol = "grey30", 
          # Variable labels
          vlcex = 1.2, vlabels = colnames(rv$radar_df),
          caxislabels = round(seq(from = 0, to = rv$radar_df["max",1], length.out = 5)), 
          title = input$sel_vp_radar
        )
        legend(
          x = "bottomleft", legend = rownames(rv$radar_df[-c(1,2),]), horiz = FALSE,
          bty = "n", pch = 20 , col = c("#438ECC", "#E97778", "#00BA38"),
          text.col = "black", cex = 1, pt.cex = 1.5
        )
        
        dev.off()
        
        shinyalert(text = paste("Saved!", file, sep = "\n"), 
                   type = "success",
                   showCancelButton = FALSE,
                   size = "m")
      })
      
      
      # Download Radar (pdf)
      observeEvent(input$download_radar_pdf, {
        if(!(input$sel_cond3 %in% c("none", input$sel_cond1, input$sel_cond2))){ # 3 conditions
          dir.create(file.path(out_folder(), 
                               paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]], "VS", output_tags()[[input$sel_cond3]])), 
                     showWarnings = FALSE)
          file <- file.path(out_folder(), paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]], "VS", output_tags()[[input$sel_cond3]]),
                            paste(input$sel_vp_radar, "radarplot.pdf", sep = "_"))
        } else {
          dir.create(file.path(out_folder(), 
                               paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]])), 
                     showWarnings = FALSE)
          file <- file.path(out_folder(), paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]]),
                            paste(input$sel_vp_radar, "radarplot.pdf", sep = "_"))
        }
        pdf(file)
        fmsb::radarchart(
          rv$radar_df, axistype = 1,
          # Customize the polygon
          pcol = c("#438ECC", "#E97778", "#00BA38"), 
          pfcol = scales::alpha(c("#438ECC", "#E97778", "#00BA38"), 0.5), plwd = 2, plty = 1,
          # Customize the grid
          cglcol = "grey", cglty = 1, cglwd = 0.8,
          # Customize the axis
          axislabcol = "grey30", 
          # Variable labels
          vlcex = 1.2, vlabels = colnames(rv$radar_df),
          caxislabels = round(seq(from = 0, to = rv$radar_df["max",1], length.out = 5)), 
          title = input$sel_vp_radar
        )
        legend(
          x = "bottomleft", legend = rownames(rv$radar_df[-c(1,2),]), horiz = FALSE,
          bty = "n", pch = 20 , col = c("#438ECC", "#E97778", "#00BA38"),
          text.col = "black", cex = 1, pt.cex = 1.5
        )
        dev.off()
        
        shinyalert(text = paste("Saved!", file, sep = "\n"), 
                   type = "success",
                   showCancelButton = FALSE,
                   size = "m")
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
          lab_c3 <- rv$db.names[[input$sel_cond3]]
        } 
        
          
        
        
        # Get table with only unique couplets
        rv$uni_couplets_tab <- getDistinctCouplets(data_cond1 = data_cond1,
                                                data_cond2 = data_cond2,
                                                data_cond3 = data_cond3, 
                                                lab_c1 = rv$db.names[[isolate({input$sel_cond1})]],
                                                lab_c2 = rv$db.names[[isolate({input$sel_cond2})]],
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
        

        
        
        
        
          
      })
      
      # Download table
      observeEvent(input$download_geneTab, {
        if(!(input$sel_cond3 %in% c("none", input$sel_cond1, input$sel_cond2))){ # 3 conditions
          dir.create(file.path(out_folder(), 
                               paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]], "VS", output_tags()[[input$sel_cond3]])), 
                     showWarnings = FALSE)
          file <- file.path(out_folder(), paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]], "VS", output_tags()[[input$sel_cond3]]),
                            "Unique_IntPair-ClustPair_couplets_table.csv")
        } else {
          dir.create(file.path(out_folder(), 
                               paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]])), 
                     showWarnings = FALSE)
          file <- file.path(out_folder(), paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]]),
                            "Unique_IntPair-ClustPair_couplets_table.csv")
        }
        
        write.csv(rv$uni_couplets_tab, file, quote = TRUE, row.names = FALSE)
        
        shinyalert(text = paste("Saved!", file, sep = "\n"), 
                   type = "success",
                   showCancelButton = FALSE,
                   size = "m")
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


          # Download Dotplot (tiff)
          observeEvent(input$download_dotplot_tiff, {
            if(!(input$sel_cond3 %in% c("none", input$sel_cond1, input$sel_cond2))){ # 3 conditions
              dir.create(file.path(out_folder(), 
                                   paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]], "VS", output_tags()[[input$sel_cond3]])), 
                         showWarnings = FALSE)
              file <- file.path(out_folder(), paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]], "VS", output_tags()[[input$sel_cond3]]),
                                "Unique_IntPair-ClustPair_coupl_dotplot.tiff")
            } else {
              dir.create(file.path(out_folder(), 
                                   paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]])), 
                         showWarnings = FALSE)
              file <- file.path(out_folder(), paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]]),
                                "Unique_IntPair-ClustPair_coupl_dotplot.tiff")
            }
            
            
            tiff(file, height = max(500, 30*n_rows_dot()))
            plot(unique_dotplot())
            dev.off()
            
            shinyalert(text = paste("Saved!", file, sep = "\n"), 
                       type = "success",
                       showCancelButton = FALSE,
                       size = "m")
          })
          
          
          # Download Dotplot (pdf)
          observeEvent(input$download_dotplot_pdf, {
            if(!(input$sel_cond3 %in% c("none", input$sel_cond1, input$sel_cond2))){ # 3 conditions
              dir.create(file.path(out_folder(), 
                                   paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]], "VS", output_tags()[[input$sel_cond3]])), 
                         showWarnings = FALSE)
              file <- file.path(out_folder(), paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]], "VS", output_tags()[[input$sel_cond3]]),
                                "Unique_IntPair-ClustPair_coupl_dotplot.pdf")
            } else {
              dir.create(file.path(out_folder(), 
                                   paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]])), 
                         showWarnings = FALSE)
              file <- file.path(out_folder(), paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]]),
                                "Unique_IntPair-ClustPair_coupl_dotplot.pdf")
            }
            
            ggsave(filename = file,
                   plot = unique_dotplot(),
                   device = "pdf", width = 12, height = 20, units = "cm", scale = 2)
            
            shinyalert(text = paste("Saved!", file, sep = "\n"), 
                       type = "success",
                       showCancelButton = FALSE,
                       size = "m")
          })
          
        

          
          ####--------------------- Generate Pie Chart
          output$piechart.ui <- renderUI({
            plotOutput(session$ns("unique.piechart"))
          })
          
          output$unique.piechart <- renderPlot({
            req(data.dotplot.filt())
            getPieChart(data.dotplot.filt())
          })
          
          # Download Pie (tiff)
          observeEvent(input$download_pie_tiff, {
            if(!(input$sel_cond3 %in% c("none", input$sel_cond1, input$sel_cond2))){ # 3 conditions
              dir.create(file.path(out_folder(), 
                                   paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]], "VS", output_tags()[[input$sel_cond3]])), 
                         showWarnings = FALSE)
              file <- file.path(out_folder(), paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]], "VS", output_tags()[[input$sel_cond3]]),
                                "Unique_IntPair-ClustPair_coupl_piechart.tiff")
            } else {
              dir.create(file.path(out_folder(), 
                                   paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]])), 
                         showWarnings = FALSE)
              file <- file.path(out_folder(), paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]]),
                                "Unique_IntPair-ClustPair_coupl_piechart.tiff", sep = "_")
            }
            
            
            tiff(file, height = max(500, 30*n_rows_dot()))
            plot(getPieChart(data.dotplot.filt()))
            dev.off()
            
            shinyalert(text = paste("Saved!", file, sep = "\n"), 
                       type = "success",
                       showCancelButton = FALSE,
                       size = "m")
          })
          
          # Download Pie (pdf)
          observeEvent(input$download_pie_pdf, {
            if(!(input$sel_cond3 %in% c("none", input$sel_cond1, input$sel_cond2))){ # 3 conditions
              dir.create(file.path(out_folder(), 
                                   paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]], "VS", output_tags()[[input$sel_cond3]])), 
                         showWarnings = FALSE)
              file <- file.path(out_folder(), paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]], "VS", output_tags()[[input$sel_cond3]]),
                                "Unique_IntPair-ClustPair_coupl_piechart.pdf")
            } else {
              dir.create(file.path(out_folder(), 
                                   paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]])), 
                         showWarnings = FALSE)
              file <- file.path(out_folder(), paste0("InterCellar_results_", output_tags()[[input$sel_cond1]], "VS", output_tags()[[input$sel_cond2]]),
                                "Unique_IntPair-ClustPair_coupl_piechart.pdf", sep = "_")
            }
            
            
            ggsave(filename = file,
                   plot = getPieChart(data.dotplot.filt()),
                   device = "pdf", width = 12, height = 20, units = "cm", scale = 2)
            
            shinyalert(text = paste("Saved!", file, sep = "\n"), 
                       type = "success",
                       showCancelButton = FALSE,
                       size = "m")
          })
          

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
        
        lab_c1 <- rv$db.names[[isolate({input$sel_cond1})]]
        lab_c2 <- rv$db.names[[isolate({input$sel_cond2})]]
        
        data_cond3 <- NULL
        lab_c3 <- NULL
        annot_cond3 <- NULL
        ranked_terms_cond3 <- NULL
        if(!(input$sel_cond3 %in% c("none", isolate({input$sel_cond1}), isolate({input$sel_cond2})))){
          # Data condition 3
          data_cond3 <- filt.data.list()[[input$sel_cond3]]
          lab_c3 <- rv$db.names[[input$sel_cond3]]
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
          message("err")
        },
        warning = function(cond){
          message("war")
        })
        
        
        output$chooseCond_signF_ui <- renderUI({
          req(signFunc_table_unique)
          cond_list <- as.list(unique(signFunc_table_unique$condition))
            selectInput(ns("chooseCond_signF"),
                        label = "Choose Condition:",
                        choices = cond_list,
                        multiple = FALSE
            )
        })
        
        filt_signFunc_tab <- reactive({
          req(signFunc_table_unique)
          signFunc_table_unique %>%
            filter(p_value <= input$maxPval) %>%
            filter(condition %in% input$chooseCond_signF) %>%
            arrange(`p_value`)
        })
        
        
        
        
        
        
        output$funcTab <- DT::renderDT({
          filt_signFunc_tab()
        }, filter = list(position = 'top', clear = FALSE),
        options = list(scrollX= TRUE, scrollCollapse = TRUE, processing = FALSE),
        escape = FALSE, selection = 'single')
        
        output$download_funcTab <- downloadHandler(
          filename = function() {
            paste0("MC_", input$chooseCond_signF, "signFuncTerms_table.csv")
          },
          content = function(file) {
            write.csv(filt_signFunc_tab(), file, quote = TRUE, row.names = FALSE)
          }
        )
        
        ### Plot Sunburst
        func_selected <- reactive({
          as.character(filt_signFunc_tab()$functionalTerm[
            input$funcTab_rows_selected])})
        
        int_p_fun <- reactive({
          int_list <- as.character(filt_signFunc_tab()$int_pair_list[
            input$funcTab_rows_selected])
          int_list <- unlist(strsplit(int_list, split=","))
          return(int_list)
        })
        
        sel.data.sunburst <- reactive({
          db <- filt.data.list()[[db.list()[[which(names(db.list()) == input$chooseCond_signF)]]]]
          db <- db %>%
            filter(int_pair %in% int_p_fun())
          return(db)
        })
          
         
          
        output$sunburst.ui <- renderUI({
          if(length(input$funcTab_rows_selected) == 0){
            NULL
          } else{
            sidebarLayout(
              sidebarPanel(width = 4,
                           h4("Selected Functional Term:"),
                           textOutput(ns("sel_fun_text")),
                           br(),
                           radioButtons(session$ns("num_or_weight_radio"),
                                        label = "Show",
                                        choices = list("Number of interactions" = "n_int",
                                                       "Weighted number of interactions (by score)" = "weighted"),
                           ),
                           br(),
                           h4("Annotated unique IntPairs:"),
                           br(),
                           DT::DTOutput(ns("annot_intp_table")),
              ),
              mainPanel(width = 8,
                        downloadButton(ns("download_sunburst"),
                                       "Download Plot"),
                        plotlyOutput(ns("sunburst.plot")) %>% withSpinner()
                        
              )
            )
          }
        })
        
        output$sel_fun_text <- renderText({
          func_selected()
        })
        output$annot_intp_table <- DT::renderDT({
          data.frame(int_pair = int_p_fun())
        }, options = list(scrollX= TRUE,
                          scrollCollapse = TRUE,
                          processing = FALSE), escape = FALSE, selection = 'none')
        

        output$sunburst.plot <- renderPlotly({
          getSunburst(sel.data.sunburst(), 
                      func_selected(), 
                      int_p_fun(), 
                      cluster.colors()[[db.list()[[which(names(db.list()) == input$chooseCond_signF)]]]],
                      input$num_or_weight_radio)
        })
        
        # Download sunburst
        output$download_sunburst <- downloadHandler(
          filename = function() {
            paste0("MC_", func_selected(), "_", input$chooseCond_signF, "_sunburst.html")},
          content = function(file) {
            fig <- getSunburst(sel.data.sunburst(), func_selected(), int_p_fun(), 
                               cluster.colors()[[db.list()[[which(names(db.list()) == input$chooseCond_signF)]]]],
                               input$num_or_weight_radio)
            htmlwidgets::saveWidget(fig, file = file, selfcontained = TRUE)
          }
        )
        
        
        
        
        
        
        })
      



      
      
      
      
      
   
    
    
    
    
    
    
    
    
    
  })
}
    

