#' gene_verse UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList uiOutput downloadButton
#' @importFrom shinydashboard infoBoxOutput
#' @importFrom DT DTOutput
#' @importFrom shinycssloaders withSpinner
mod_gene_verse_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(width = 3,
             h2("Gene-verse")
             )
    ),
    fluidRow(
      infoBoxOutput(ns("prot_total_info")),
      infoBoxOutput(ns("ligand_info")),
      infoBoxOutput(ns("receptor_info"))
      
    ),
    fluidRow(
      box(title = "Filters",
          width = 12,
          status = "success",
          solidHeader = TRUE,
          collapsible = TRUE,
          column(width = 5,
                 uiOutput(ns("ann_strategy_checkbox_ui"))
          )
      )
    ),
    fluidRow(
      tabBox(
        id = ns('gene_verse_tabbox'),
        width = 12,
        height = "auto",
        tabPanel(h4("Table"),
                 column(2,
                        downloadButton(ns("download_geneTab"), "Download Table"),
                        ),
                 column(2,
                        actionButton(ns("clear_rows"), "Clear Rows")
                        ),
                 br(),
                 br(),
                 DT::DTOutput(ns("gene_table")) %>% withSpinner()
        ),
        tabPanel(h4("Dot Plot"),
                 uiOutput(ns("dotplot.text.ui")),
                 uiOutput(ns("dotplot.ui")),
                 
        ),
        tabPanel(h4("All vs all"),
                 h4("Here you can compare dot plots generated for multiple conditions."),
                 p("Only the int-pairs/cluster-pairs that are unique to a 
                    certain condition are shown!"),
                 p("For each condition, please download the table output of 
                   the dot plot (in the previous tab) and upload it below."),
                 p("Ideally, similar sets of int-pairs and cluster-pairs should 
                   be considered. At least two conditions are required."),
                 sidebarLayout(
                   sidebarPanel(width = 3,
                                textInput(ns("cond1_lab"), "Condition #1 label"),
                                textInput(ns("cond2_lab"), "Condition #2 label"),
                                textInput(ns("cond3_lab"), "Condition #3 label"),
                                hr(),
                                fileInput(ns("csv_cond1"), 
                                          "Dotplot data, condition #1 (csv)", 
                                          multiple = FALSE, 
                                          accept = ".csv"),
                                fileInput(ns("csv_cond2"), 
                                          "Dotplot data, condition #2 (csv)", 
                                          multiple = FALSE, 
                                          accept = ".csv"),
                                fileInput(ns("csv_cond3"), 
                                          "Dotplot data, condition #3 (csv)", 
                                          multiple = FALSE, 
                                          accept = ".csv"),
                                actionButton(ns("plot_dotplot"), "Plot!"),
                                hr(),
                                downloadButton(ns("download_dotplot_pdf"), 
                                               "Download Dotplot (pdf)"),
                                downloadButton(ns("download_dotplot_tiff"), 
                                               "Download Dotplot (tiff)")
                                
                                
                   ),
                   mainPanel(width = 9,
                             plotOutput(ns("dotplot_unique"), height = 1000) 
                   )
                 )
     
        )
      )
    )
 
  )
}
    
#' gene_verse Server Functions
#'
#' @noRd 
#' @importFrom shiny downloadHandler
#' @importFrom shinydashboard renderInfoBox infoBox
#' @importFrom DT renderDT dataTableProxy selectRows
#' @importFrom utils write.csv
#' @importFrom colourpicker colourInput
#' @importFrom tidyr unite
#' @importFrom grDevices tiff pdf dev.off
#' @importFrom ggplot2 ggsave

mod_gene_verse_server <- function(id, filt.data){
  moduleServer( id, function(input, output, session){
    
    
    rv <- reactiveValues(gene.filt.data = NULL, gene.table = NULL)
  
    
    observeEvent(filt.data(), {
      # Initialize filters 
      if("annotation_strategy" %in% colnames(filt.data())){
        # List of sources from which the interactions are annotated 
        sources.list <- as.list(unique(unlist(strsplit(
          as.character(filt.data()$annotation_strategy), ","))))
        names(sources.list) <- unlist(sources.list)
        output$ann_strategy_checkbox_ui <- renderUI(
          checkboxGroupInput(session$ns("ann_strategy_checkbox"),
                             label = h4("Annotation Sources for Interaction Pairs"),
                             choices = sources.list,
                             selected = names(sources.list),
                             inline = TRUE)
        )
      } else {
        # No filtering options available 
        output$ann_strategy_checkbox_ui <- renderUI(
          h4(textOutput(session$ns("no_filters")))
        )
        output$no_filters <- renderText({
          "There are no filtering options on genes available for your dataset!"
        })
      }
    })
    
    # Update gene.filt.data when filtering + create gene table to display
    observeEvent(c(filt.data(), input$ann_strategy_checkbox), {
      req(filt.data())
      if("annotation_strategy" %in% colnames(filt.data())){
        # create gene table to display
        gene.tab <- getGeneTable(filt.data())
        rv$gene.table <- gene.tab[grep(paste(input$ann_strategy_checkbox, 
                                                 collapse = "|"), 
                                           gene.tab$annotation_strategy),]
        
        # Update filtered data matrix to return
        rv$gene.filt.data <- filt.data()[grep(
          paste(input$ann_strategy_checkbox, collapse = "|"), 
          filt.data()$annotation_strategy), ]
      }
      else {
        rv$gene.table <- getGeneTable(filt.data())
      }
      

    })
    
    
    
    # unique proteins (and complexes) that participate in an interaction
    prot.unique <- reactive({
      req(rv$gene.table)
      unique(unlist(strsplit(as.character(rv$gene.table$int_pair), " & ")))
      })
    output$prot_total_info <- renderInfoBox({
      req(rv$gene.table)
      infoBox(h4("Proteins & Complexes"), 
              value = length(prot.unique()), 
              icon = icon("dna"), 
              fill = FALSE, 
              color = "light-blue")
    })
    output$ligand_info <- renderInfoBox({
      req(rv$gene.table)
      infoBox(h4("Ligands"), 
              value = getNumLR(rv$gene.table, type = "L"), 
              icon = icon("shapes"), 
              fill = FALSE, 
              color = "orange")
    })
    output$receptor_info <- renderInfoBox({
      req(rv$gene.table)
      infoBox(h4("Receptors"), 
              value = getNumLR(rv$gene.table, type = "R"), 
              icon = icon("hands"), 
              fill = FALSE, 
              color = "purple")
    })
    
 
    
    ####--- Gene Table ---####
    
    # Plot table
    output$gene_table <- DT::renderDT({
      req(rv$gene.table)
      rv$gene.table
    }, filter = list(position = 'top', clear = FALSE),
    options = list(scrollX= TRUE, scrollCollapse = TRUE, processing = FALSE),
    escape = FALSE)
    
    # Using a datatable proxy to manipulate the object
    proxy <- DT::dataTableProxy("gene_table")
    
    # Clear rows button
    observeEvent(input$clear_rows, {
      proxy %>% selectRows(NULL)
    })
    
    # Download table
    output$download_geneTab <- downloadHandler(
      filename = function() {
        "Gene-verse_table.csv"
      },
      content = function(file) {
        write.csv(rv$gene.table, file, quote = TRUE, row.names = FALSE)
      }
    )
    
    
    
    
    ####--- Dotplot ---####
    output$no_genes_selected <- renderText({
      "Select the int-pairs from the Table to see them in a Dot Plot!"
      })

    output$dotplot.text.ui <- renderUI({
      if(length(input$gene_table_rows_selected) == 0){
        h3(textOutput(session$ns("no_genes_selected")))
      } else{
        NULL
      }
    })




    observeEvent(input$gene_table_rows_selected, {
      if(length(input$gene_table_rows_selected) > 0){
        intpair_selected <- reactive({
          as.character(rv$gene.table$int_pair[input$gene_table_rows_selected])
          })
        data.dotplot <- reactive({
          rv$gene.filt.data %>%
            filter(int_pair %in% intpair_selected())
        })
        cluster.list.dot <- reactive({getClusterNames(data.dotplot())})


        # generate UI
        output$dotplot.ui <- renderUI({
          sidebarLayout(
            sidebarPanel(width = 3,
                         checkboxGroupInput(session$ns("cluster_selected_dotplot"),
                                            label = "Clusters:",
                                            choices = cluster.list.dot(),
                                            selected = names(cluster.list.dot()),
                                            inline = FALSE),
                         colourInput(session$ns("col_high"), 
                                     label = "Color high score:", 
                                     value = "#131780"),
                         colourInput(session$ns("col_low"), 
                                     label = "Color low score:", 
                                     value = "aquamarine"),
                         hr(),
                         downloadButton(session$ns("download_dotplot_tiff"),
                                        "Download DotPlot (tiff)"),
                         downloadButton(session$ns("download_dotplot_data"),
                                        "Download data (csv)"),

            ),
            mainPanel(width = 9,
                      uiOutput(session$ns("gene.dotplot.ui"))

            )
          )

        })
        
        # React to checkbox
        data.dotplot.filt <- reactive({
          req(data.dotplot())
          data.dotplot() %>%
            filter(clustA %in% input$cluster_selected_dotplot & 
                     clustB %in% input$cluster_selected_dotplot)
        })
        # get dotplot
        dot_list <- reactive({
          req(data.dotplot.filt())
          getDotPlot_selInt(data.dotplot.filt(), 
                            clust.order = unique(data.dotplot.filt()$clustA),
                            low_color = input$col_low, 
                            high_color = input$col_high)
        })
        # get height size for dotplot
        n_rows_dot <- reactive({
          req(data.dotplot.filt())
          clust_p <- unite(data.dotplot.filt(), col = "clust_p", clustA:clustB)
          n_rows_dot <- length(unique(clust_p$clust_p))
          n_rows_dot
        })





        # generate UI plot
        output$gene.dotplot.ui <- renderUI({
          plotOutput(session$ns("gene.dotplot"), 
                     height = max(500, 30*n_rows_dot())) %>% withSpinner()
        })
        # generate plot
        output$gene.dotplot <- renderPlot({
          req(dot_list())
          dot_list()$p
        })
        # generate download button handler
        output$download_dotplot_tiff <- downloadHandler(
          filename = function() {
            "Gene-verse_dotplot.tiff"
          },
          content = function(file) {
            tiff(file, height = max(500, 30*n_rows_dot()))
            plot(dot_list()$p)
            dev.off()
          }
        )
        
        # generate download button handler
        output$download_dotplot_data <- downloadHandler(
          filename = function() {
            "Gene-verse_dotplot_data.csv"
          },
          content = function(file) {
            write.csv(dot_list()$data_dot, file, quote = TRUE, row.names = FALSE)
          }
        )



      }
    })
    
    ########---------- all vs all ------#######
    
    ### Dotplot of unique int-pairs/cluster-pairs
    observeEvent(input$plot_dotplot, {
      if(is.null(input$csv_cond1)){
        shinyalert(text = "Please select a csv file for condition 1", 
                   type = "error",
                   showCancelButton = FALSE)
      } 
      if(is.null(input$csv_cond2)){
        shinyalert(text = "Please select a csv file for condition 2", 
                   type = "error",
                   showCancelButton = FALSE)
      } 
      if(input$cond1_lab == ""){
        shinyalert(text = "Please specify a label for condition 1", 
                   type = "error",
                   showCancelButton = FALSE)
      }
      if(input$cond2_lab == ""){
        shinyalert(text = "Please specify a label for condition 2", 
                   type = "error",
                   showCancelButton = FALSE)
      }
      if(input$cond3_lab == "" & !is.null(input$csv_cond3)){
        shinyalert(text = "Please specify a label for condition 3", 
                   type = "error",
                   showCancelButton = FALSE)
      }
      
      
      
      req(input$csv_cond1, input$csv_cond2,
          input$cond1_lab,input$cond2_lab)
      file_c1 <- input$csv_cond1
      tab_c1 <- read.csv(file_c1$datapath)
      
      file_c2 <- input$csv_cond2
      tab_c2 <- read.csv(file_c2$datapath)
      
      if(!is.null(input$csv_cond3)){
        file_c3 <- input$csv_cond3
        tab_c3 <- read.csv(file_c3$datapath)
      } else {
        tab_c3 <- NULL
      }
      
      if(!all(c("int_pair", "cluster_pair") %in% colnames(tab_c1))){
        shinyalert(text = "Looks like the csv file for condition 1 isn't the right one!", 
                   type = "error",
                   showCancelButton = FALSE)
        tab_c1 <- NULL
      }
      if(!all(c("int_pair", "cluster_pair") %in% colnames(tab_c2))){
        shinyalert(text = "Looks like the csv file for condition 2 isn't the right one!", 
                   type = "error",
                   showCancelButton = FALSE)
        tab_c2 <- NULL
      }
      if(!is.null(tab_c3) & !all(c("int_pair", "cluster_pair") %in% colnames(tab_c3))){
        shinyalert(text = "Looks like the csv file for condition 3 isn't the right one!", 
                   type = "error",
                   showCancelButton = FALSE)
        tab_c3 <- NULL
      }
      
      req(tab_c1, tab_c2)
      
      tab_c1$condition <- input$cond1_lab
      tab_c2$condition <- input$cond2_lab
      
      data_dotplot <- rbind(tab_c1, tab_c2)
      data_dotplot$condition <- factor(data_dotplot$condition, 
                                       levels = c(input$cond1_lab,input$cond2_lab))
      if(!is.null(tab_c3)){
        tab_c3$condition <- input$cond3_lab
        data_dotplot <- rbind(data_dotplot, tab_c3)
        data_dotplot$condition <- factor(data_dotplot$condition, 
                                         levels = c(input$cond1_lab,input$cond2_lab, input$cond3_lab))
      }
      

      
      unique_dotplot <- getUniqueDotplot(data_dotplot)
      
      
      output$dotplot_unique <- renderPlot({
        unique_dotplot
      })
      
      # Download dotplot (tiff)
      output$download_dotplot_tiff <- downloadHandler(
        filename = function() {
          paste0("Gene-verse_allvsall_unique_dotplot.tiff")
        },
        content = function(file) {
          
          tiff(file, width = 700, height = 1000)
          plot(unique_dotplot)
          dev.off()
        }
      )
      # Download dotplot (pdf)
      output$download_dotplot_pdf <- downloadHandler(
        filename = function() {
          paste0("Gene-verse_allvsall_unique_dotplot.pdf")
        },
        content = function(file) {
          
          ggsave(filename = file, 
                 plot = unique_dotplot,
                 device = "pdf", width = 12, height = 20, units = "cm", scale = 2)
        }
      )
      
    })

    
    return(rv)
 
  })
}
    
