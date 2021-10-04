#' gene_verse UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList uiOutput actionButton
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
          fluidRow(
            column(width = 12,
                   uiOutput(ns("geneverse_filters_ui"))
            )
          )
          
          
      )
    ),
    fluidRow(
      tabBox(
        id = ns('gene_verse_tabbox'),
        width = 12,
        height = "auto",
        tabPanel(h4("Table"),
                 h4("Select int-pairs from the Table to generate Dot Plot and Network!"),
                 column(2,
                        actionButton(ns("download_geneTab"),
                                     "Table",
                                     icon = icon("download"))
                        ),
                 column(2,
                        actionButton(ns("clear_rows"), "Clear Rows")
                        ),
                 br(),
                 br(),
                 DT::DTOutput(ns("gene_table")) %>% withSpinner()
        ),
        tabPanel(h4("Dot Plot"),
                 uiOutput(ns("dotplot.ui")),
                 
        ),
        
        tabPanel(h4("Network"),
                 sidebarLayout(
                   sidebarPanel(width = 3,
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
                                textInput(ns("net_tag"), label = "File tag for saving:",
                                          placeholder = "CXCL_family"),
                                actionButton(ns("download_network"),
                                             "Network (html)",
                                             icon = icon("download")),
                                hr(),
                                h4("Selected Int-Pair(s):"),
                                htmlOutput(ns("sel_intpair_text"))
                                
                   ),
                   mainPanel(width = 9,
                             visNetworkOutput(ns("gene.net"), 
                                              height = "550px") 
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
#' @importFrom shinyalert shinyalert

mod_gene_verse_server <- function(id, input_sidebarmenu, input.data, gene.table, out_folder){
  moduleServer( id, function(input, output, session){
    
    
    rv <- reactiveValues(okay_flag = FALSE,
                         gene.filt.data = NULL, 
                         gene.table_out = NULL, 
                         input.tool = NULL)
  
    observeEvent(input_sidebarmenu(), {
      if(input_sidebarmenu() == "gene-verse"){
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

    
    observeEvent({
      req(rv$okay_flag)
      input.data()}, {
        
        # Get input tool that was used
        if("annotation_strategy" %in% colnames(input.data())){
          rv$input.tool <- "cpdb"
        } else if("scSignalR_specific" %in% colnames(input.data())) {
          rv$input.tool <- "scsr"
        } else if("pathway_cellchat" %in% colnames(input.data())){
          rv$input.tool <- "cellchat"
        } else {
          rv$input.tool <- "custom"
        }
        
        
      })
    
    
    

    observeEvent(rv$input.tool, {
      req(input.data())
      # Update filters
      if(rv$input.tool == "cpdb"){
        # List of sources from which the interactions are annotated
        sources.list <- as.list(unique(unlist(strsplit(
          as.character(input.data()$annotation_strategy), ","))))
        names(sources.list) <- unlist(sources.list)
        output$geneverse_filters_ui <- renderUI(
          checkboxGroupInput(session$ns("ann_strategy_checkbox"),
                             label = h4("Annotation Sources for Interaction Pairs"),
                             choices = sources.list,
                             selected = names(sources.list),
                             inline = TRUE)
        )


      } else if(rv$input.tool == "scsr"){
        output$geneverse_filters_ui <- renderUI(
          radioButtons(session$ns("scsr_radio"),
                       label = h4("Select only int-pairs labelled as 'specific' by scSignalR:"),
                       choices = c("true", "false"),
                       selected = "false",
                       inline = TRUE)
        )
      } else if(rv$input.tool == "cellchat"){
        # List of pathways annotated by cellchat
        pathway.list <- as.list(unique(unlist(
          as.character(input.data()$pathway_cellchat))))
        names(pathway.list) <- unlist(pathway.list)
        # List of annotation for cellchat
        sources.list <- as.list(unique(unlist(
          as.character(input.data()$annotation_cellchat))))
        names(sources.list) <- unlist(sources.list)

        output$geneverse_filters_ui <- renderUI(
          tagList(
            column(width = 4,
                   selectInput(session$ns("cellchat_exclude_pathway"),
                               label = h4("Exclude selected Pathways"),
                               choices = c(list("None" = "none"), pathway.list),
                               selected = "none",
                               multiple = TRUE)
            ),
            column(width = 4,
                   checkboxGroupInput(session$ns("cellchat_ann_checkbox"),
                                      label = h4("Annotation"),
                                      choices = sources.list,
                                      selected = names(sources.list),
                                      inline = FALSE)
            ),
            column(width = 4,
                   br(),
                   actionButton(session$ns("apply_filt_cellchat"),
                                label = h4("Filter!"),
                                class = "btn-success")
            )

          )

        )
      } else if(rv$input.tool == "custom"){
        # No filtering options available
        output$geneverse_filters_ui <- renderUI(
          h4(textOutput(session$ns("no_filters")))
        )
        output$no_filters <- renderText({
          "There are no filtering options on genes available for your dataset!"
        })
        
        # Generate table for "no filters" cases
        rv$gene.table_out <- getGeneTable(input.data())
        rv$gene.filt.data <- input.data()
      

      }

    })

    # React to filters for CPDB
    observeEvent(input$ann_strategy_checkbox, {
      req(input.data())
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message= "Computing Gene Table", value = 0.5)

      # create gene table to display
      gene.tab <- getGeneTable(input.data())
      rv$gene.table_out <- gene.tab[grep(paste(input$ann_strategy_checkbox,
                                           collapse = "|"),
                                     gene.tab$annotation_strategy),]
      # Update filtered data matrix to return
      rv$gene.filt.data <- input.data()[grep(
        paste(input$ann_strategy_checkbox, collapse = "|"),
        input.data()$annotation_strategy), ]
    })

    # React to filters for SCSR
    observeEvent(input$scsr_radio, {
      req(input.data())

      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message= "Computing Gene Table", value = 0.5)

      # Update filtered data matrix to return
      if(input$scsr_radio == "true"){
        rv$gene.filt.data <- input.data() %>%
          filter(scSignalR_specific == "specific")
      } else {
        rv$gene.filt.data <- input.data()
      }

      rv$gene.table_out <- getGeneTable(rv$gene.filt.data)

    })

    # React to filters for cellchat
    observeEvent(input$apply_filt_cellchat, {
      req(input.data())
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message= "Computing Gene Table", value = 0.5)


      # Update filtered matrix
      if(length(input$cellchat_exclude_pathway) == 1 & input$cellchat_exclude_pathway == "none"){
        rv$gene.filt.data <- input.data() %>%
          filter(annotation_cellchat %in% input$cellchat_ann_checkbox)
      } else {
        rv$gene.filt.data <- input.data() %>%
          filter(!(pathway_cellchat %in% input$cellchat_exclude_pathway)) %>%
          filter(annotation_cellchat %in% input$cellchat_ann_checkbox)

      }
      rv$gene.table_out <- getGeneTable(rv$gene.filt.data)
    })
    
    
    



    # unique proteins (and complexes) that participate in an interaction
    prot.unique <- reactive({
      req(gene.table())
      unique(unlist(strsplit(as.character(gene.table()$int_pair), " & ")))
      })
    output$prot_total_info <- renderInfoBox({
      req(gene.table())
      infoBox(h4("Proteins & Complexes"),
              value = length(prot.unique()),
              icon = icon("dna"),
              fill = FALSE,
              color = "light-blue")
    })
    output$ligand_info <- renderInfoBox({
      req(gene.table())
      infoBox(h4("Ligands"),
              value = getNumLR(gene.table(), type = "L"),
              icon = icon("shapes"),
              fill = FALSE,
              color = "orange")
    })
    output$receptor_info <- renderInfoBox({
      req(gene.table())
      infoBox(h4("Receptors"),
              value = getNumLR(gene.table(), type = "R"),
              icon = icon("hands"),
              fill = FALSE,
              color = "purple")
    })



    ####--- Gene Table ---####

    # Plot table
    output$gene_table <- DT::renderDT({
      req(gene.table())
      gene.table()
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
    observeEvent(input$download_geneTab, {
      dir.create(file.path(out_folder(), "gene_verse"), showWarnings = FALSE)
      file <- file.path(out_folder(), "gene_verse", 
                        "IntPairs_table.csv")
      write.csv(gene.table(), file, quote = TRUE, row.names = FALSE)
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
    })
    
    


  
  
  ####--- Dotplot ---####
 
  observeEvent(input$gene_table_rows_selected, {
    if(length(input$gene_table_rows_selected) > 0){
      intpair_selected <- reactive({
        as.character(gene.table()$int_pair[input$gene_table_rows_selected])
        })
      data.dotplot <- reactive({
        rv$gene.filt.data %>%
          filter(int_pair %in% intpair_selected())
      })
      cluster.list.dot <- reactive({getClusterA_Names(data.dotplot())})
  
  
      # generate UI
      output$dotplot.ui <- renderUI({
        sidebarLayout(
          sidebarPanel(width = 3,
                       checkboxGroupInput(session$ns("cluster_selected_dotplot"),
                                          label = "Sender clusters:",
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
                       textInput(session$ns("dot_tag"), label = "File tag for saving:",
                                 placeholder = "CXCL_family"),
                       actionButton(session$ns("download_dotplot_tiff"),
                                      "DotPlot (tiff)", icon = icon("download")),
                       actionButton(session$ns("download_dotplot_pdf"),
                                      "Dotplot (pdf)", icon = icon("download")),
                       actionButton(session$ns("download_dotplot_data"),
                                      "Dotplot (csv)", icon = icon("download")),
  
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
          filter(clustA %in% input$cluster_selected_dotplot)
      })
      
      # get dotplot
      
      dot_list <- reactive({
        req(data.dotplot.filt(), input$col_low, input$col_high)
        getDotPlot_selInt(data.dotplot.filt(),
                          clust.order = unique(data.dotplot.filt()$clustA),
                          low_color = input$col_low,
                          high_color = input$col_high)
      })
      
      
      # get height size for dotplot
      
      n_rows_dot <- reactive({
        clust_p <- unite(data.dotplot.filt(), col = "clust_p", clustA:clustB)
        length(unique(clust_p$clust_p))
      })
      
  
 
  
      # generate UI plot
      output$gene.dotplot.ui <- renderUI({
        plotOutput(session$ns("gene.dotplot"),
                   height = max(500, 30*n_rows_dot())) %>% withSpinner()
      })
      
      # generate plot
      output$gene.dotplot <- renderPlot({
        req(dot_list())
        plot(dot_list()$p)
      })
      
    rv$n_rows_dot <- n_rows_dot()
    rv$dotplot <- dot_list()$p
    rv$dot_data <- dot_list()$data_dot
    }
    
  })
  
    # Download Dotplot (tiff)
    observeEvent(input$download_dotplot_tiff, {
      validate(need(input$dot_tag, "Please specify file tag!"))
      dir.create(file.path(out_folder(), "gene_verse"), showWarnings = FALSE)
      file <- file.path(out_folder(), "gene_verse", 
                        paste("IntPairs_selected", input$dot_tag, "dotplot.tiff", sep = "_"))
      tiff(file, height = max(500, 30*rv$n_rows_dot))
      plot(rv$dotplot)
      dev.off()
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "l")
    })
    
    # Download dotplot (pdf)
    observeEvent(input$download_dotplot_pdf, {
      validate(need(input$dot_tag, "Please specify file tag!"))
      dir.create(file.path(out_folder(), "gene_verse"), showWarnings = FALSE)
      file <- file.path(out_folder(), "gene_verse", 
                        paste("IntPairs_selected", input$dot_tag, "dotplot.pdf", sep = "_"))
      ggsave(filename = file,
             plot = rv$dotplot,
             device = "pdf", width = 12, height = 20, units = "cm", scale = 2)
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "l")
    })
    
    
    # Download dotplot (csv)
    observeEvent(input$download_dotplot_data, {
      validate(need(input$dot_tag, "Please specify file tag!"))
      dir.create(file.path(out_folder(), "gene_verse"), showWarnings = FALSE)
      file <- file.path(out_folder(), "gene_verse", 
                        paste("IntPairs_selected", input$dot_tag, "dotplot.csv", sep = "_"))
      write.csv(rv$dot_data, file, quote = TRUE, row.names = FALSE)
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "l")
    })


    ####--- Network ---####

    observeEvent(input$gene_table_rows_selected, {
      if(length(input$gene_table_rows_selected) > 0){
        intpair_selected <- reactive({
          as.character(gene.table()$int_pair[input$gene_table_rows_selected])
        })

        output$sel_intpair_text <- renderText({
          paste(intpair_selected(), collapse = "<br>")
        })

        data.filt.net <- reactive({
          d <- rv$gene.filt.data %>%
            filter(int_pair %in% intpair_selected()) %>%
            filter(int.type %in% tolower(input$autocrine_checkbox_net))

        })


        rv$net <- createNetwork(data.filt.net(), input$num_or_weight_radio, input$edge_weight)

        # Plot network
        output$gene.net <- renderVisNetwork({
          validate(
            need(!is.null(input$autocrine_checkbox_net), 'Check at least one interaction type!')
          )

          req(rv$net)
          if(any("circle" %in% rv$net$nodes$shape)){
            # cluster names are numbers -> no background
            visNetwork(rv$net$nodes, rv$net$edges, width = "100%") %>%
              visNodes(font = list(size = 18),
                       scaling = list(min = 10, max = 40)) %>%
              visIgraphLayout(smooth = TRUE)
          } else {
            visNetwork(rv$net$nodes, rv$net$edges, width = "100%") %>%
              visNodes(font = list(size = 18, background = "#ffffff"),
                       scaling = list(min = 10, max = 40)) %>%
              visIgraphLayout(smooth = TRUE)
          }

        })
        
        

        

      }
      })
    
    
    # download network
    observeEvent(input$download_network, {
      validate(need(input$net_tag, "Please specify file tag!"))
      dir.create(file.path(out_folder(), "gene_verse"), showWarnings = FALSE)
      file <- file.path(out_folder(), "gene_verse", 
                        paste("IntPairs_selected", input$net_tag,
                              input$num_or_weight_radio,
                              input$edge_weight,  "network.html", sep = "_"))
      network <- visNetwork(rv$net$nodes, rv$net$edges, width = "100%") %>%
        visNodes(font = list(size = 18, background = "#ffffff"),
                 scaling = list(min = 10, max = 40)) %>%
        visIgraphLayout(smooth = TRUE)
      htmlwidgets::saveWidget(network, file = file, selfcontained = TRUE)
      
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "l")
    })
    
    return(rv)
 
  })
}
    
