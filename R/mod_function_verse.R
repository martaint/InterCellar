#' function_verse UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList checkboxInput conditionalPanel selectInput 
#' checkboxGroupInput actionButton uiOutput
#' @importFrom shinydashboard valueBoxOutput
#' @importFrom DT dataTableOutput
#' @importFrom plotly plotlyOutput
#' @importFrom shinycssloaders withSpinner
mod_function_verse_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(width = 8,
             h2("Function-verse")),
      valueBoxOutput(ns("tot_functions"))
    ),
    fluidRow(
      box(title = "Annotation",
          width = 12,
          status = "warning",
          solidHeader = TRUE,
          collapsible = TRUE,
          
          column(width = 6,
                 checkboxInput(ns("go_checkbox"),
                               label = h4("Gene Ontology"),
                               value = FALSE),
                 conditionalPanel(
                   condition = "input.go_checkbox",
                   ns = ns,
                   selectInput(ns("select_ensembl"), 
                               label = "Select Ensembl version",
                               choices = list("Ensembl Genes 102" = "102"),
                               multiple = FALSE),
                   column(width = 6,
                          checkboxGroupInput(ns("go_sources_checkbox"),
                                             label = "Ontologies",
                                             choices = list(`GO:Biological Process`= "biological_process",
                                                            `GO:Cellular Component`= "cellular_component",
                                                            `GO:Molecular Function`= "molecular_function"),
                                             selected = c("biological_process", 
                                                          "cellular_component", 
                                                          "molecular_function"),
                                             inline = FALSE)
                   ),
                   column(width = 6,
                          selectInput(ns("go_evidence_exclude"),
                                      label = HTML("<p>Exclude <a href='http://geneontology.org/docs/guide-go-evidence-codes/' target='_blank'>Evidence Code</a></p>"),
                                      choices = list("EXP", "IDA", "IPI", "IMP", 
                                                     "IGI", "IEP", "HTP", "HDA", 
                                                     "HMP", "HGI", "HEP", "IBA", 
                                                     "IBD", "IKR", "IRD", "ISS", 
                                                     "ISO", "ISA", "ISM", "IGC", 
                                                     "RCA", "TAS", "NAS", "IC", 
                                                     "ND", "IEA"),
                                      multiple = TRUE)
                   )
                 )
                 
          ),
          column(width = 6,
                 checkboxInput(ns("pathways_checkbox"),
                               label = h4("Pathways"),
                               value = FALSE),
                 conditionalPanel(
                   condition = "input.pathways_checkbox",
                   ns = ns,
                   checkboxGroupInput(ns("pathways_sources_checkbox"),
                                      label = "Databases",
                                      choices = list(`BioCarta`= "biocarta",
                                                     `KEGG`= "kegg",
                                                     `NCI-Nature`= "nci",
                                                     `PANTHER`= "panther",
                                                     `PathBank` ="pathbank",
                                                     `PharmGKB`="pharmgkb",
                                                     `Reactome`= "reactome",
                                                     `SMPDB`= "smpdb"),
                                      selected = c("biocarta","kegg",
                                                   "nci", "panther","pharmgkb",
                                                   "reactome"),
                                      inline = TRUE)
                 ),
                 column(width = 3, offset=9,
                        actionButton(ns("annotate"), 
                                     label = h4("Annotate!"), 
                                     class="btn-warning")
                 )
          )
      )
    ),
    fluidRow(
      tabBox(
        id = ns('function_verse_tabbox'),
        width = 12,
        tabPanel(h4("Table"),
                 uiOutput(ns("download_funcverse_tab_ui")),
                 br(),
                 br(),
                 DT::dataTableOutput(ns("function_table")) 
        ),
        tabPanel(h4("Barplot"),
                 plotlyOutput(ns("function_bar")) %>% withSpinner()),
        tabPanel(h4("Ranking"),
                 downloadButton(ns("download_rankTab"), "Download Table"),
                 br(),
                 br(),
                 DT::dataTableOutput(ns("function_rank_table"))), 
        tabPanel(h4("Sunburst"),
                 uiOutput(ns("sunburst.text.ui")),
                 uiOutput(ns("sunburst.ui")))
        
        
        
      )
    )
 
  )
}
    
#' function_verse Server Functions
#'
#' @noRd 
#' @importFrom shiny Progress
#' @importFrom xlsx write.xlsx
#' @importFrom dplyr mutate group_by summarise arrange n
#' @importFrom plotly renderPlotly plot_ly layout config
#' @importFrom htmlwidgets JS
mod_function_verse_server <- function(id, filt.data, gene.table){
  moduleServer( id, function(input, output, session){
    
    rv <- reactiveValues(nTermsBYdataset = NULL, 
                         genePairs_func_mat = NULL,
                         rank.terms = NULL)
    #### Annotate!
    
    observeEvent(input$annotate, {
      if(!input$go_checkbox & !input$pathways_checkbox){
        shinyalert(text = "Please select at least one source of annotation!",
                   type = "error",
                   showCancelButton = FALSE)
        source <- input$go_checkbox & input$pathways_checkbox
      } else{
        source <- input$go_checkbox | input$pathways_checkbox
      }
      
      req(source)
        
        
     
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message= "Performing Annotation", value = 0)
      
      # Gene Ontology
      if(input$go_checkbox){
        progress$set(value= 0.2, detail = "GO")
        GO_annotation <- annotateGO(input$select_ensembl, 
                                    input$go_evidence_exclude, 
                                    input$go_sources_checkbox,
                                    filt.data())
        if(!input$pathways_checkbox){
          pathways_annotation <- NULL
        }
      }
      if(input$pathways_checkbox){
        if(input$go_checkbox){
          progress$set(value= 0.5, detail = "Pathways")
        } else{
          progress$set(value= 0.2, detail = "Pathways")
          GO_annotation <- NULL
        }
        
        pathways_annotation <- annotatePathways(species = "hsapiens", 
                                                input$pathways_sources_checkbox, 
                                                filt.data())
      }
      
      
      progress$set(value= 0.8, detail = "Creating Table")
      # Combine GO and pathways
      if(!is.null(GO_annotation) & !is.null(pathways_annotation)){
        data.fun.annot <- reactive({
          combineAnnotations(GO_annotation, pathways_annotation)
        })
        nTermsBYdatasetGO <- getNtermsBYdb(GO_annotation)
        nTermsBYdatasetPath <- getNtermsBYdb(pathways_annotation)
        rv$nTermsBYdataset <- rbind(nTermsBYdatasetGO, nTermsBYdatasetPath)
      } else if(!is.null(GO_annotation) & is.null(pathways_annotation)){
        data.fun.annot <- reactive({
          GO_annotation
        })
        rv$nTermsBYdataset <- getNtermsBYdb(GO_annotation)
      } else if(is.null(GO_annotation) & !is.null(pathways_annotation)){
        data.fun.annot <- reactive({ 
          pathways_annotation
        })
        rv$nTermsBYdataset <- getNtermsBYdb(pathways_annotation)
      } 
      
      
      
      rv$genePairs_func_mat <- buildPairsbyFunctionMatrix(data.fun.annot())
      
      output$download_funcverse_tab_ui <- renderUI({
        req(data.fun.annot())
        downloadButton(session$ns("download_funcTab"), "Download Table")
      })
      
      output$download_funcTab <- downloadHandler(
        filename = function() {
          "Function-verse_table.xlsx"
        },
        content = function(file) {
          write.xlsx(data.fun.annot(), file)
        }
      )
      
      # Plot table
      output$function_table <- DT::renderDataTable({
        if("GO_id" %in% colnames(data.fun.annot())){
          data.fun.annot() %>%
            mutate(GO_id = goLink(GO_id))
        } else{data.fun.annot()}
        
      }, options = list(scrollX= TRUE, 
                        scrollCollapse = TRUE, 
                        processing = FALSE), escape = FALSE)
      
      ## Ranking table of functional terms
      rv$rank.terms <- getRankedTerms(data.fun.annot(), gene.table())
    })   
      
      
      # Plot barplot
      output$function_bar <- renderPlotly({
        req(rv$nTermsBYdataset)
        fig <- plot_ly(rv$nTermsBYdataset, 
                       x = ~source, y = ~n_terms, type = "bar")
        fig <- fig %>% layout(title = "Total number of functional Terms by Source",
                              xaxis = list(title = "Source DB"),
                              yaxis = list(title = "# Terms"))
        fig <- fig %>% config(modeBarButtonsToRemove = c(
                                'sendDataToCloud', 'autoScale2d', 'resetScale2d', 
                                'hoverClosestCartesian', 'hoverCompareCartesian',
                                'zoom2d','pan2d','select2d','lasso2d'
                              ))
        fig
        
      })
      
     
      
      output$function_rank_table <- DT::renderDataTable({
        req(rv$rank.terms)
        rv$rank.terms
      }, options = list(scrollX= TRUE,
                        scrollCollapse = TRUE,
                        processing = FALSE, 
                        columnDefs = list(list(
        targets = 3,
        render = htmlwidgets::JS(
          "function(data, type, row, meta) {",
          "return type === 'display' && data.length > 20 ?",
          "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
          "}")
      ))), escape = FALSE, selection = 'single')
      
      
      output$download_rankTab <- downloadHandler(
        filename = function() {
          "Function-verse_Rank_table.xlsx"
        },
        content = function(file) {
          write.xlsx(rv$rank.terms, file)
        }
      )
      

    
    # Ranking table selection and generation of sunburst plot
    output$no_func_selected <- renderText({
      "Select one functional term from the Ranking to see the cluster enrichment
      in a sunburst plot!"})
    
    output$sunburst.text.ui <- renderUI({
      if(length(input$function_rank_table_rows_selected) == 0){
        textOutput(session$ns("no_func_selected"))
      } else{
        NULL
      }
    })
    
    observeEvent(input$function_rank_table_rows_selected, {
      req(rv$rank.terms)
      
      func_selected <- reactive({
        as.character(rv$rank.terms$functional_term[
          input$function_rank_table_rows_selected])})
      int_p_fun <- reactive({
        int_list <- as.character(rv$rank.terms$int_pair_list[
          input$function_rank_table_rows_selected])
        int_list <- unlist(strsplit(int_list, split=","))
        int_list
      })
      
      sel.data <- filt.data() %>%
        filter(int_pair %in% int_p_fun())
      
      trace1 <- sel.data %>%
        group_by(clustA, int_pair) %>%
        dplyr::summarise(n = length(unique(int_pair)))
      
      trace1 <- trace1 %>%
        group_by(clustA) %>%
        dplyr::summarise(n_tot = n())
      
      trace2 <- sel.data %>%
        group_by(clustA, clustB) %>%
        dplyr::summarise(n = n(), int_pair = paste(int_pair, collapse = ","))
      
      # generate UI
      output$sunburst.ui <- renderUI({
        sidebarLayout(
          sidebarPanel(width = 4,
                       h4("Selected Functional Term:"),
                       textOutput(session$ns("sel_fun_text")),
                       br(),
                       h4("Annotated IntPairs:"),
                       br(),
                       DT::dataTableOutput(session$ns("annot_intp_table")),
          ),
          mainPanel(width = 8,
                    downloadButton(session$ns("download_sunburst"),
                                   "Download Plot"),
                    plotlyOutput(session$ns("sunburst.plot"))
                    
          )
        )
        
      })
      
      output$sel_fun_text <- renderText({
        func_selected()
      })
      output$annot_intp_table <- DT::renderDataTable({
        data.frame(int_pair = int_p_fun())
      }, options = list(scrollX= TRUE,
                        scrollCollapse = TRUE,
                        processing = FALSE), escape = FALSE, selection = 'none')
      
      sunburst.df <- data.frame(
        parents = c(rep(func_selected(), nrow(trace1)), paste0("tr1_",trace2$clustA)), 
        ids = c(paste0("tr1_", trace1$clustA), paste0("tr2_",trace2$clustA, trace2$clustB)),
        labels = c(trace1$clustA, trace2$clustB), 
        values = c(trace1$n_tot, trace2$n),
        text = c(paste(trace1$n_tot, length(int_p_fun()), sep="/"), trace2$int_pair),
        stringsAsFactors = FALSE)
      
      cluster.list <- getClusterNames(filt.data())
      # assign a color to each cluster 
      cluster.colors <- hue_pal(c = 80, l = 80)(length(names(cluster.list)))
      names(cluster.colors) <- names(cluster.list)
      
      
      sunburst.df$parents <- gsub(" ", "<br>", sunburst.df$parents)
      sunburst.df$text <- gsub(",", "<br>", sunburst.df$text)
      sunburst.df$colors <- cluster.colors[sunburst.df$labels]
      
      output$sunburst.plot <- renderPlotly({
        plot_ly(sunburst.df, ids = ~ids, labels = ~labels, 
                parents = ~parents, values = ~values, 
                hovertext = ~text, type = 'sunburst', 
                hoverinfo = "label+text", marker = list(colors= ~colors))
      })
      
      # Download sunburst
      output$download_sunburst <- downloadHandler(
        filename = function() {
          paste0("Function-verse_sunburst_", func_selected(), ".html")},
        content = function(file) {
          fig <- plot_ly(sunburst.df, ids = ~ids, labels = ~labels, 
                         parents = ~parents, values = ~values, 
                         hovertext = ~text, type = 'sunburst', 
                         hoverinfo = "label+text", 
                         marker = list(colors= ~colors))
          htmlwidgets::saveWidget(fig, file = file, selfcontained = TRUE)
        }
      )
      
      
      
      
      
    })
    
    return(rv)
 
  })
}
    
