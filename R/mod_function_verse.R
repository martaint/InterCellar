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
#' @importFrom DT DTOutput
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
                               value = TRUE),
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
                               value = TRUE),
                 conditionalPanel(
                   condition = "input.pathways_checkbox",
                   ns = ns,
                   checkboxGroupInput(ns("pathways_sources_checkbox"),
                                      label = "Databases",
                                      choices = list(`BioCarta`= "biocarta",
                                                     `KEGG`= "kegg",
                                                     `NCI-Nature`= "nci",
                                                     `PANTHER`= "panther",
                                                     `PharmGKB`="pharmgkb",
                                                     `Reactome`= "reactome"
                                                     ),
                                      selected = c("biocarta","kegg",
                                                   "nci", "panther","pharmgkb",
                                                   "reactome"),
                                      inline = TRUE)
                 ),
                 column(width = 3, offset=9,
                        actionButton(ns("annotate"), 
                                     label = h4("Annotate!"), 
                                     class="btn-warning")
                        #uiOutput(ns("action_btn"))
                        
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
                 DT::DTOutput(ns("function_table")) 
        ),
        tabPanel(h4("Barplot"),
                 plotlyOutput(ns("function_bar")) %>% withSpinner()),
        tabPanel(h4("Ranking"),
                 downloadButton(ns("download_rankTab"), "Download Table"),
                 br(),
                 br(),
                 DT::DTOutput(ns("function_rank_table"))), 
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
#' @importFrom utils write.csv
#' @importFrom dplyr mutate group_by summarise arrange n
#' @importFrom plotly renderPlotly plot_ly layout config
#' @importFrom htmlwidgets JS
#' @importFrom DT renderDT DTOutput
#' @importFrom shinyalert shinyalert
mod_function_verse_server <- function(id, filt.data, function_table, nTermsBYdataset, rank.terms){
  moduleServer( id, function(input, output, session){
    
    rv <- reactiveValues(function_table_out = NULL,
                         nTermsBYdataset_out = NULL, 
                         genePairs_func_mat_out = NULL,
                         rank.terms_out = NULL)
    
    
    # output$action_btn <- renderUI({
    #   actionButton(session$ns("annotate"), 
    #                label = h4("Annotate!"), 
    #                class="btn-warning")
    # })
      
    
    #### Annotate!
    # check that at least one source is selected
    ann_source <- reactive({input$go_checkbox | input$pathways_checkbox})
    observeEvent(input$annotate, {
      #req(ann_source())
      if(!ann_source()){
        shinyalert(text = "Please select at least one source of annotation!",
                   type = "error",
                   showCancelButton = FALSE)
      }
    })
    
    
    
    
    data.fun.annot <- eventReactive(input$annotate, {
      req(ann_source())
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message= "Performing Annotation", value = 0)
      
      
      # Gene Ontology
      if(input$go_checkbox){
        progress$set(value= 0.2, detail = "GO")
        GO_annotation <- suppressWarnings(annotateGO(input$select_ensembl, 
                                                     input$go_evidence_exclude, 
                                                     input$go_sources_checkbox,
                                                     filt.data()))
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
        
        pathways_annotation <- annotatePathways(input$pathways_sources_checkbox,
                                                filt.data())
      }
      
      
      progress$set(value= 0.8, detail = "Creating Table")
      # Combine GO and pathways
      if(!is.null(GO_annotation) & !is.null(pathways_annotation)){
        data.fun.annot <- combineAnnotations(GO_annotation, pathways_annotation)
        
        nTermsBYdatasetGO <- getNtermsBYdb(GO_annotation)
        nTermsBYdatasetPath <- getNtermsBYdb(pathways_annotation)
        rv$nTermsBYdataset_out <- rbind(nTermsBYdatasetGO, nTermsBYdatasetPath)
      } else if(!is.null(GO_annotation) & is.null(pathways_annotation)){
        data.fun.annot <- GO_annotation
        rv$nTermsBYdataset_out <- getNtermsBYdb(GO_annotation)
      } else if(is.null(GO_annotation) & !is.null(pathways_annotation)){
        data.fun.annot <- pathways_annotation
        rv$nTermsBYdataset_out <- getNtermsBYdb(pathways_annotation)
      }
      return(data.fun.annot)
    }, ignoreInit = TRUE)
    
    observeEvent(data.fun.annot(), {
      rv$function_table_out <- data.fun.annot()
    })
    
    observeEvent(data.fun.annot(), {
      rv$genePairs_func_mat_out <- buildPairsbyFunctionMatrix(data.fun.annot())
      # Check how many int-pairs could not be annotated
      num_notAnn <- sum(!(unique(filt.data()$int_pair) %in%
                            rownames(rv$genePairs_func_mat_out)))
      if(num_notAnn > 0){
        shinyalert(text = paste0("Warning! With the current choice of functional
                                 databases, ", num_notAnn, " int-pairs could not
                                 be annotated. They will be excluded from further
                                 analysis."),
                   type = "warning",
                   showCancelButton = FALSE)
      }
      # Ranking table of functional terms
      rv$rank.terms_out <- getRankedTerms(data.fun.annot())
    })


      
      
      
      
      output$download_funcverse_tab_ui <- renderUI({
        req(data.fun.annot())
        downloadButton(session$ns("download_funcTab"), "Download Table")
      })
      
      output$download_funcTab <- downloadHandler(
        filename = function() {
          "Function-verse_table.csv"
        },
        content = function(file) {
          write.csv(data.fun.annot(), quote = TRUE, file)
        }
      )
      
      
      
      # Plot table
      output$function_table <- DT::renderDT({
        req(function_table())
        dt <- function_table()
        
        if("GO_id" %in% colnames(dt)){
          dt %>%
            mutate(GO_id = goLink(GO_id))
          } else{dt}


      }, filter = list(position = 'top', clear = FALSE),
      options = list(scrollX= TRUE,
                     scrollCollapse = TRUE,
                     processing = FALSE), escape = FALSE)

      
      
      
    
    
    
    
      # Plot barplot
      output$function_bar <- renderPlotly({
        req(nTermsBYdataset())
        fig <- plot_ly(nTermsBYdataset(),
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



      output$function_rank_table <- DT::renderDT({
        req(rank.terms())
        rank.terms()
      }, filter = list(position = 'top', clear = FALSE),
      options = list(scrollX= TRUE,
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
          "Function-verse_Rank_table.csv"
        },
        content = function(file) {
          write.csv(rank.terms(), file, quote = TRUE)
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
      req(rank.terms())

      func_selected <- reactive({
        as.character(rank.terms()$functional_term[
          input$function_rank_table_rows_selected])})
      int_p_fun <- reactive({
        int_list <- as.character(rank.terms()$int_pair_list[
          input$function_rank_table_rows_selected])
        int_list <- unlist(strsplit(int_list, split=","))
        int_list
      })

      sel.data <- filt.data() %>%
        filter(int_pair %in% int_p_fun())



      # generate UI
      output$sunburst.ui <- renderUI({
        if(length(input$function_rank_table_rows_selected) == 0){
          NULL
        } else{
          sidebarLayout(
            sidebarPanel(width = 4,
                         h4("Selected Functional Term:"),
                         textOutput(session$ns("sel_fun_text")),
                         br(),
                         h4("Annotated IntPairs:"),
                         br(),
                         DT::DTOutput(session$ns("annot_intp_table")),
            ),
            mainPanel(width = 8,
                      downloadButton(session$ns("download_sunburst"),
                                     "Download Plot"),
                      plotlyOutput(session$ns("sunburst.plot")) %>% withSpinner()

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

      cluster.list <- getClusterNames(filt.data())
      # assign a color to each cluster
      cluster.colors <- hue_pal(c = 80, l = 80)(length(names(cluster.list)))
      names(cluster.colors) <- names(cluster.list)

      output$sunburst.plot <- renderPlotly({
        getSunburst(sel.data, func_selected(), int_p_fun(), cluster.colors)
      })

      # Download sunburst
      output$download_sunburst <- downloadHandler(
        filename = function() {
          paste0("Function-verse_sunburst_", func_selected(), ".html")},
        content = function(file) {
          fig <- getSunburst(sel.data, func_selected(), int_p_fun(),
                             cluster.colors)
          htmlwidgets::saveWidget(fig, file = file, selfcontained = TRUE)
        }
      )





    })
    
    return(rv)
 
  })
}
    
