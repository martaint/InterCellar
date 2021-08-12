#' upload_custom UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList fluidRow column fileInput radioButtons
#' @importFrom shinydashboard box
#' @importFrom DT DTOutput
#' 
mod_upload_custom_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(width = 5,
             box(width = 12,
                 title = "Required columns",
                 status = "danger",
                 collapsible = TRUE,
                 solidHeader = TRUE,
                 div(tags$ul(
                   tags$li(p(tags$b("int_pair"), ": Interaction pair 
                                      representing the names of two components 
                                      that participate in the interaction 
                                      (e.g. ", tags$em("compA_compB"), "). 
                                      Valid separators: '_', '&', ':'.")),
                   tags$li(p(tags$b("clustA/clustB"), ": Cluster name or 
                                      number where the expression of ", 
                             tags$em("compA/compB"), "is tested.")),
                   tags$li(p(tags$b("typeA/typeB"), ": molecular type 
                                      of ", em("compA/compB"), ", either L 
                                      (ligand) or R (receptor).")),
                   tags$li(p(tags$b("value"), ": A numeric value 
                                      representing a ", tags$em("score"), 
                             "for the interaction (e.g. average 
                                      expression of ", em("compA_compB"), 
                             "over clustA and clustB)."))
                 ))
             ),
             box(width = 12,
                 title = "Additional columns",
                 status = "primary",
                 collapsible = TRUE,
                 collapsed = TRUE,
                 solidHeader = TRUE,
                 div(tags$ul(
                   tags$li(p(tags$b("geneA/geneB"), ": HGNC symbol of 
                                      ", em("compA/compB"), ". When ", 
                             em("compA/compB"), "is a complex, 
                             multiple HGNC symbols are expected, comma separated.
                             If these columns are NOT provided, InterCellar 
                             will use ", em("compA/compB"), " from column ", 
                             tags$b("int_pair"), ".")),
                   tags$li(p(tags$b("p_value"), ": a statistical p 
                                      value for the interaction.")),
                 ))
             ),
      ),
      column(width = 7,
             h3("Example input table"),
             DT::DTOutput(ns("custom_input"))
      )
    ),
    fluidRow(
      hr(),
      column(width = 12, 
             sidebarLayout(
               sidebarPanel(width = 5,
                            
                            fileInput(ns("custom_input_file"), 
                                      "Choose Input File (.csv/.tsv/.xlsx)", 
                                      multiple = FALSE, 
                                      accept = c(".csv", ".tsv", ".xlsx")),
                            radioButtons(ns("custom_input_sepInt"), 
                                         label = "int_pair separator:",
                                         choices = list("_" = "_",
                                                        "&" = "&",
                                                        ":" = ":")),
                            
                            fluidRow(
                              column(width = 2, offset = 10,
                                     actionButton(ns("custom_button"), 
                                                  label = tags$b("Upload"))
                              )
                            )
                            
                            
               ),
               mainPanel(width = 7,  
                         h3("Uploaded input table"),
                         DT::DTOutput(ns("custom_table_user"))
               )
             )
             )
      
      
    )
             
    
 
  )
}
    
#' upload_custom Server Function
#' @importFrom DT renderDT datatable
#' @importFrom utils read.csv read.table
#' @importFrom tools file_ext 
#' @importFrom readxl read_excel
#' @importFrom shinyFeedback showFeedbackSuccess
#' @noRd 
mod_upload_custom_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    rv <- reactiveValues(data = list(db1_c = NULL,
                                     db2_c = NULL,
                                     db3_c = NULL), 
                         db_names = list(db1_c = NULL,
                                         db2_c = NULL,
                                         db3_c = NULL))
    
    output$custom_input <- DT::renderDT({
      custom_tab <- read.csv(app_sys("app", "extdata", "custom_input.csv"), 
                             header = TRUE)
      colnames(custom_tab) <- c("<span style='color:red'>int_pair</span>",
                                "<span style='color:blue'>geneA</span>",
                                "<span style='color:blue'>geneB</span>",
                                "<span style='color:red'>typeA</span>",
                                "<span style='color:red'>typeB</span>",
                                "<span style='color:red'>clustA</span>",
                                "<span style='color:red'>clustB</span>",
                                "<span style='color:red'>value</span>",
                                "<span style='color:blue'>p_value</span>")
      DT::datatable(custom_tab, escape = FALSE, 
                    options = list(scrollX= TRUE, scrollCollapse = TRUE))
      })
    
    observeEvent(input$custom_input_file, {
      shinyFeedback::showFeedbackSuccess(
        inputId = "custom_input_file", 
        text = "Great! Choose separator and click Upload!"
      )
    })
    
    observeEvent(input$custom_button, {
      if(is.null(input$custom_input_file)){
        shinyalert(text = "Please select a file to upload!", type = "error",
                   showCancelButton = FALSE)
      }
      
      file <- input$custom_input_file
      req(file)
      ext <- tools::file_ext(file$datapath)
      
      validate(need(ext %in% c("csv", "tsv", "xlsx"), "Please choose a file
                    with the required extension (.csv/.tsv/.xlsx)."))
      switch (ext,
              csv = {tryCatch({
                tab <- read.csv(file$datapath, header = TRUE)
                },
                error = function(e) {
                  # return a safeError if a parsing error occurs
                  stop(safeError("Error reading input file"))
                }
              )},
              tsv = {tryCatch({
                tab <- read.table(file$datapath, sep = "\t", 
                                  header = TRUE)
              },
              error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError("Error reading input file"))
              }
              )},
              xlsx = {tryCatch({
                tab <- read_excel(file$datapath, col_names = TRUE)
              },
              error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError("Error reading input file"))
              }
              )}
      )
      
      output$custom_table_user <- DT::renderDT({
        req(tab)
        
        tab
      }, options = list(scrollX= TRUE, scrollCollapse = TRUE, pageLength = 5))
      
      # Checks on required columns
      req.columns <- c("int_pair", "typeA", "typeB", "clustA", "clustB", "value")
      if(!all(req.columns %in% colnames(tab))){
        missing.col <- req.columns[!(req.columns %in% colnames(tab))]
        shinyalert(text = paste("Looks like these required columns are missing:", 
                                paste(missing.col, collapse = " "), sep = " "), 
                   type = "error",
                   showCancelButton = FALSE)
      } else if(length(grep(input$custom_input_sepInt, tab$int_pair)) != nrow(tab)){
        shinyalert(text = "Looks like the chosen separator does not match with 
                   the one found in the uploaded file", 
                   type = "error",
                   showCancelButton = FALSE)
      } else {
        data <- read.customInput(tab, input$custom_input_sepInt)
        ## Update input.data with ordered L-R interactions
        rv$data$db1_c <- updateInputLR(data)
        shinyalert(text = "Your data was successfully loaded and preprocessed! 
             Check it out at table view!", type = "success",
                   showCancelButton = FALSE)
      }
       
    })
    
    
    return(rv)
  })
}
    

 
