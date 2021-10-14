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
                 status = "success",
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
      box(width = 12,
          title = "Custom CCI data #1",
          status = "primary",
          collapsible = TRUE,
          solidHeader = TRUE,
          column(6,
                 textInput(ns("db1_name"),
                           label = h4("CCI data ID"),
                           placeholder = "my_CCI_data1"),
                 textInput(ns("db1_out_folder"),
                           label = h4("Output folder tag"),
                           placeholder = "my_out_folder1")
          ),
          column(6,
                 radioButtons(ns("db1_sepInt"), 
                              label = "Int_pair separator:",
                              choices = list("_" = "_",
                                             "&" = "&",
                                             ":" = ":"),
                              inline = TRUE),
                 fileInput(ns("db1_file"), 
                           "Choose Input File (.csv/.tsv/.xlsx)", 
                           multiple = FALSE, 
                           accept = c(".csv", ".tsv", ".xlsx"))
                 
                 
          ),
          fluidRow(column(width = 1, offset = 11,
                          actionButton(ns("input_file_button1"), 
                                       label = h4("GO!"),
                                       class = "btn-primary")
          ))
      )), #fluidrow
    fluidRow(
      box(width = 12,
          title = "Custom CCI data #2",
          status = "warning",
          collapsible = TRUE,
          collapsed = TRUE,
          solidHeader = TRUE,
          column(6,
                 textInput(ns("db2_name"),
                           label = h4("CCI data ID"),
                           placeholder = "my_CCI_data2"),
                 textInput(ns("db2_out_folder"),
                           label = h4("Output folder tag"),
                           placeholder = "my_out_folder2")
          ),
          column(6,
                 radioButtons(ns("db2_sepInt"), 
                              label = "Int_pair separator:",
                              choices = list("_" = "_",
                                             "&" = "&",
                                             ":" = ":"),
                              inline = TRUE),
                 fileInput(ns("db2_file"), 
                           "Choose Input File (.csv/.tsv/.xlsx)", 
                           multiple = FALSE, 
                           accept = c(".csv", ".tsv", ".xlsx"))
                 
                 
          ),
          fluidRow(column(width = 1, offset = 11,
                          actionButton(ns("input_file_button2"), 
                                       label = h4("GO!"),
                                       class = "btn-warning")
          ))
      )), # fluidrow
    fluidRow(
      box(width = 12,
          title = "Custom CCI data #3",
          status = "danger",
          collapsible = TRUE,
          collapsed = TRUE,
          solidHeader = TRUE,
          column(6,
                 textInput(ns("db3_name"),
                           label = h4("CCI data ID"),
                           placeholder = "my_CCI_data3"),
                 textInput(ns("db3_out_folder"),
                           label = h4("Output folder tag"),
                           placeholder = "my_out_folder3")
          ),
          column(6,
                 radioButtons(ns("db3_sepInt"), 
                              label = "Int_pair separator:",
                              choices = list("_" = "_",
                                             "&" = "&",
                                             ":" = ":"),
                              inline = TRUE),
                 fileInput(ns("db3_file"), 
                           "Choose Input File (.csv/.tsv/.xlsx)", 
                           multiple = FALSE, 
                           accept = c(".csv", ".tsv", ".xlsx"))
                 
                 
          ),
          fluidRow(column(width = 1, offset = 11,
                          actionButton(ns("input_file_button3"), 
                                       label = h4("GO!"),
                                       class = "btn-danger")
          ))
      )),
             
    
 
  )
}
    
#' upload_custom Server Function
#' @importFrom DT renderDT datatable
#' @importFrom utils read.csv read.table
#' @importFrom tools file_ext 
#' @importFrom readxl read_excel
#' @importFrom shinyFeedback showFeedbackSuccess
#' @noRd 
mod_upload_custom_server <- function(id, output_folder) {
  moduleServer(id, function(input, output, session) {
    rv <- reactiveValues(data = list(db1_c = NULL,
                                     db2_c = NULL,
                                     db3_c = NULL), 
                         db_names = list(db1_c = NULL,
                                         db2_c = NULL,
                                         db3_c = NULL),
                         output_folders_path = list(db1_c = NULL,
                                                    db2_c = NULL,
                                                    db3_c = NULL),
                         output_tags = list(db1_c = NULL,
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
    
    
    
    ##### Custom data #1
    
    observeEvent(input$input_file_button1, {
      
      if(input$db1_name == ""){
        shinyalert(text = "Please specify an ID for your CCI data #1!", type = "error",
                   showCancelButton = FALSE)
      }
      if(input$db1_out_folder == ""){
        shinyalert(text = "Please specify an output folder tag for your CCI data #1!", type = "error",
                   showCancelButton = FALSE)
      }
      if(identical(output_folder(), character(0))){
        shinyalert(text = "Please select an output folder!", 
                   type = "error",
                   showCancelButton = FALSE)
      }
      req(output_folder())
      req(input$db1_name, input$db1_out_folder)
      
      #check that tags are not repeated
      tags <- c(input$db1_out_folder, input$db2_out_folder, input$db3_out_folder)
      tags <- tags[!(tags == "")]
      if(any(duplicated(tags))){
        shinyalert(text = "It looks like tags for output folders are not unique! Please re-upload your data after changing repeated tags 
                   to avoid overwriting results!", 
                   type = "error",
                   showCancelButton = FALSE)
      } else {
        # create directory with tags
        dir.create(file.path(output_folder(), paste0("InterCellar_results_", input$db1_out_folder)))
        shinyalert(text = paste0("Directory created: ", file.path(output_folder(), paste0("InterCellar_results_", input$db1_out_folder))), 
                   type = "success",
                   showCancelButton = FALSE,
                   timer = 4000)
        rv$output_folders_path$db1_c <- file.path(output_folder(), paste0("InterCellar_results_", input$db1_out_folder))
      }
      
      rv$db_names$db1_c <- as.character(input$db1_name)
      rv$output_tags$db1_c <- as.character(input$db1_out_folder)
      
      if(is.null(input$db1_file)){
        shinyalert(text = "Please select a file to upload!", type = "error",
                   showCancelButton = FALSE)
      }
      
      file <- input$db1_file
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
      
      
      # Checks on required columns
      req.columns <- c("int_pair", "typeA", "typeB", "clustA", "clustB", "value")
      if(!all(req.columns %in% colnames(tab))){
        missing.col <- req.columns[!(req.columns %in% colnames(tab))]
        shinyalert(text = paste("Looks like these required columns are missing:", 
                                paste(missing.col, collapse = " "), sep = " "), 
                   type = "error",
                   showCancelButton = FALSE)
      } else if(length(grep(input$db1_sepInt, tab$int_pair)) != nrow(tab)){
        shinyalert(text = "Looks like the chosen separator does not match with 
                   the one found in the uploaded file", 
                   type = "error",
                   showCancelButton = FALSE)
      } else {
        data <- read.customInput(tab, input$db1_sepInt)
        ## Update input.data with ordered L-R interactions
        rv$data$db1_c <- updateInputLR(data)
        shinyalert(text = "Your data was successfully loaded and preprocessed!", 
                   type = "success",
                   showCancelButton = FALSE)
      }
       
    })
    
    ##### Custom data #2
    
    observeEvent(input$input_file_button2, {
      
      if(input$db2_name == ""){
        shinyalert(text = "Please specify an ID for your CCI data #2!", type = "error",
                   showCancelButton = FALSE)
      }
      if(input$db2_out_folder == ""){
        shinyalert(text = "Please specify an output folder tag for your CCI data #2!", type = "error",
                   showCancelButton = FALSE)
      }
      if(identical(output_folder(), character(0))){
        shinyalert(text = "Please select an output folder!", 
                   type = "error",
                   showCancelButton = FALSE)
      }
      req(output_folder())
      req(input$db2_name, input$db2_out_folder)
      
      #check that tags are not repeated
      tags <- c(input$db1_out_folder, input$db2_out_folder, input$db3_out_folder)
      tags <- tags[!(tags == "")]
      if(any(duplicated(tags))){
        shinyalert(text = "It looks like tags for output folders are not unique! Please re-upload your data after changing repeated tags 
                   to avoid overwriting results!", 
                   type = "error",
                   showCancelButton = FALSE)
      } else {
        # create directory with tags
        dir.create(file.path(output_folder(), paste0("InterCellar_results_", input$db2_out_folder)))
        shinyalert(text = paste0("Directory created: ", file.path(output_folder(), paste0("InterCellar_results_", input$db2_out_folder))), 
                   type = "success",
                   showCancelButton = FALSE,
                   timer = 4000)
        rv$output_folders_path$db2_c <- file.path(output_folder(), paste0("InterCellar_results_", input$db2_out_folder))
      }
      
      rv$db_names$db2_c <- as.character(input$db2_name)
      rv$output_tags$db2_c <- as.character(input$db2_out_folder)
      
      if(is.null(input$db2_file)){
        shinyalert(text = "Please select a file to upload!", type = "error",
                   showCancelButton = FALSE)
      }
      
      file <- input$db2_file
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
      
      
      # Checks on required columns
      req.columns <- c("int_pair", "typeA", "typeB", "clustA", "clustB", "value")
      if(!all(req.columns %in% colnames(tab))){
        missing.col <- req.columns[!(req.columns %in% colnames(tab))]
        shinyalert(text = paste("Looks like these required columns are missing:", 
                                paste(missing.col, collapse = " "), sep = " "), 
                   type = "error",
                   showCancelButton = FALSE)
      } else if(length(grep(input$db2_sepInt, tab$int_pair)) != nrow(tab)){
        shinyalert(text = "Looks like the chosen separator does not match with 
                   the one found in the uploaded file", 
                   type = "error",
                   showCancelButton = FALSE)
      } else {
        data <- read.customInput(tab, input$db2_sepInt)
        ## Update input.data with ordered L-R interactions
        rv$data$db2_c <- updateInputLR(data)
        shinyalert(text = "Your data was successfully loaded and preprocessed!", type = "success",
                   showCancelButton = FALSE)
      }
      
    })
    
    ##### Custom data #3
    
    observeEvent(input$input_file_button3, {
      
      if(input$db3_name == ""){
        shinyalert(text = "Please specify an ID for your CCI data #3!", type = "error",
                   showCancelButton = FALSE)
      }
      if(input$db3_out_folder == ""){
        shinyalert(text = "Please specify an output folder tag for your CCI data #3!", type = "error",
                   showCancelButton = FALSE)
      }
      if(identical(output_folder(), character(0))){
        shinyalert(text = "Please select an output folder!", 
                   type = "error",
                   showCancelButton = FALSE)
      }
      req(output_folder())
      req(input$db3_name, input$db3_out_folder)
      
      #check that tags are not repeated
      tags <- c(input$db1_out_folder, input$db2_out_folder, input$db3_out_folder)
      tags <- tags[!(tags == "")]
      if(any(duplicated(tags))){
        shinyalert(text = "It looks like tags for output folders are not unique! Please re-upload your data after changing repeated tags 
                   to avoid overwriting results!", 
                   type = "error",
                   showCancelButton = FALSE)
      } else {
        # create directory with tags
        dir.create(file.path(output_folder(), paste0("InterCellar_results_", input$db3_out_folder)))
        shinyalert(text = paste0("Directory created: ", file.path(output_folder(), paste0("InterCellar_results_", input$db3_out_folder))), 
                   type = "success",
                   showCancelButton = FALSE,
                   timer = 4000)
        rv$output_folders_path$db3_c <- file.path(output_folder(), paste0("InterCellar_results_", input$db3_out_folder))
      }
      
      rv$db_names$db3_c <- as.character(input$db3_name)
      rv$output_tags$db3_c <- as.character(input$db3_out_folder)
      
      
      if(is.null(input$db3_file)){
        shinyalert(text = "Please select a file to upload!", type = "error",
                   showCancelButton = FALSE)
      }
      
      file <- input$db3_file
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
      
      
      # Checks on required columns
      req.columns <- c("int_pair", "typeA", "typeB", "clustA", "clustB", "value")
      if(!all(req.columns %in% colnames(tab))){
        missing.col <- req.columns[!(req.columns %in% colnames(tab))]
        shinyalert(text = paste("Looks like these required columns are missing:", 
                                paste(missing.col, collapse = " "), sep = " "), 
                   type = "error",
                   showCancelButton = FALSE)
      } else if(length(grep(input$db3_sepInt, tab$int_pair)) != nrow(tab)){
        shinyalert(text = "Looks like the chosen separator does not match with 
                   the one found in the uploaded file", 
                   type = "error",
                   showCancelButton = FALSE)
      } else {
        data <- read.customInput(tab, input$db3_sepInt)
        ## Update input.data with ordered L-R interactions
        rv$data$db3_c <- updateInputLR(data)
        shinyalert(text = "Your data was successfully loaded and preprocessed!", type = "success",
                   showCancelButton = FALSE)
      }
      
    })
    
    
    return(rv)
  })
}
    

 
