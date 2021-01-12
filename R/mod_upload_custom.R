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
#' @importFrom DT dataTableOutput
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
                   tags$li(p(tags$b("clustA"), ": Cluster name or 
                                      number where the expression of ", 
                             tags$em("compA"), "is tested.")),
                   tags$li(p(tags$b("clustB"), ": Cluster name or 
                                      number where the expression of ", 
                             tags$em("compB"), "is tested.")),
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
                                      multiple HGNC symbols are expected, ',' 
                                      separated.")),
                   tags$li(p(tags$b("typeA/typeB"), ": molecular type 
                                      of ", em("compA/compB"), ", either L 
                                      (ligand) or R (receptor).")),
                   tags$li(p(tags$b("pvalue"), ": a statistical p 
                                      value for the interaction.")),
                 ))
             ),
      ),
      column(width = 6,
             h3("Exemplary input table"),
             DT::dataTableOutput(ns("custom_input"))
      )
    ),
    fluidRow(
      hr(),
      column(width = 4, offset = 1,
             fileInput(ns("custom_input_file"), 
                       "Choose Input File (.csv/.tsv/.xlsx)", 
                       multiple = FALSE, 
                       accept = c(".csv", ".tsv", ".xlsx")),
             radioButtons(ns("custom_input_sepInt"), 
                          label = "int_pair separator:",
                          choices = list("_" = "_",
                                         "&" = "&",
                                         ":" = ":")),
      ),
      column(width = 7,
             h3("Uploaded input table"),
             DT::dataTableOutput(ns("custom_table_user")))
    )
             
    
 
  )
}
    
#' upload_custom Server Function
#'
#' @noRd 
mod_upload_custom_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
  })
}
    
## To be copied in the UI
# 
    
## To be copied in the server
# mod_upload_custom_server("upload_custom_ui_1")
 
