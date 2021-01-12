#' table_view UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList fluidRow column downloadButton 
#' @importFrom shinydashboard box

#' @importFrom DT dataTableOutput
#' @importFrom shinycssloaders withSpinner

mod_table_view_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(width = 12,
             h2("Interactions: Table View")
      )
    ),
    fluidRow(
      box(width = 12,
          status = "primary",
          column(width = 8,
                 h4("Here you can see your uploaded cell-cell 
                            interactions, pre-processed  by InterCellar."),
          ),
          column(width=4,
                 downloadButton(ns("download_table_view"), "Download Table"),
          ),
          column(width=12, 
                 br(),
                 DT::dataTableOutput(ns("input_data")) %>% withSpinner()
          )
          
          
          
      )
    )
    
    
            
    
 
  )
}
    
#' table_view Server Function
#'
#' @noRd 
mod_table_view_server <- function(id, input.data) {
  moduleServer(id, function(input, output, session) {
    output$input_data <- DT::renderDataTable({input.data()}, 
                                             options = list(scrollX= TRUE, 
                                                            scrollCollapse = TRUE, 
                                                            processing = FALSE))
    
  })
}
    

 
