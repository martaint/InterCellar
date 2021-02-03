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

#' @importFrom DT DTOutput
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
                 br(),
                 br(),
                 downloadButton(ns("download_table_view"), "Download Table")
          ),
          
          column(width=12, 
                 br(),
                 DT::DTOutput(ns("input_data")) %>% withSpinner()
          )
          
          
          
      )
    )
    
    
            
    
 
  )
}
    
#' table_view Server Function
#' @importFrom utils write.csv
#' @noRd 
mod_table_view_server <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    observeEvent(data(), {
      output$input_data <- DT::renderDT({data()}, 
                                               options = list(scrollX= TRUE, 
                                                              scrollCollapse = TRUE, 
                                                              processing = FALSE))
    })
    
    # Download Table
    output$download_table_view <- downloadHandler(
      filename = function() {
        "TabView_preprocessed_table.csv"
      },
      content = function(file) {
        write.csv(data(), file, quote = TRUE, row.names = FALSE)
      }
    )
    
    
  })
}
    

 
