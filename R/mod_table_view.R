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
                            communication, pre-processed  by InterCellar."),
                 br(),
                 br(),
                 actionButton(ns("download_table_view"), "Table (csv)", icon = icon("download"))
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
mod_table_view_server <- function(id, data, out_folder) {
  moduleServer(id, function(input, output, session) {
    
    
    output$input_data <- DT::renderDT({
      d <- data()
      d$clustA <- as.factor(d$clustA)
      d$clustB <- as.factor(d$clustB)
      d
    }, filter = list(position = 'top', clear = FALSE), 
    options = list(scrollX= TRUE, 
                   scrollCollapse = TRUE, processing = FALSE))
    
    
    # Download Table
    observeEvent(input$download_table_view, {
      dir.create(file.path(out_folder(), "table_view"), showWarnings = FALSE)
      file <- file.path(out_folder(), "table_view", "preprocessed_table.csv")
      write.csv(data(), file = file, quote = TRUE, row.names = FALSE)
      shinyalert(text = paste("Saved!", file, sep = "\n"), 
                 type = "success",
                 showCancelButton = FALSE,
                 size = "m")
    })
    

    
    
    
    
 })
}
    

 
