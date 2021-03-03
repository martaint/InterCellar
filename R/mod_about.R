#' about UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList imageOutput
mod_about_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(width = 6,
             img(src = "www/About_shiny.png", width = "650px")
             ),
      column(width = 6,
             h2("InterCellar: interactive exploration of cellular interactions")
      )
    )
    
    
 
  )
}
    
#' about Server Functions
#'
#' @noRd 
#' @importFrom shiny renderImage
mod_about_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    
    # output$workflow_img <- renderImage({
    #   list(src = "www/About_shiny.png",
    #        contentType = "image/png",
    #        width = 400,
    #        height = 500,
    #        alt = "InterCellar workflow")
    # }, deleteFile = FALSE)
  })
}
    
