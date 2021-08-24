#' multi_cond UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_multi_cond_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' multi_cond Server Functions
#'
#' @noRd 
mod_multi_cond_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}
    

