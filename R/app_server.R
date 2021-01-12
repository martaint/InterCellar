#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @noRd
app_server <- function( input, output, session ) {
  # List the first level callModules here
    # Upload
    input.data <- mod_upload_server("upload_ui_1")
    # Table view
    mod_table_view_server("table_view_ui_1", input.data)

}
