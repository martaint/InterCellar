#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @noRd
app_server <- function( input, output, session ) {
  # List the first level callModules here
  rv <- reactiveValues(input.data = NULL)
    # Upload
    upload.data <- mod_upload_server("upload_ui_1")
    # Upload custom
    upload.data.custom <- mod_upload_custom_server("upload_custom_ui_1")
    
    # Assign input.data to unique data object
    observeEvent(upload.data$data, {
      rv$input.data <- upload.data$data
    })
    observeEvent(upload.data.custom$data, {
      rv$input.data <- upload.data.custom$data
    })
      
    # Table view
    mod_table_view_server("table_view_ui_1", reactive(rv$input.data))
    
    # Cluster-verse
    mod_cluster_verse_server("cluster_verse_ui_1")

}
