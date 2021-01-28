#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @noRd
app_server <- function( input, output, session ) {
  # List the first level callModules here
  rv <- reactiveValues(input.data = NULL, filt.data = NULL)
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
    clust.data <- mod_cluster_verse_server("cluster_verse_ui_1", reactive(rv$input.data))
    observeEvent(clust.data$filt.data, {
      rv$filt.data <- clust.data$filt.data
    })
    
    # Gene-verse
    gene.data <- mod_gene_verse_server("gene_verse_ui_1", reactive(rv$filt.data))
    observeEvent(gene.data$gene.filt.data, {
      rv$filt.data <- gene.data$gene.filt.data
    })
    
    # Function-verse
    func.data <- mod_function_verse_server("function_verse_ui_1", reactive(rv$filt.data), 
                              reactive(gene.data$gene.table))
    
    # Int-pair modules
    mod_int_pair_modules_server("int_pair_modules_ui_1", 
                                reactive(rv$filt.data), 
                                reactive(func.data$genePairs_func_mat),
                                reactive(gene.data$gene.table),
                                reactive(func.data$rank.terms))

}
