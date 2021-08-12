#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @noRd
#' @import shiny
app_server <- function( input, output, session ) {
  
  
  if(golem::get_golem_options("reproducible")){
    seed <- 123
  } else {
    seed <- NULL
  }
  
  # Trying to solve an error from biomaRt
  httr::set_config(httr::config(ssl_verifypeer = FALSE))
  
  mod_about_server("about_ui_1")
  
  
  
  rv <- reactiveValues(input.data = list(db1 = NULL,
                                        db2 = NULL,
                                        db3 = NULL), 
                       filt.data = list(db1 = NULL,
                                        db2 = NULL,
                                        db3 = NULL), 
                       
                       gene.table = list(db1 = NULL,
                                        db2 = NULL,
                                        db3 = NULL),
                       genePairs_func_mat = list(db1 = NULL,
                                                 db2 = NULL,
                                                 db3 = NULL),
                       rank.terms = list(db1 = NULL,
                                         db2 = NULL,
                                         db3 = NULL))
    # Upload
    upload.data <- mod_upload_server("upload_ui_1")
    # Upload custom
    upload.data.custom <- mod_upload_custom_server("upload_custom_ui_1")
    
    # generate select widget in sidebar when data are uploaded
    observeEvent(c(upload.data$db_names, upload.data.custom$db_names), {
      db.names <- c(upload.data$db_names, upload.data.custom$db_names)
      # remove NULL elements from the list
      db.names[sapply(db.names, is.null)] <- NULL
      
      # reverse list elements <-> names
      db.list <- as.list(names(db.names))
      names(db.list) <- unlist(as.character(db.names))
      
      output$select_db <- renderUI({
        selectInput("selected_db", label = h4("Active CCC data:"),
                    choices = db.list,
                    multiple = FALSE)
      })
    }, ignoreInit = TRUE)
    
    #### Integrate uploaded data
    observeEvent(c(upload.data$data, upload.data.custom$data), {
      data <- c(upload.data$data, upload.data.custom$data)
      # remove NULL elements from the list
      data[sapply(data, is.null)] <- NULL
      rv$input.data <- data
    }, ignoreInit = TRUE)
    
    
    
    clust.data = list(db1 = NULL,
                      db2 = NULL,
                      db3 = NULL)
    gene.data = list(db1 = NULL,
                     db2 = NULL,
                     db3 = NULL)
      
    observeEvent(input$selected_db, {
      
      # Table view
      output$table_view <- renderUI({
        mod_table_view_ui(paste0("table_view_ui_1",input$selected_db))
      })
      mod_table_view_server(paste0("table_view_ui_1", input$selected_db), 
                            reactive(rv$input.data[[input$selected_db]]))
      
      # Cluster-verse
      output$cluster_verse <- renderUI({
        mod_cluster_verse_ui(paste0("cluster_verse_ui_1",input$selected_db))
      })
      clust.data[[input$selected_db]] <- mod_cluster_verse_server(id = paste0("cluster_verse_ui_1",input$selected_db),
                                             input.data = reactive(rv$input.data[[input$selected_db]]))
      
     
      # Gene-verse
      output$gene_verse <- renderUI({
        mod_gene_verse_ui(paste0("gene_verse_ui_1",input$selected_db))
      })
      gene.data[[input$selected_db]] <- mod_gene_verse_server(paste0("gene_verse_ui_1",input$selected_db),
                                         reactive(rv$input.data[[input$selected_db]]))

      # Get the saved filtered CCCdata from cluster and gene verse and create filt.data for function-verse
      observeEvent(c(clust.data[[input$selected_db]]$filt.data, gene.data[[input$selected_db]]$gene.filt.data), {
        rv$filt.data[[input$selected_db]] <-  dplyr::intersect(clust.data[[input$selected_db]]$filt.data,
                                                               gene.data[[input$selected_db]]$gene.filt.data)
      })


      observeEvent(gene.data[[input$selected_db]]$gene.table, {
        rv$gene.table[[input$selected_db]] <- gene.data[[input$selected_db]]$gene.table
      })

      

      # Function-verse
      output$function_verse <- renderUI({
        req(rv$filt.data)
        mod_function_verse_ui(paste0("function_verse_ui_1",input$selected_db))
      })
      func.data <- mod_function_verse_server(paste0("function_verse_ui_1",input$selected_db),
                                             reactive(rv$filt.data[[input$selected_db]]),
                                             reactive(rv$gene.table[[input$selected_db]]))
      observeEvent(func.data$genePairs_func_mat, {
        rv$genePairs_func_mat[[input$selected_db]] <- func.data$genePairs_func_mat
      })
      observeEvent(func.data$rank.terms, {
        rv$rank.terms[[input$selected_db]] <- func.data$rank.terms
      })

      # Int-pair modules
      output$int_pair_modules <- renderUI({
        mod_int_pair_modules_ui(paste0("int_pair_modules_ui_1",input$selected_db))
      })
      mod_int_pair_modules_server(paste0("int_pair_modules_ui_1",input$selected_db),
                                  reactive(seed),
                                  reactive(input$sidebarmenu),
                                  reactive(rv$filt.data[[input$selected_db]]),
                                  reactive(rv$genePairs_func_mat[[input$selected_db]]),
                                  reactive(rv$gene.table[[input$selected_db]]),
                                  reactive(rv$rank.terms[[input$selected_db]]))





                     
  
     
        
    })
    
    
    
    
   
    

}
