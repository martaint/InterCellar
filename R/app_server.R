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
  
  null.list <- list(db1 = NULL,
                    db2 = NULL,
                    db3 = NULL,
                    db1_c = NULL,
                    db2_c = NULL,
                    db3_c = NULL)
  
  rv <- reactiveValues(input.data = null.list, 
                       filt.data = null.list, 
                       cluster.colors = null.list,
                       clust_checkbox_selected = null.list,
                       clust_minScore = null.list,
                       clust_maxPval = null.list,
                       cluster.list = null.list,
                       function_table = null.list,
                       nTermsBYdataset = null.list,
                       genePairs_func_mat = null.list,
                       rank.terms = null.list
                       )
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
      
      
      rv$db.list <- db.list
      
    }, ignoreInit = TRUE)
    
    #### Integrate uploaded data
    observeEvent(c(upload.data$data, upload.data.custom$data), {
      data <- c(upload.data$data, upload.data.custom$data)
      # remove NULL elements from the list
      #data[sapply(data, is.null)] <- NULL
      rv$input.data <- data
      
      
      # Check that filt.data objects have been generated in universes, if not -> input.data
      for(l in seq_along(rv$filt.data)){
        if(is.null(rv$filt.data[[l]]) & !is.null(rv$input.data[[l]])){
          rv$filt.data[[l]] <- rv$input.data[[l]]
        }
      }
      
      ### Assign colors to clusters of each input data
      
      for(l in seq_along(rv$input.data)){
        if(!is.null(rv$input.data[[l]])){
          rv$cluster.colors[[l]] <- getClusterColors(rv$input.data[[l]]) 
        }
      }
      
      ### Assign cluster-verse initial filters
      
      for(l in seq_along(rv$input.data)){
        if(!is.null(rv$input.data[[l]])){
          cluster.list <- getClusterNames(rv$input.data[[l]])
          rv$clust_checkbox_selected[[l]] <- names(cluster.list)
          rv$clust_minScore[[l]] <- min(rv$input.data[[l]]$score) 
          rv$clust_maxPval[[l]] <- 0.05
          rv$cluster.list[[l]] <- cluster.list
        }
      }
      
    }, ignoreInit = TRUE)
    
    
    
    # initialize lists to hold results of modules
    clust.data = null.list
    gene.data = null.list
    func.data = null.list
      
    
      # Table view
      mod_table_view_server("table_view_ui_1", 
                            reactive(rv$input.data[[input$selected_db]]))
      # Cluster-verse
      clust.data <- mod_cluster_verse_server(id = "cluster_verse_ui_1",
                                             reactive(input$sidebarmenu),
                                             input.data = reactive(rv$input.data[[input$selected_db]]),
                                             cluster.list = reactive(rv$cluster.list[[input$selected_db]]),
                                             checkbox_selected = reactive(rv$clust_checkbox_selected[[input$selected_db]]),
                                             minScore = reactive(rv$clust_minScore[[input$selected_db]]),
                                             maxPval = reactive(rv$clust_maxPval[[input$selected_db]]))
      
      observeEvent(clust.data$checkbox_selected_out, {
        req(input$selected_db)
        rv$clust_checkbox_selected[[input$selected_db]] <- clust.data$checkbox_selected_out
      })
      
      observeEvent(clust.data$minScore_out, {
        req(input$selected_db)
        rv$clust_minScore[[input$selected_db]] <- clust.data$minScore_out
      })
      observeEvent(clust.data$maxPval_out, {
        req(input$selected_db)
        rv$clust_maxPval[[input$selected_db]] <- clust.data$maxPval_out
      })
      
      output$debug <- renderUI({
        verbatimTextOutput("debug_text")
      })
      
      output$debug_text <- renderPrint({
        print(rv$cluster.list)
      })
      
      
      
      
      # Function-verse
      func.data <- mod_function_verse_server("function_verse_ui_1",
                                             reactive(rv$filt.data[[input$selected_db]]),
                                             reactive(rv$function_table[[input$selected_db]]),
                                             reactive(rv$nTermsBYdataset[[input$selected_db]]),
                                             reactive(rv$genePairs_func_mat[[input$selected_db]]),
                                             reactive(rv$rank.terms[[input$selected_db]]))
      
      observeEvent(func.data$function_table_out, {
        req(input$selected_db)
        rv$function_table[[input$selected_db]] <- func.data$function_table_out
      })
      observeEvent(func.data$nTermsBYdataset_out, {
        req(input$selected_db)
        rv$nTermsBYdataset[[input$selected_db]] <- func.data$nTermsBYdataset_out
      })
      observeEvent(func.data$genePairs_func_mat_out, {
        req(input$selected_db)
        rv$genePairs_func_mat[[input$selected_db]] <- func.data$genePairs_func_mat_out
      })
      observeEvent(func.data$rank.terms_out, {
        req(input$selected_db)
        rv$rank.terms[[input$selected_db]] <- func.data$rank.terms_out
      })
      
      
    
    #observeEvent(input$selected_db, {
      
      #input_selected_db <- isolate({input$selected_db})
      
      # # Table view
      # output$table_view <- renderUI({
      #   mod_table_view_ui(paste0("table_view_ui_1", input_selected_db))
      # })
      # mod_table_view_server(paste0("table_view_ui_1", input_selected_db), 
      #                       reactive(rv$input.data[[input_selected_db]]))
      # Table view
      # output$table_view <- renderUI({
      #   mod_table_view_ui("table_view_ui_1")
      # })
    
    
      
      
      # output$cluster_verse <- renderUI({
      #   req(input$selected_db)
      #   mod_cluster_verse_ui("cluster_verse_ui_1")
      # })
      

      # # Gene-verse
      # output$gene_verse <- renderUI({
      #   mod_gene_verse_ui(paste0("gene_verse_ui_1",input_selected_db))
      # })
      # gene.data[[input_selected_db]] <- mod_gene_verse_server(paste0("gene_verse_ui_1",input_selected_db),
      #                                    reactive(rv$input.data[[input_selected_db]]))
      # 
      # # Get the saved filtered CCCdata from cluster and gene verse and create filt.data for function-verse
      # # observeEvent(c(clust.data[[input_selected_db]]$filt.data, gene.data[[input_selected_db]]$gene.filt.data), {
      # #   rv$filt.data[[input_selected_db]] <-  dplyr::intersect(clust.data[[input_selected_db]]$filt.data,
      # #                                                          gene.data[[input_selected_db]]$gene.filt.data)
      # # })
      # 
      # 
      # 
      # 
      # 
      # #output$function_verse <- renderUI({
      # #  mod_function_verse_ui("function_verse_ui_1")
      # #})
      # 
      # Saving each func.data in separate objects to be retrieved by the multi-condition comparison
      
      # 
      # 
      # 
      # # Int-pair modules
      # output$int_pair_modules <- renderUI({
      #   mod_int_pair_modules_ui(paste0("int_pair_modules_ui_1",input_selected_db))
      # })
      # # int-pair modules has no returned values until now
      # mod_int_pair_modules_server(paste0("int_pair_modules_ui_1",input_selected_db),
      #                             reactive(seed),
      #                             reactive(input$sidebarmenu),
      #                             reactive(rv$filt.data[[input_selected_db]]),
      #                             reactive(rv$genePairs_func_mat[[input_selected_db]]),
      #                             reactive(rv$rank.terms[[input_selected_db]]))
      # 
     
      

      

      
      
      
      
   # })
    
    
    
    
    
    
    observeEvent(rv$db.list, {
      
      # Generate multi condition UI and server
      # Multiple conditions
      output$multi_cond <- renderUI({
        mod_multi_cond_ui("multi_cond_ui_1")
      })
      
      
      mod_multi_cond_server("multi_cond_ui_1",
                            reactive(input$sidebarmenu),
                            reactive(rv$db.list),
                            reactive(rv$cluster.colors),
                            reactive(rv$filt.data),
                            reactive(rv$genePairs_func_mat),
                            reactive(rv$rank.terms)
      )
      
    })
    
    
   
    

}
