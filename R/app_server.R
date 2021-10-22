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
  requireNamespace("httr", quietly = TRUE)
  httr::set_config(httr::config(ssl_verifypeer = FALSE))
  
  
  
  mod_about_server("about_ui_1")
  
  null.list <- list(db1 = NULL,
                    db2 = NULL,
                    db3 = NULL,
                    db1_c = NULL,
                    db2_c = NULL,
                    db3_c = NULL)
  
  rv <- reactiveValues(output_folders_path = null.list,
                       input.data = null.list, 
                       filt.data = null.list, 
                       cluster.colors = null.list,
                       clust_checkbox_selected = null.list,
                       clust_minScore = null.list,
                       clust_maxPval = null.list,
                       cluster.list = null.list,
                       clust.filt.data = null.list,
                       gene.filt.data = null.list,
                       gene.table = null.list,
                       function_table = null.list,
                       nTermsBYdataset = null.list,
                       genePairs_func_mat = null.list,
                       rank.terms = null.list,
                       ipM_vp_selected = null.list,
                       ipM_flow_selected = null.list
                       )
    # Upload
    upload.data <- mod_upload_server("upload_ui_1")
    # Upload custom
    upload.data.custom <- mod_upload_custom_server("upload_custom_ui_1",
                                                   output_folder = reactive(upload.data$output_folder))
    
    # Output folders
    
    observeEvent(c(upload.data$output_folders_path, upload.data.custom$output_folders_path), {

      # make one list of output folders
      output_tags <- c(upload.data$output_folders_path, upload.data.custom$output_folders_path)
      # remove NULL elements from the list
      output_tags[sapply(output_tags, is.null)] <- NULL

      #check that tags are not repeated between supported and custom
      if(any(duplicated(output_tags))){
        shinyalert(text = "It looks like some tags of output folders are not unique! Please re-upload your data after changing repeated tags
                   to avoid overwriting results!",
                   type = "error",
                   showCancelButton = FALSE)
      } else {
        # create list to be passed to server modules
        rv$output_folders_path <- c(upload.data$output_folders_path, upload.data.custom$output_folders_path)
        
      }


    })
    
    # Output tags
    
    
    observeEvent(c(upload.data$output_tags, upload.data.custom$output_tags), {
      
      # make one list of output tags
      rv$output_tags <- c(upload.data$output_tags, upload.data.custom$output_tags)
      
      
      
    })
    
    
    # generate select widget in sidebar when data are uploaded
    observeEvent(c(upload.data$db_names, upload.data.custom$db_names), {
      db.names <- c(upload.data$db_names, upload.data.custom$db_names)
      # remove NULL elements from the list
      db.names[sapply(db.names, is.null)] <- NULL
      
      # reverse list elements <-> names
      db.list <- as.list(names(db.names))
      names(db.list) <- unlist(as.character(db.names))
      
      output$select_db <- renderUI({
        selectInput("selected_db", label = h4("Active CCI data:"),
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
    ip.modules.data = null.list
    
      # Table view
      mod_table_view_server("table_view_ui_1", 
                            reactive(rv$input.data[[input$selected_db]]),
                            out_folder = reactive(rv$output_folders_path[[input$selected_db]]))
      # Cluster-verse
      clust.data <- mod_cluster_verse_server(id = "cluster_verse_ui_1",
                                             reactive(input$sidebarmenu),
                                             input.data = reactive(rv$input.data[[input$selected_db]]),
                                             cluster.list = reactive(rv$cluster.list[[input$selected_db]]),
                                             checkbox_selected = reactive(rv$clust_checkbox_selected[[input$selected_db]]),
                                             minScore = reactive(rv$clust_minScore[[input$selected_db]]),
                                             maxPval = reactive(rv$clust_maxPval[[input$selected_db]]),
                                             out_folder = reactive(rv$output_folders_path[[input$selected_db]]))
      
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
      observeEvent(clust.data$filt.data, {
        req(input$selected_db)
        rv$clust.filt.data[[input$selected_db]] <- clust.data$filt.data
      })
      
      output$debug <- renderUI({
        verbatimTextOutput("debug_text")
      })
      
      output$debug_text <- renderPrint({
        print(gene.data$n_rows_dot)
      })
      
      # Gene-verse
      
      gene.data <- mod_gene_verse_server("gene_verse_ui_1",
                                         reactive(input$sidebarmenu),
                                         reactive(rv$input.data[[input$selected_db]]),
                                         reactive(rv$gene.table[[input$selected_db]]),
                                         out_folder = reactive(rv$output_folders_path[[input$selected_db]]))
      
      observeEvent(gene.data$gene.table_out, {
        req(input$selected_db)
        rv$gene.table[[input$selected_db]] <- gene.data$gene.table_out
      })
      
      observeEvent(gene.data$gene.filt.data, {
        req(input$selected_db)
        rv$gene.filt.data[[input$selected_db]] <- gene.data$gene.filt.data
      })
      
      
      # Get the saved filtered CCIdata from cluster and gene verse and create filt.data for function-verse
      observeEvent(c(clust.data$filt.data, gene.data$gene.filt.data), {
        req(input$selected_db)
        if(!is.null(rv$clust.filt.data[[input$selected_db]]) & !is.null(rv$gene.filt.data[[input$selected_db]])){
          rv$filt.data[[input$selected_db]] <-  dplyr::intersect(rv$clust.filt.data[[input$selected_db]],
                                                                 rv$gene.filt.data[[input$selected_db]])
        } else if(!is.null(rv$clust.filt.data[[input$selected_db]]) & is.null(rv$gene.filt.data[[input$selected_db]])){
          rv$filt.data[[input$selected_db]] <-  rv$clust.filt.data[[input$selected_db]]
        } else if(is.null(rv$clust.filt.data[[input$selected_db]]) & !is.null(rv$gene.filt.data[[input$selected_db]])){
          rv$filt.data[[input$selected_db]] <-  rv$gene.filt.data[[input$selected_db]]
        }
        
      })
      
      
      
      # Function-verse
      func.data <- mod_function_verse_server("function_verse_ui_1",
                                             reactive(rv$filt.data[[input$selected_db]]),
                                             reactive(rv$function_table[[input$selected_db]]),
                                             reactive(rv$nTermsBYdataset[[input$selected_db]]),
                                             reactive(rv$rank.terms[[input$selected_db]]),
                                             out_folder = reactive(rv$output_folders_path[[input$selected_db]]))
      
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
      
  

      # Int-pair modules
      
      
      ip.modules.data <- mod_int_pair_modules_server("int_pair_modules_ui_1",
                                  reactive(seed),
                                  reactive(input$sidebarmenu),
                                  reactive(rv$filt.data[[input$selected_db]]),
                                  reactive(rv$genePairs_func_mat[[input$selected_db]]),
                                  reactive(rv$rank.terms[[input$selected_db]]),
                                  reactive(rv$ipM_vp_selected[[input$selected_db]]),
                                  reactive(rv$ipM_flow_selected[[input$selected_db]]),
                                  out_folder = reactive(rv$output_folders_path[[input$selected_db]]))

     
      observeEvent(ip.modules.data$ipM_vp_selected, {
        req(input$selected_db)
        rv$ipM_vp_selected[[input$selected_db]] <- ip.modules.data$ipM_vp_selected
      })
      observeEvent(ip.modules.data$ipM_flow_selected, {
        req(input$selected_db)
        rv$ipM_flow_selected[[input$selected_db]] <- ip.modules.data$ipM_flow_selected
      })

      

      
      # Multiple conditions
     multi.cond <-  mod_multi_cond_server("multi_cond_ui_1",
                            reactive(input$sidebarmenu),
                            reactive(rv$db.list),
                            reactive(rv$cluster.colors),
                            reactive(rv$filt.data),
                            reactive(rv$genePairs_func_mat),
                            reactive(rv$rank.terms),
                            out_folder = reactive(upload.data$output_folder),
                            output_tags = reactive(rv$output_tags)
      )
      
   
    
   
    

}
