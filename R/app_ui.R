#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @importFrom shinydashboard dashboardPage dashboardHeader dropdownMenu 
#' notificationItem dashboardSidebar sidebarMenu menuItem menuSubItem 
#' dashboardBody tabItems tabItem
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # List the first level UI elements here 
    dashboardPage(skin = "black",
                  dashboardHeader(title = span(img(src = "www/img_header.png", width = "50px"), "InterCellar"),
                                  dropdownMenu(type="tasks",
                                               icon=icon("github"),
                                               badgeStatus=NULL,
                                               headerText="Check out",
                                               notificationItem(
                                                 text="GitHub Repo",
                                                 icon=icon("code-branch"), 
                                                 status="primary",
                                                 href="https://github.com/martaint/InterCellar"
                                               )
                                  )
                                  ),
                  dashboardSidebar(sidebarMenu(id = "sidebarmenu",
                                               menuItem("About", tabName = "about", icon = icon("rocket")),
                                               menuItem("1. Data",
                                                        menuSubItem("Upload", tabName = "upload", icon = icon("file-import")),
                                                        menuSubItem("Table view", tabName = "data_table", icon = icon("table")),
                                                        icon = icon("database")),
                                               menuItem("2. Universes", icon = icon("meteor"),
                                                        menuSubItem("Cluster-verse", tabName = "cluster-verse",
                                                                    icon = icon("project-diagram")),
                                                        menuSubItem("Gene-verse", tabName = "gene-verse", icon = icon("dna")),
                                                        menuSubItem("Function-verse", tabName = "function-verse", icon = icon("comments"))),
                                               
                                               menuItem("3. Analyze", tabName = "analyze", icon = icon("arrows-alt"),
                                                        menuSubItem("Int-Pair Modules", tabName = "ipModules",
                                                                    icon = icon("handshake")))
                  )),
                  dashboardBody(
                    tabItems(
                      tabItem(tabName = "about",
                              h2("InterCellar: functional visualization of cellular interactions")
                      ),
                      tabItem(tabName = "upload",
                              mod_upload_ui("upload_ui_1")
                      ),
                      tabItem(tabName = "data_table",
                              mod_table_view_ui("table_view_ui_1")
                      ),
                      tabItem(tabName = "cluster-verse",
                              mod_cluster_verse_ui("cluster_verse_ui_1")
                      ),
                      tabItem(tabName = "gene-verse",
                              mod_gene_verse_ui("gene_verse_ui_1")
                      ),
                      tabItem(tabName = "function-verse",
                              mod_function_verse_ui("function_verse_ui_1")
                      ),
                      tabItem(tabName = "ipModules",
                              mod_int_pair_modules_ui("int_pair_modules_ui_1")
                      )
                    )
                  )
    )
  )#
}#

#' Add external Resources to the Application
#' 
#' This function is internally used to add external 
#' resources inside the Shiny application. 
#' 
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @importFrom shinyalert useShinyalert
#' @importFrom shinyFeedback useShinyFeedback
#' @noRd
golem_add_external_resources <- function(){
  
  add_resource_path(
    'www', app_sys('app/www')
  )
  add_resource_path(
    'extdata', app_sys('app/extdata')
  )
 
  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'InterCellar'
    ),
    # Add here other external resources
    
    shinyalert::useShinyalert(),
    shinyFeedback::useShinyFeedback()
  )
}

