# Building a Prod-Ready, Robust Shiny Application.
# 
# README: each step of the dev files is optional, and you don't have to 
# fill every dev scripts before getting started. 
# 01_start.R should be filled at start. 
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
# 
# 
###################################
#### CURRENT FILE: DEV SCRIPT #####
###################################

# Engineering

## Dependencies ----
## Add one line by package you want to add as dependency
usethis::use_package( "shinydashboard" )
usethis::use_package( "shinyFiles" )
usethis::use_package( "shinycssloaders" )
usethis::use_package( "data.table" )
usethis::use_package( "fs" )
usethis::use_package( "dplyr" )
usethis::use_package( "tidyr" )
usethis::use_package( "circlize" )
usethis::use_package( "colourpicker" )
usethis::use_package( "dendextend" )
usethis::use_package( "factoextra" )
usethis::use_package( "ggplot2" )
usethis::use_package( "graphite" )
usethis::use_package( "htmlwidgets" , type = "Suggests")
usethis::use_package( "plotly" )
usethis::use_package( "plyr" )
usethis::use_package( "readxl" )
usethis::use_package( "scales" , type = "Suggests")
usethis::use_package( "shinyFeedback" )
usethis::use_package( "shinyalert" )
usethis::use_package( "tibble" )
usethis::use_package( "umap" )
usethis::use_package( "visNetwork" )
usethis::use_package( "wordcloud2" )
usethis::use_package( "xlsx" )
usethis::use_package( "colorspace" , type = "Suggests")
usethis::use_package( "signal" , type = "Suggests")
usethis::use_package( "igraph" )
usethis::use_package( "ComplexHeatmap")

## Add modules ----
## Create a module infrastructure in R/
golem::add_module( name = "upload" ) # Name of the module
golem::add_module( name = "upload_custom" ) 
golem::add_module( name = "table_view" ) 
golem::add_module( name = "cluster_verse" ) 
golem::add_module( name = "gene_verse" ) 
golem::add_module( name = "function_verse" ) 
golem::add_module( name = "int_pair_modules" ) 

## Add helper functions ----
## Creates ftc_* and utils_*
golem::add_fct( "upload" ) 
golem::add_fct( "upload_custom" ) 
golem::add_fct( "cluster_verse" ) 
golem::add_fct( "gene_verse" ) 
golem::add_fct( "function_verse" ) 
golem::add_fct( "int_pair_modules" ) 

golem::add_utils( "upload" )


## External resources
## Creates .js and .css files at inst/app/www
golem::add_js_file( "script" )
golem::add_js_handler( "handlers" )
golem::add_css_file( "custom" )

## Add internal datasets ----
## If you have data in your package
usethis::use_data_raw( name = "internal_data" ) 

## Tests ----
## Add one line by test you want to create
usethis::use_test( "app" )

# Documentation

## Vignette ----
usethis::use_vignette("InterCellar")
usethis::use_vignette("example_workflow")


## Code coverage ----

# Travis CI
usethis::use_travis() 
usethis::use_travis_badge() 


usethis::use_appveyor()

# You're now set! ----
# go to dev/03_deploy.R
rstudioapi::navigateToFile("dev/03_deploy.R")

