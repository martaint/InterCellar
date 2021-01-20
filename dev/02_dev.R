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



## Add modules ----
## Create a module infrastructure in R/
golem::add_module( name = "upload" ) # Name of the module
golem::add_module( name = "upload_custom" ) 
golem::add_module( name = "table_view" ) 
golem::add_module( name = "cluster_verse" ) 

## Add helper functions ----
## Creates ftc_* and utils_*
golem::add_fct( "upload" ) 
golem::add_fct( "filters" ) 

golem::add_utils( "upload" )
golem::add_fct( "upload_custom" ) 


## External resources
## Creates .js and .css files at inst/app/www
golem::add_js_file( "script" )
golem::add_js_handler( "handlers" )
golem::add_css_file( "custom" )

## Add internal datasets ----
## If you have data in your package
usethis::use_data_raw( name = "intercellarDB", open = FALSE ) 

## Tests ----
## Add one line by test you want to create
usethis::use_test( "app" )

# Documentation

## Vignette ----
usethis::use_vignette("InterCellar")
devtools::build_vignettes()

## Code coverage ----

# Travis CI
usethis::use_travis() 
usethis::use_travis_badge() 


usethis::use_appveyor()

# You're now set! ----
# go to dev/03_deploy.R
rstudioapi::navigateToFile("dev/03_deploy.R")

