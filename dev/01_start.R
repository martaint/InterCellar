# Building a Prod-Ready, Robust Shiny Application.
# 
# README: each step of the dev files is optional, and you don't have to 
# fill every dev scripts before getting started. 
# 01_start.R should be filled at start. 
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
# 
# 
########################################
#### CURRENT FILE: ON START SCRIPT #####
########################################

## Fill the DESCRIPTION ----
## Add meta data about your application
golem::fill_desc(
  pkg_name = "InterCellar", # The Name of the package containing the App 
  pkg_title = "InterCellar: an R-Shiny app for interactive analysis of cell-cell
  interactions", # The Title of the package containing the App 
  pkg_description = "InterCellar is implemented as an R/Bioconductor Package
  containing a Shiny app that allows users to interactively analyze the results
  of cell-cell interactions for scRNA-seq data.", # The Description of the package containing the App 
  author_first_name = "Marta", # Your First Name
  author_last_name = "Interlandi", # Your Last Name
  author_email = "marta.interlandi01@gmail.com", # Your Email
  repo_url = "https://github.com/martaint/InterCellar" # The URL of the GitHub Repo (optional) 
)     

## Set {golem} options ----
golem::set_golem_options()

## Create Common Files ----
## See ?usethis for more information
usethis::use_mit_license(copyright_holder = "Marta Interlandi")
usethis::use_readme_rmd( open = FALSE )
usethis::use_code_of_conduct()
usethis::use_lifecycle_badge( "Stable" )
usethis::use_news_md( open = FALSE )

## Use git ----
usethis::use_git()

## Vignette
usethis::use_vignette("user_guide")

## Init Testing Infrastructure ----
## Create a template for tests
golem::use_recommended_tests()

## Use Recommended Packages ----
golem::use_recommended_deps()

## Favicon ----
# If you want to change the favicon (default is golem's one)
golem::remove_favicon()
golem::use_favicon() # path = "path/to/ico". Can be an online file. 

## Add helper functions ----
golem::use_utils_ui()
golem::use_utils_server()

# You're now set! ----

# go to dev/02_dev.R
rstudioapi::navigateToFile( "dev/02_dev.R" )

