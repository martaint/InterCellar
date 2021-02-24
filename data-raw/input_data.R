## code to prepare `input.data` dataset goes here

source("../R/fct_upload.R")
folder <- "/marta_home/InterCellar-docs/validation_data/Tirosh_melanoma/CPDB_out/"
input.data <- read.CPDBv2(folder)

usethis::use_data(input.data, compress = "bzip2")


