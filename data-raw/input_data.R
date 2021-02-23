## code to prepare `input_data` dataset goes here

source("../R/fct_upload.R")
folder <- "/marta_home/InterCellar-docs/validation_data/Tirosh_melanoma/CPDB_out/"
input_data <- read.CPDBv2(folder)

usethis::use_data(input_data, compress = "bzip2")


