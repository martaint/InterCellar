## code to prepare `internal_data` dataset goes here

### Attention!!! usethis::use_data(internal = TRUE) will create /R/sysdata.rda object by 
### overwriting the existing data -> first load the existing datasets!!

#load("./R/sysdata.rda")
library(dplyr)


####--- Pathways from graphite db ---####


library(graphite)

# list of databases available for human and mouse
path.list <- pathwayDatabases()
hsapiens.db <- as.character(path.list$database[path.list$species == "hsapiens"])

getDB <- function(species, database){
    db <- graphite::pathways(species = species, database = database)
    db <- graphite::convertIdentifiers(db, to = "SYMBOL")
    db <- lapply(db, function(x) graphite::nodes(x))
    return(db)
}

hs_biocarta <- getDB(species = "hsapiens", database = "biocarta")
#biocarta: 31.135 sec elapsed

hs_kegg <- getDB(species = "hsapiens", database = "kegg")
#kegg: 91.647 sec elapsed

hs_nci <- getDB(species = "hsapiens", database = "nci")
#nci: 71.819 sec elapsed

hs_panther <- getDB(species = "hsapiens", database = "panther")
#panther: 17.56 sec elapsed

#hs_pathbank <- getDB(species = "hsapiens", database = "pathbank")
# too long

hs_pharmgkb <- getDB(species = "hsapiens", database = "pharmgkb")
#pharmgkb: 10.627 sec elapsed

hs_reactome <- getDB(species = "hsapiens", database = "reactome")
#reactome: 783.533 sec elapsed

#hs_smpdb <- getDB(species = "hsapiens", database = "smpdb")
# too long



##### Saving only datasets that take long to download
usethis::use_data(hs_biocarta,
                  hs_kegg,
                  hs_nci,
                  hs_panther,
                  hs_pharmgkb,
                  hs_reactome,
                  overwrite = TRUE,
                  internal = TRUE,
                  compress = "bzip2")


