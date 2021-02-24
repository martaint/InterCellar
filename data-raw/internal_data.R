## code to prepare `internal_data` dataset goes here

### Attention!!! usethis::use_data(internal = TRUE) will create /R/sysdata.rda object by 
### overwriting the existing data -> first load the existing datasets!!

#load("./R/sysdata.rda")
library(dplyr)


####--- GO annotation ---####
## Gene Ontology annotation from biomaRt 

#####----- Ensembl version

ensembl.version.current <- biomaRt::listEnsembl()
ensembl.version.current <- gsub(" ", "", ensembl.version.current[
    ensembl.version.current$biomart == "genes", "version"])

##----- Download the entire GO annotations from biomaRt
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                         dataset = "hsapiens_gene_ensembl")
genes <- biomaRt::getBM(mart = mart,
                        attributes = c("hgnc_symbol"))
GO.biomart <-  biomaRt::getBM(mart = mart,
                              attributes = c("hgnc_symbol","go_id", "name_1006", 
                                             "namespace_1003", "go_linkage_type"), 
                              values = genes$hgnc_symbol,
                              filters = "hgnc_symbol")
GO.biomart <- GO.biomart %>% 
    dplyr::mutate_all(~ifelse(. %in% c("N/A", "null", ""), NA, .)) %>% 
    na.omit()

# Rename columns 
colnames(GO.biomart) <- c("gene_symbol", "go_id", "go_term", 
                          "domain", "go_linkage_type")

# Name it with version
GO_ensembl_hs_102 <- GO.biomart




####--- Pathways from graphite db ---####


library(graphite)

# list of databases available for human and mouse
path.list <- pathwayDatabases()
hsapiens.db <- as.character(path.list$database[path.list$species == "hsapiens"])

getDB <- function(species, database){
    tictoc::tic(database)
    db <- graphite::pathways(species = species, database = database)
    db <- graphite::convertIdentifiers(db, to = "SYMBOL")
    tictoc::toc()
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
usethis::use_data(hs_reactome,
                  overwrite = TRUE,
                  internal = TRUE,
                  compress = "bzip2")

