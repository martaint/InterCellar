## code to prepare `internal_data` dataset goes here

### Attention!!! usethis::use_data() will create /R/sysdata.rda object by 
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
    db <- graphite::pathways(species = species, database = database)
    db <- graphite::convertIdentifiers(db, to = "SYMBOL")
    return(db)
}

hs_biocarta <- getDB(species = "hsapiens", database = "biocarta")
hs_kegg <- getDB(species = "hsapiens", database = "kegg")
hs_nci <- getDB(species = "hsapiens", database = "nci")
hs_panther <- getDB(species = "hsapiens", database = "panther")
hs_pathbank <- getDB(species = "hsapiens", database = "pathbank")
hs_pharmgkb <- getDB(species = "hsapiens", database = "pharmgkb")
hs_reactome <- getDB(species = "hsapiens", database = "reactome")
hs_smpdb <- getDB(species = "hsapiens", database = "smpdb")




##### Saving all data together
usethis::use_data(GO_ensembl_hs_102,
                  hs_biocarta,
                  hs_kegg, 
                  hs_nci,
                  hs_panther,
                  #hs_pathbank,
                  hs_pharmgkb,
                  hs_reactome,
                  #hs_smpdb,
                  overwrite = TRUE,
                  internal = TRUE,
                  compress = "bzip2")

