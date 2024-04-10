rm(list = ls())

# Set Reproducibility -----------------------------------------------------
set.seed(1234)



# Load Packages -----------------------------------------------------------
packs <- c('ggplot2', 'jsonlite')
lapply(packs, function(pack) {
  if (!requireNamespace(pack, quietly = TRUE)) {install.packages(pack)}})
lapply(packs, require, character.only = TRUE)



# Load Data ---------------------------------------------------------------
biospecimen <- fromJSON("data/biospecimen.json")
clinical <- fromJSON("data/clinical.json")




