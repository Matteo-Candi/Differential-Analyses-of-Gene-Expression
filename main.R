rm(list = ls())

# Set Reproducibility -----------------------------------------------------
set.seed(1234)



# Load Packages -----------------------------------------------------------
packs <- c('ggplot2', 'jsonlite', 'readr')
lapply(packs, function(pack) {
  if (!requireNamespace(pack, quietly = TRUE)) {install.packages(pack)}})
lapply(packs, require, character.only = TRUE)



# Load and Cleaning Data --------------------------------------------------

# JSON
biospecimen <- fromJSON("data/json/biospecimen.json")
clinical <- fromJSON("data/json/clinical.json")

filtered_biospecimen <- fromJSON("data/json/filtered_biospecimen.json")
filtered_clinical <- fromJSON("data/json/filtered_clinical.json")


# TSV
bio_aliquot <- read_tsv("data/tsv/biospecimen/aliquot.tsv")
bio_analyte <- read_tsv("data/tsv/biospecimen/analyte.tsv")
bio_portion <- read_tsv("data/tsv/biospecimen/portion.tsv")
bio_sample <- read_tsv("data/tsv/biospecimen/sample.tsv")
bio_slide <- read_tsv("data/tsv/biospecimen/slide.tsv")

clin_clinical <- read_tsv("data/tsv/clinical/clinical.tsv")
clin_exposure <- read_tsv("data/tsv/clinical/exposure.tsv")
clin_family_history <- read_tsv("data/tsv/clinical/family_history.tsv")
clin_follow_up <- read_tsv("data/tsv/clinical/follow_up.tsv")
clin_pathology_detail <- read_tsv("data/tsv/clinical/pathology_detail.tsv")


# Replace '-- with NA.
replace_with_na <- function(data) {
  data[data == "'--"] <- NA
  return(data)
}

bio_aliquot <- replace_with_na(bio_aliquot)
bio_analyte <- replace_with_na(bio_analyte)
bio_portion <- replace_with_na(bio_portion)
bio_sample <- replace_with_na(bio_sample)
bio_slide <- replace_with_na(bio_slide)

clin_clinical <- replace_with_na(clin_clinical)
clin_exposure <- replace_with_na(clin_exposure)
clin_family_history <- replace_with_na(clin_family_history)
clin_follow_up <- replace_with_na(clin_follow_up)
clin_pathology_detail <- replace_with_na(clin_pathology_detail)
