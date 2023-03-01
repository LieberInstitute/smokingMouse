
library("here")
library("sessioninfo")


##################
##   Raw data
##################

## Load

load(here("~/Desktop/smokingMouse_Indirects/raw-data/rse_gene_smoking_mouse_n208.Rdata"), verbose = TRUE)
load(here("~/Desktop/smokingMouse_Indirects/raw-data/rse_tx_smoking_mouse_n208.Rdata"), verbose = TRUE)
load(here("~/Desktop/smokingMouse_Indirects/raw-data/rse_jx_smoking_mouse_n208.Rdata"), verbose = TRUE)
load(here("~/Desktop/smokingMouse_Indirects/raw-data/rse_exon_smoking_mouse_n208.Rdata"), verbose = TRUE)
load(here("~/Desktop/smokingMouse_Indirects/raw-data/Genes_DE_sva.rda"), verbose = TRUE)
Maternal_Smoking_pheno <- read.table("~/Desktop/smokingMouse_Indirects/raw-data/Maternal_Smoking_pheno.txt")


## Add metadata to human data from Semick SA et al. (2018)
metadata(fetalGene) <- list(
  "Downloaded_from"="https://github.com/LieberInstitute/Smoking_DLPFC_Devel",
  "Cite_this_paper"="https://www.nature.com/articles/s41380-018-0223-1"
)

metadata(adultGene) <- list(
  "Downloaded_from"="https://github.com/LieberInstitute/Smoking_DLPFC_Devel",
  "Cite_this_paper"="https://www.nature.com/articles/s41380-018-0223-1"
)


## Save
save(rse_gene, file="rse_gene_smoking_mouse_n208.Rdata")
save(rse_tx, file="rse_tx_smoking_mouse_n208.Rdata")
save(rse_jx, file="rse_jx_smoking_mouse_n208.Rdata")
save(rse_exon, file="rse_exon_smoking_mouse_n208.Rdata")
save(fetalGene, file="fetalGene.Rdata")
save(adultGene, file="adultGene.Rdata")
save(Maternal_Smoking_pheno, file="Maternal_Smoking_pheno.txt")





## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# setting  value
# version  R version 4.2.2 (2022-10-31)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Mexico_City
# date     2023-02-28
# rstudio  2022.12.0+353 Elsbeth Geranium (desktop)
# pandoc   NA
