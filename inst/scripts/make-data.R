
## This was run locally but can be ran on JHPCE at. ....????

library(here)
library(sessioninfo)
library(SummarizedExperiment)
library(GenomicRanges)

## Load rse objects and human data

load(here("~/Desktop/smokingMouse_Indirects/smokingMouse_pkg/inst/extdata/rse_gene_complete.Rdata"), verbose = TRUE)
# rse_gene
load(here("~/Desktop/smokingMouse_Indirects/smokingMouse_pkg/inst/extdata/rse_exon_complete.Rdata"), verbose = TRUE)
# rse_exon
load(here("~/Desktop/smokingMouse_Indirects/smokingMouse_pkg/inst/extdata/rse_tx_complete.Rdata"), verbose = TRUE)
# rse_tx
load(here("~/Desktop/smokingMouse_Indirects/smokingMouse_pkg/inst/extdata/rse_jx_complete.Rdata"), verbose = TRUE)
# rse_jx
load(here("~/Desktop/smokingMouse_Indirects/raw-data/Genes_DE_sva.rda"), verbose = TRUE)
# fetalGene
# adultGene

## Make GRanges from human data
de_genes_prenatal_human_brain_smoking<- makeGRangesFromDataFrame(fetalGene, keep.extra.columns = TRUE)
de_genes_adult_human_brain_smoking <- makeGRangesFromDataFrame(adultGene, keep.extra.columns = TRUE)


## Add metadata to rse objects
metadata(rse_gene) <- list(
  "Obtained_from"="https://github.com/LieberInstitute/smokingMouse_Indirects/blob/main/smokingMouse_pkg/inst/scripts/make-data_smokingMouse.R"
)

metadata(rse_exon) <- list(
  "Obtained_from"="https://github.com/LieberInstitute/smokingMouse_Indirects/blob/main/smokingMouse_pkg/inst/scripts/make-data_smokingMouse.R"
)

metadata(rse_tx) <- list(
  "Obtained_from"="https://github.com/LieberInstitute/smokingMouse_Indirects/blob/main/smokingMouse_pkg/inst/scripts/make-data_smokingMouse.R"
)

metadata(rse_jx) <- list(
  "Obtained_from"="https://github.com/LieberInstitute/smokingMouse_Indirects/blob/main/smokingMouse_pkg/inst/scripts/make-data_smokingMouse.R"
)

## Add metadata to human data from Semick SA et al. (2018)
metadata(de_genes_prenatal_human_brain_smoking) <- list(
  "Downloaded_from"="https://github.com/LieberInstitute/Smoking_DLPFC_Devel",
  "Cite_this_paper"="https://www.nature.com/articles/s41380-018-0223-1"
)

metadata(de_genes_adult_human_brain_smoking) <- list(
  "Downloaded_from"="https://github.com/LieberInstitute/Smoking_DLPFC_Devel",
  "Cite_this_paper"="https://www.nature.com/articles/s41380-018-0223-1"
)



## Save rse objects under more informative names
save(rse_gene, file="~/Desktop/smokingMouse/inst/extdata/rse_gene_mouse_RNAseq_nic-smo.Rdata", version = 2)
save(rse_tx, file="~/Desktop/smokingMouse/inst/extdata/rse_tx_mouse_RNAseq_nic-smo.Rdata", version = 2)
save(rse_jx, file="~/Desktop/smokingMouse/inst/extdata/rse_jx_mouse_RNAseq_nic-smo.Rdata", version = 2)
save(rse_exon, file="~/Desktop/smokingMouse/inst/extdata/rse_exon_mouse_RNAseq_nic-smo.Rdata", version = 2)
save(de_genes_prenatal_human_brain_smoking, file="~/Desktop/smokingMouse/inst/extdata/de_genes_prenatal_human_brain_smoking.Rdata", version = 2)
save(de_genes_adult_human_brain_smoking, file="~/Desktop/smokingMouse/inst/extdata/de_genes_adult_human_brain_smoking.Rdata", version = 2)







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
