
library(here)
library(sessioninfo)


outdir <- "smokingMouse_files"
pkgname <- "smokingMouse"


meta <- data.frame(
  Title = c(
    "rse_gene_mouse_RNAseq_nic-smo",
    "rse_tx_mouse_RNAseq_nic-smo",
    "rse_jx_mouse_RNAseq_nic-smo",
    "rse_exon_mouse_RNAseq_nic-smo",
    "de_genes_prenatal_human_brain_smoking",
    "de_genes_adult_human_brain_smoking"
  ),
  Description = c(
    "RangedSummarizedExperiment with bulk RNA expression data of mouse genes in brain/blood of nicotine/smoking-exposed pup/adult mice; generated at the Lieber Institute for Brain Development (LIBD) and available through the smokingMouse Bioconductor package.",
    "RangedSummarizedExperiment with bulk RNA expression data of mouse transcripts (tx) in brain/blood of nicotine/smoking-exposed pup/adult mice; generated at the Lieber Institute for Brain Development (LIBD) and available through the smokingMouse Bioconductor package.",
    "RangedSummarizedExperiment with bulk RNA expression data of mouse exon-exon junctions (jx) in brain/blood of nicotine/smoking-exposed pup/adult mice; generated at the Lieber Institute for Brain Development (LIBD) and available through the smokingMouse Bioconductor package.",
    "RangedSummarizedExperiment with bulk RNA expression data of mouse exons in brain/blood of nicotine/smoking-exposed pup/adult mice; generated at the Lieber Institute for Brain Development (LIBD) and available through the smokingMouse Bioconductor package.",
    "GRanges with the information of differential expression of human genes in the cigarette smoke-exposed human prenatal brain, generated in Semick et al. 2020 (DOI: https://doi.org/10.1038/s41380-018-0223-1) and obtained from https://github.com/LieberInstitute/Smoking_DLPFC_Devel.",
    "GRanges with the information of differential expression of human genes in the cigarette smoke-exposed human adult brain, generated in Semick et al. 2020 (DOI: https://doi.org/10.1038/s41380-018-0223-1) and obtained from https://github.com/LieberInstitute/Smoking_DLPFC_Devel."
  ),
  BiocVersion = "3.18",
  Genome = rep(c("GRCm38", "GRCh37"), c(4,2)),
  SourceType = "GTF",
  SourceUrl = "https://bioconductor.org/packages/smokingMouse",
  SourceVersion = "Feb 28 2023",
  Species = rep(c("Mus musculus", "Homo sapiens"), c(4,2)),
  TaxonomyId = rep(c("10090", "9606"),c(4,2)),
  Coordinate_1_based = TRUE,
  DataProvider = "Lieber Institute for Brain Development (LIBD)",
  Maintainer = "Daianna Gonzalez-Padilla <glezdaianna@gmail.com>",
  RDataClass = rep(c("RangedSummarizedExperiment", "GenomicRanges"), c(4,2)),
  DispatchClass = "Rda",
  RDataPath = file.path(
    pkgname,
    outdir,
    c(
      "rse_gene_mouse_RNAseq_nic-smo.Rdata",
      "rse_tx_mouse_RNAseq_nic-smo.Rdata",
      "rse_jx_mouse_RNAseq_nic-smo.Rdata",
      "rse_exon_mouse_RNAseq_nic-smo.Rdata",
      "de_genes_prenatal_human_brain_smoking.Rdata",
      "de_genes_adult_human_brain_smoking.Rdata"
    )
  ),
  Tags = "smokingMouse:LIBD:mouse:smoking:nicotine:brain:blood:development",
  row.names = NULL,
  stringsAsFactors = FALSE
)

write.csv(
  meta,
  file = here::here("inst", "extdata", "metadata.csv"),
  row.names = FALSE
)

## Check
if (FALSE) {
  AnnotationHubData::makeAnnotationHubMetadata(here::here(), fileName = "metadata.csv")
}





## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.0 (2023-04-21)
# os       macOS Monterey 12.5.1
# system   aarch64, darwin20
# ui       RStudio
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       America/Mexico_City
# date     2024-02-08
# rstudio  2023.06.1+524 Mountain Hydrangea (desktop)
# pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# !  package                * version   date (UTC) lib source
# AnnotationDbi            1.63.2    2023-07-03 [1] Bioconductor
# AnnotationForge          1.43.3    2023-06-30 [1] Bioconductor
# AnnotationHub            3.9.1     2023-06-14 [1] Bioconductor
# AnnotationHubData        1.31.0    2023-07-07 [1] Bioconductor
# Biobase                * 2.61.0    2023-06-02 [1] Bioconductor
# BiocCheck                1.37.7    2023-07-18 [1] Bioconductor
# BiocFileCache            2.9.1     2023-07-14 [1] Bioconductor
# BiocGenerics           * 0.48.1    2023-11-02 [1] Bioconductor
# BiocIO                   1.11.0    2023-06-02 [1] Bioconductor
# BiocManager            * 1.30.21.1 2023-07-18 [1] CRAN (R 4.3.0)
# BiocParallel             1.35.3    2023-07-07 [1] Bioconductor
# BiocStyle                2.29.1    2023-07-19 [1] Bioconductor
# biocthis                 1.11.4    2023-07-05 [1] Bioconductor
# BiocVersion              3.18.0    2023-05-11 [1] Bioconductor
# biocViews                1.69.1    2023-06-02 [1] Bioconductor
# biomaRt                  2.57.1    2023-06-14 [1] Bioconductor
# Biostrings               2.69.2    2023-07-05 [1] Bioconductor
# bit                      4.0.5     2022-11-15 [1] CRAN (R 4.3.0)
# bit64                    4.0.5     2020-08-30 [1] CRAN (R 4.3.0)
# bitops                   1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
# blob                     1.2.4     2023-03-17 [1] CRAN (R 4.3.0)
# brio                     1.1.3     2021-11-30 [1] CRAN (R 4.3.0)
# cachem                   1.0.8     2023-05-01 [1] CRAN (R 4.3.0)
# callr                    3.7.3     2022-11-02 [1] CRAN (R 4.3.0)
# cli                      3.6.1     2023-03-23 [1] CRAN (R 4.3.0)
# codetools                0.2-19    2023-02-01 [1] CRAN (R 4.3.0)
# commonmark               1.9.0     2023-03-17 [1] CRAN (R 4.3.0)
# crayon                   1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# curl                     5.0.1     2023-06-07 [1] CRAN (R 4.3.0)
# DBI                      1.1.3     2022-06-18 [1] CRAN (R 4.3.0)
# dbplyr                   2.3.3     2023-07-07 [1] CRAN (R 4.3.0)
# DelayedArray             0.26.6    2023-07-02 [1] Bioconductor
# desc                     1.4.2     2022-09-08 [1] CRAN (R 4.3.0)
# devtools                 2.4.5     2022-10-11 [1] CRAN (R 4.3.0)
# digest                   0.6.33    2023-07-07 [1] CRAN (R 4.3.0)
# dplyr                    1.1.2     2023-04-20 [1] CRAN (R 4.3.0)
# ellipsis                 0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
# evaluate                 0.21      2023-05-05 [1] CRAN (R 4.3.0)
# ExperimentHub            2.9.1     2023-07-14 [1] Bioconductor
# fansi                    1.0.5     2023-10-08 [1] CRAN (R 4.3.1)
# fastmap                  1.1.1     2023-02-24 [1] CRAN (R 4.3.0)
# filelock                 1.0.2     2018-10-05 [1] CRAN (R 4.3.0)
# formatR                  1.14      2023-01-17 [1] CRAN (R 4.3.0)
# fs                       1.6.3     2023-07-20 [1] CRAN (R 4.3.0)
# futile.logger          * 1.4.3     2016-07-10 [1] CRAN (R 4.3.0)
# futile.options           1.0.1     2018-04-20 [1] CRAN (R 4.3.0)
# generics                 0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb           * 1.37.2    2023-06-21 [1] Bioconductor
# GenomeInfoDbData         1.2.10    2023-05-28 [1] Bioconductor
# GenomicAlignments        1.37.0    2023-07-07 [1] Bioconductor
# GenomicFeatures          1.53.1    2023-06-22 [1] Bioconductor
# GenomicRanges          * 1.54.1    2023-10-30 [1] Bioconductor
# glue                     1.6.2     2022-02-24 [1] CRAN (R 4.3.0)
# graph                    1.79.0    2023-06-02 [1] Bioconductor
# here                   * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
# hms                      1.1.3     2023-03-21 [1] CRAN (R 4.3.0)
# htmltools                0.5.5     2023-03-23 [1] CRAN (R 4.3.0)
# htmlwidgets              1.6.2     2023-03-17 [1] CRAN (R 4.3.0)
# httpuv                   1.6.11    2023-05-11 [1] CRAN (R 4.3.0)
# httr                     1.4.6     2023-05-08 [1] CRAN (R 4.3.0)
# interactiveDisplayBase   1.39.0    2023-06-02 [1] Bioconductor
# IRanges                * 2.36.0    2023-10-26 [1] Bioconductor
# jsonlite                 1.8.8     2023-12-04 [1] CRAN (R 4.3.1)
# KEGGREST                 1.41.0    2023-07-07 [1] Bioconductor
# knitr                    1.43      2023-05-25 [1] CRAN (R 4.3.0)
# lambda.r                 1.2.4     2019-09-18 [1] CRAN (R 4.3.0)
# later                    1.3.1     2023-05-02 [1] CRAN (R 4.3.0)
# lattice                  0.21-8    2023-04-05 [1] CRAN (R 4.3.0)
# lifecycle                1.0.3     2022-10-07 [1] CRAN (R 4.3.0)
# magrittr                 2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# Matrix                   1.6-4     2023-11-30 [1] CRAN (R 4.3.1)
# MatrixGenerics         * 1.13.0    2023-05-20 [1] Bioconductor
# matrixStats            * 1.0.0     2023-06-02 [1] CRAN (R 4.3.0)
# memoise                  2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
# mime                     0.12      2021-09-28 [1] CRAN (R 4.3.0)
# miniUI                   0.1.1.1   2018-05-18 [1] CRAN (R 4.3.0)
# OrganismDbi              1.43.0    2023-06-14 [1] Bioconductor
# pillar                   1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgbuild                 1.4.2     2023-06-26 [1] CRAN (R 4.3.0)
# pkgconfig                2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# pkgload                  1.3.2.1   2023-07-08 [1] CRAN (R 4.3.0)
# png                      0.1-8     2022-11-29 [1] CRAN (R 4.3.0)
# prettyunits              1.1.1     2020-01-24 [1] CRAN (R 4.3.0)
# processx                 3.8.2     2023-06-30 [1] CRAN (R 4.3.0)
# profvis                  0.3.8     2023-05-02 [1] CRAN (R 4.3.0)
# progress                 1.2.2     2019-05-16 [1] CRAN (R 4.3.0)
# promises                 1.2.0.1   2021-02-11 [1] CRAN (R 4.3.0)
# ps                       1.7.5     2023-04-18 [1] CRAN (R 4.3.0)
# purrr                    1.0.1     2023-01-10 [1] CRAN (R 4.3.0)
# R.cache                  0.16.0    2022-07-21 [1] CRAN (R 4.3.0)
# R.methodsS3              1.8.2     2022-06-13 [1] CRAN (R 4.3.0)
# R.oo                     1.25.0    2022-06-12 [1] CRAN (R 4.3.0)
# R.utils                  2.12.2    2022-11-11 [1] CRAN (R 4.3.0)
# R6                       2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# rappdirs                 0.3.3     2021-01-31 [1] CRAN (R 4.3.0)
# RBGL                     1.77.1    2023-06-02 [1] Bioconductor
# Rcpp                     1.0.11    2023-07-06 [1] CRAN (R 4.3.0)
# RCurl                    1.98-1.12 2023-03-27 [1] CRAN (R 4.3.0)
# remotes                  2.4.2.1   2023-07-18 [1] CRAN (R 4.3.0)
# restfulr                 0.0.15    2022-06-16 [1] CRAN (R 4.3.0)
# rjson                    0.2.21    2022-01-09 [1] CRAN (R 4.3.0)
# rlang                    1.1.1     2023-04-28 [1] CRAN (R 4.3.0)
# rmarkdown                2.23      2023-07-01 [1] CRAN (R 4.3.0)
# roxygen2                 7.2.3     2022-12-08 [1] CRAN (R 4.3.0)
# rprojroot                2.0.3     2022-04-02 [1] CRAN (R 4.3.0)
# Rsamtools                2.17.0    2023-07-07 [1] Bioconductor
# RSQLite                  2.3.1     2023-04-03 [1] CRAN (R 4.3.0)
# rstudioapi               0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
# rtracklayer              1.61.0    2023-07-07 [1] Bioconductor
# RUnit                    0.4.32    2018-05-18 [1] CRAN (R 4.3.0)
# S4Arrays                 1.1.4     2023-06-02 [1] Bioconductor
# S4Vectors              * 0.40.2    2023-11-25 [1] Bioconductor 3.18 (R 4.3.2)
# sessioninfo            * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
# shiny                    1.7.4.1   2023-07-06 [1] CRAN (R 4.3.0)
# VP smokingMouse           * 0.99.91   2023-06-28 [?] Github (LieberInstitute/smokingMouse@2e7640c) (on disk 0.99.5)
# stringdist               0.9.10    2022-11-07 [1] CRAN (R 4.3.0)
# stringi                  1.7.12    2023-01-11 [1] CRAN (R 4.3.0)
# stringr                  1.5.0     2022-12-02 [1] CRAN (R 4.3.0)
# styler                   1.10.1    2023-06-05 [1] CRAN (R 4.3.0)
# SummarizedExperiment   * 1.30.2    2023-06-06 [1] Bioconductor
# testthat                 3.2.1     2023-12-02 [1] CRAN (R 4.3.1)
# tibble                   3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyselect               1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# urlchecker               1.0.1     2021-11-30 [1] CRAN (R 4.3.0)
# usethis                  2.2.2     2023-07-06 [1] CRAN (R 4.3.0)
# utf8                     1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                    0.6.4     2023-10-12 [1] CRAN (R 4.3.1)
# withr                    2.5.2     2023-10-30 [1] CRAN (R 4.3.1)
# xfun                     0.39      2023-04-20 [1] CRAN (R 4.3.0)
# XML                      3.99-0.14 2023-03-19 [1] CRAN (R 4.3.0)
# xml2                     1.3.5     2023-07-06 [1] CRAN (R 4.3.0)
# xtable                   1.8-4     2019-04-21 [1] CRAN (R 4.3.0)
# XVector                  0.41.1    2023-06-02 [1] Bioconductor
# yaml                     2.3.8     2023-12-11 [1] CRAN (R 4.3.1)
# zlibbioc                 1.47.0    2023-05-20 [1] Bioconductor
# 
# [1] /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library
# 
# V ── Loaded and on-disk version mismatch.
# P ── Loaded and on-disk path mismatch.
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
