
<!-- README.md is generated from README.Rmd. Please edit that file -->

# smokingMouse

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/LieberInstitute/smokingMouse)](https://github.com/LieberInstitute/smokingMouse/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/LieberInstitute/smokingMouse)](https://github.com/LieberInstitute/smokingMouse/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Bioc release
status](http://www.bioconductor.org/shields/build/release/data-experiment/smokingMouse.svg)](https://bioconductor.org/checkResults/release/data-experiment-LATEST/smokingMouse)
[![Bioc devel
status](http://www.bioconductor.org/shields/build/devel/data-experiment/smokingMouse.svg)](https://bioconductor.org/checkResults/devel/data-experiment-LATEST/smokingMouse)
[![Bioc downloads
rank](https://bioconductor.org/shields/downloads/release/smokingMouse.svg)](http://bioconductor.org/packages/stats/bioc/smokingMouse/)
[![Bioc
support](https://bioconductor.org/shields/posts/smokingMouse.svg)](https://support.bioconductor.org/tag/smokingMouse)
[![Bioc last
commit](https://bioconductor.org/shields/lastcommit/devel/data-experiment/smokingMouse.svg)](http://bioconductor.org/checkResults/devel/data-experiment-LATEST/smokingMouse/)
[![Bioc
dependencies](https://bioconductor.org/shields/dependencies/release/smokingMouse.svg)](https://bioconductor.org/packages/release/bioc/html/smokingMouse.html#since)
[![R-CMD-check-bioc](https://github.com/LieberInstitute/smokingMouse/actions/workflows/R-CMD-check-bioc.yaml/badge.svg)](https://github.com/LieberInstitute/smokingMouse/actions/workflows/R-CMD-check-bioc.yaml)

<!-- badges: end -->

Welcome to the `smokingMouse` project! Here you are able to access the
mouse and human datasets used for the analysis of the smoking mouse LIBD
project.

## Overview

This project consisted of a differential expression analysis (DEA)
involving 4 features: genes, exons, transcripts and junctions. The main
goal of this study was to explore the effects of smoking and nicotine
exposures on the developing brain of mice pups. As secondary objectives,
this work evaluated the affected genes by each substance on adult brain
in order to compare pup and adult results, and the effects of smoking on
adult blood and brain to search for overlapping biomarkers in both
tissues. Finally, differentially expressed (DE) genes in mice were
compared against results of a previous study on human to see how
consistent the results were.

The next table summarizes the analyses done at each level.

<figure>
<img src="man/figures/Table_of_Analyses.png" align="center" width="800px" />
<figcaption style="color: gray; line-height: 0.94; text-align: justify">

<font size="-1.5"><b>Table 1</b>: <b>1. Data preparation</b>: in this
first step, counts of genes, exons and junctions were normalized to CPM
and scaled; transcript expression values were only scaled since they
were already in TMP. Then, low-expression features were removed using
the indicated methods and samples were separated by tissue and age in
order to create subsets of the data for downstream analyses. <b>2.
Exploratory Data Analysis</b>: QC metrics of the samples were examined
and used to filter them; sample level effects were explored through
dimensionality reduction methods and rare samples in PCA plots were
manually removed from the datasets; gene level effects were evaluated
with analyses of explanatory variables and variance partition. <b>3.
Differential Expression Analysis</b>: with the relevant variables
identified in the previous steps, the DEA was performed at the gene
level for nicotine and smoking, adult and pup, and blood and brain
samples, and for 3 models: the naive one modeled <span
style="font-family: monospace">\~Group + batch effects</span>, the
adjusted model modeled <span style="font-family: monospace">\~Group +
Pregnancy + batch effects</span> for adults and <span
style="font-family: monospace">\~Group + Sex + batch effects</span> for
pups, and the interaction model <span
style="font-family: monospace">\~Group\*Pregnancy + batch effects</span>
for adults and <span style="font-family: monospace">\~Group\*Sex + batch
effects</span> for pups; DEA on the rest of the levels was performed for
pups only and using the adjusted model. After that, signals of the
features in nicotine and smoking were compared, as well as the signals
of exons and txs vs the effects of their genes, and genes’ signals were
additionally compared in the different tissues, ages, models and species
(vs human data of a previous study). All resultant DEG and DE features
(and their genes) were quantified and compared based on their experiment
(nic/smo) and direction of regulation (up/down); DEG were further
compared against genes of DE exons and txs; mouse genes were also
compared with human genes affected by cigarette smoke or associated with
TUD. <b>4. Gene Ontology and KEGG</b>: taking the DEG and the genes of
DE txs and exons, GO & KEGG analyses were done and the expression levels
of genes that participate in brain development related processes were
explored. <b>5. DE feature visualization</b>: DEG counts were
represented in heatmaps in order to distinguish the groups of up and
down-regulated genes. <b>6. Junction annotation</b>: for novel DE jxns
of unknown gene, their nearest, preceding and following genes were
determined.</font>

<font size="0.8">Abbreviations: Jxn: junction; Tx: transcript; CPM:
counts per million; TPM: transcripts per million; TMM: Trimmed Mean of
M-Values; TMMwsp: TMM with singleton pairing; EDA: exploratory data
analysis; QC: quality control; ribo: ribosomal; mt: mitochondrial; PCA:
Principal Component Analysis; PC: principal component; MDS:
Multidimensional Scaling; DEA: differential expression analysis; DE:
differential expression/differentially expressed; DEG: differentially
expressed genes; GO: Gene Ontology; KEGG: Kyoto Encyclopedia of Genes
and Genomes; TUD: tobacco use disorder.</font>

</figcaption>
</figure>

## Study design

<figure>
<img src="man/figures/Study_design.png" align="center" width="800px" />
<figcaption style="color: gray; line-height: 0.88; text-align: justify">

<font size="-1.5"><b>Figure 1</b>: Experimental design of the study.
<b>A)</b> 36 pregnant mice and 35 not pregnant female adults were either
administered nicotine (n=12), exposed to cigarette smoke (n=24) or
controls (n=35; 11 nicotine controls and 24 smoking controls). A total
of 137 pups were born to pregnant mice: 19 were born to mice that were
administered nicotine, 46 to mice exposed to smoking and the remaining
72 to control mice (23 to nicotine controls and 49 to smoking controls).
Samples from frontal cortices of P0 pups and adults were obtained, as
well as blood samples from smoking-exposed and smoking control adults.
<b>B)</b> RNA was extracted from that samples and RNA-seq experiments
were performed, obtaining expression counts for genes, exons,
transcripts and exon-exon junctions.</font>

</figcaption>
</figure>

## smoking Mouse datasets

The mouse datasets contain the following data in a single object for
each feature (genes, exons, transcripts and exon-exon junctions):

- **Raw data**: original read counts of the features, also including the
  original information of features and samples.
- **Processed data**: normalized and scaled counts of the same features.
  In addition to the feature and sample information, the datasets
  contain information of which ones were used in downstream analyses
  (the ones that passed filtering steps) and which features were DE in
  the different experiments.

Moreover, you can find human data generated by Semick, S.A. et
al. (2018) in Mol Psychiatry, DOI:
<https://doi.org/10.1038/s41380-018-0223-1>, that contain the results of
a DEA in adult and prenatal human brain samples exposed to cigarette
smoke.

## Data specifics

- *‘rse_gene_mouse_RNAseq_nic-smo.Rdata’*: (rse_gene object) the gene
  RSE object contains raw and normalized expression data of 55401 mouse
  genes across 208 samples from brains and blood of healthy and
  nicotine/smoking-exposed pup and adult mice.
- *‘rse_tx_mouse_RNAseq_nic-smo.Rdata’*: (rse_tx object) the tx RSE
  object contains raw and normalized expression data of 142604 mouse
  transcripts across 208 samples from brains and blood of healthy and
  nicotine/smoking-exposed pup and adult mice.
- *‘rse_exon_mouse_RNAseq_nic-smo.Rdata’*: (rse_exon object) the exon
  RSE object contains raw and normalized expression data of 447670 mouse
  exons across 208 samples from brains and blood of healthy and
  nicotine/smoking-exposed pup and adult mice.
- *‘rse_jx_mouse_RNAseq_nic-smo.Rdata’*: (rse_jx object) the jx RSE
  object contains raw and normalized expression data of 1436068 mouse
  exon-exon junctions across 208 samples from brains and blood of
  healthy and nicotine/smoking-exposed pup and adult mice.

All the above datasets contain sample and feature information and
additional data of the results obtained in the filtering steps and the
DEA.

- *‘de_genes_prenatal_human_brain_smoking.Rdata’*: (object with the same
  name) data frame with DE (ctrls vs smoking-exposed samples) data of
  18067 genes in human prenatal brain samples exposed to cigarette
  smoke.
- *‘de_genes_adult_human_brain_smoking.Rdata’*: (object with the same
  name) data frame with DE (ctrls vs smoking-exposed samples) data of
  18067 genes in human adult brain samples exposed to cigarette smoke.

## R/Bioconductor package

The `smokingMouse` package contains functions for:

- Accessing the expression data from the LIBD smoking Mouse project
  ([code on
  GitHub](https://github.com/LieberInstitute/smokingMouse_Indirects)).
  The datasets are retrieved from
  [Bioconductor](http://bioconductor.org/)’s `ExperimentHub`.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `smokingMouse` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("smokingMouse")
```

And the development version from
[GitHub](https://github.com/LieberInstitute/smokingMouse) with:

``` r
BiocManager::install("LieberInstitute/smokingMouse")
```

## Example of how to access the data

Through the `smokingMouse` package you can access the mouse datasets of
the project that include the raw and processed data. Below there’s code
you can use to access the gene data but can do the same for any of the
datasets. For more details, check the documentation for
`RangedSummarizedExperiment` objects. You can also find code to access
human data.

``` r
## Connect to ExperimentHub
ehub <- ExperimentHub::ExperimentHub()
#> snapshotDate(): 2022-10-31
```

``` r

###########      The following is just provisional code for the example.      ###########
########### Once the package is approved and ready to access, the example will be correctly updated. ###########
 
## Load the package
library("smokingMouse")

########################
#      Mouse data 
########################
# ## Download the mouse gene data
# rse_gene <- function_to_access
# ## This is a RangedSummarizedExperiment object
# rse_gene

# ## Note the memory size
# lobstr::obj_size(rse_gene)
# 
# ## Check sample info 
# head(colData(rse_gene), 3)
# ## Check gene info
# head(rowData(rse_gene), 3)
# ## Access the original counts
# original_counts <- assays(rse_gene)$counts
# ## Access the log normalized counts
# logcounts <- assays(rse_gene)$logcounts



########################
#      Human data 
########################
# ## Download the human gene data
# de_genes_prenatal_human_brain_smoking <- function_to_access
# ## This is a data frame
# de_genes_prenatal_human_brain_smoking

# ## Note the memory size
# lobstr::obj_size(de_genes_prenatal_human_brain_smoking)

## Access data of human genes as normally do with data frames
```

## Citation

Below is the citation output from using `citation('smokingMouse')` in R.
Please run this yourself to check for any updates on how to cite
**smokingMouse**.

``` r
print(citation('smokingMouse'), bibtex = TRUE)
#> 
#> To cite package 'smokingMouse' in publications use:
#> 
#>   Gonzalez-Padilla D, Collado-Torres L (2023). _Provides access to
#>   smokingMouse project data _. doi:10.18129/B9.bioc.smokingMouse
#>   <https://doi.org/10.18129/B9.bioc.smokingMouse>,
#>   https://github.com/LieberInstitute/smokingMouse/smokingMouse - R
#>   package version 0.99.0,
#>   <http://www.bioconductor.org/packages/smokingMouse>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {Provides access to smokingMouse project data },
#>     author = {Daianna Gonzalez-Padilla and Leonardo Collado-Torres},
#>     year = {2023},
#>     url = {http://www.bioconductor.org/packages/smokingMouse},
#>     note = {https://github.com/LieberInstitute/smokingMouse/smokingMouse - R package version 0.99.0},
#>     doi = {10.18129/B9.bioc.smokingMouse},
#>   }
#> 
#>   Gonzalez-Padilla D, Collado-Torres L (2023). "Provides access to
#>   smokingMouse project data." _bioRxiv_. doi:10.1101/TODO
#>   <https://doi.org/10.1101/TODO>,
#>   <https://www.biorxiv.org/content/10.1101/TODO>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Provides access to smokingMouse project data},
#>     author = {Daianna Gonzalez-Padilla and Leonardo Collado-Torres},
#>     year = {2023},
#>     journal = {bioRxiv},
#>     doi = {10.1101/TODO},
#>     url = {https://www.biorxiv.org/content/10.1101/TODO},
#>   }
```

Please note that the `smokingMouse` and the [smoking
mouse](https://github.com/LieberInstitute/smokingMouse_Indirects)
project were only made possible thanks to many other R and
bioinformatics software authors, which are cited either in the vignettes
and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `smokingMouse` project is released with a
[Contributor Code of
Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.16/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation
  website](http://LieberInstitute.github.io/smokingMouse) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.16/biocthis)*.
