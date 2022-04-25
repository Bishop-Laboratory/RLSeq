# RLSeq <img src="https://rlbase-data.s3.amazonaws.com/misc/assets/whitebgRLSeq+Logo.png" align="right" alt="logo" width="240" style = "border: none; float: right;">

<!-- badges: start -->

[![BiocCheck](https://github.com/Bishop-Laboratory/RLSeq/workflows/BiocCheck/badge.svg)](https://github.com/Bishop-Laboratory/RLSeq/actions) [![Codecov test coverage](https://codecov.io/gh/Bishop-Laboratory/RLSeq/branch/main/graph/badge.svg)](https://codecov.io/gh/Bishop-Laboratory/RLSeq?branch=main) [![BioC status](http://www.bioconductor.org/shields/build/release/bioc/RLSeq.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/RLSeq)

<!-- badges: end -->

# Introduction

*RLSeq* (part of [*RLSuite*](https://gccri.bishop-lab.uthscsa.edu/rlsuite/)) is used for downstream analysis of R-loop datasets. It provides methods for data quality control and exploratory analysis within the context of the hundreds of publicly-available R-loop mapping data sets provided by [RLBase](https://github.com/Bishop-Laboratory/RLBase) and accessed via [RLHub](https://github.com/Bishop-Laboratory/RLHub). Finally, *RLSeq* provides a user-friendly HTML report that summarizes the analysis results.

**NOTE**: To run *RLSeq* in your browser, please see [*RLBase*](https://gccri.bishop-lab.uthscsa.edu/rlbase/). 

## Installation

### From Bioconductor

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RLSeq")
```

### From Github

1. Update to the `devel` version of bioconductor. 

```r
BiocManager::install(version = "devel")
```

2. Install **RLHub** and **RLSeq** with remotes.

``` r
remotes::install_github("Bishop-Laboratory/RLHub")
remotes::install_github("Bishop-Laboratory/RLSeq")
```

## Quick-start

This is an example workflow using a publicly-available R-loop mapping data set that
was reprocessed and standardized in [*RLBase*](https://gccri.bishop-lab.uthscsa.edu/rlbase/).

```r
# Peaks and coverage can be found in RLBase
rlbase <- "https://rlbase-data.s3.amazonaws.com"
pks <- file.path(rlbase, "peaks", "SRX1025890_hg38.broadPeak")
cvg <- file.path(rlbase, "coverage", "SRX1025890_hg38.bw")

# Initialize data in the RLRanges object. 
# Metadata is optional, but improves the interpretability of results
rlr <- RLRanges(
  peaks = pks,
  coverage = cvg,
  genome = "hg38",
  mode = "DRIP",
  label = "POS",
  sampleName = "TC32 DRIP-Seq"
)

# The RLSeq command performs all analyses
rlr <- RLSeq(rlr)

# Generate an html report
report(rlr)
```

The code above performs a typical analysis. It builds the `RLRanges` object, an extension of `GRanges` for
use with RLSeq. Then, it runs all core analyses using `RLSeq()`. Finally, it generates an
HTML report with `report()` (see the report 
[here](https://rlbase-data.s3.amazonaws.com/misc/rlseq_report_example.html)).

## Detail

For more information, see the package website [here](https://bishop-laboratory.github.io/RLSeq/).
