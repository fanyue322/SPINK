# SPINK

**SP**atial l**INK**age of RNA and ATAC for spatial multi-ome data

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R](https://img.shields.io/badge/R-%3E%3D%204.0.0-blue)](https://cran.r-project.org/)

SPINK is an R package for inferring **enhancer-gene regulatory interactions** from spatially resolved multi-ome data. Built on a **hierarchical linear spatial model (HLSM)** with thin-plate splines, SPINK explicitly models spatial dependence to deliver well-calibrated p-values and superior statistical power compared to existing methods.

---

## Installation

```
# Install from GitHub
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("fanyue322/SPINK")
```

## Quick Start
Analyze a spatial domain in 4 steps:


```r
library(SPINK)

# 1. Load example data (included in package)
data(object)
object

# 2. Preprocess: spatial decorrelation + cis peak-gene mapping
obj <- spink_preprocess(
  object = object,
  group.by = "domain",    # metadata column with domain labels
  domain = "R3",          # target domain (CA region)  
  refGenome = "hg38",     # or "hg19", "mm10"
  distance = 5e+05,       # search window for cis-regulatory elements (bp)      
  num.core = 2            # parallel cores for TPS fitting
)

# 3. Infer regulatory links: constrained regression + permutation testing
obj <- spink_analysis(
  object = obj,
  link.gene.thr = 0.1,        # gene-level significance threshold 
  n.permute = 1000            # permutations for empirical p-values
)

# 4. Extract results: gene-level and peak-level associations
results <- GetLinkResult(obj, domain = "R3")
head(results)
```

## Usage
SPINK analyzes spatially co-profiled transcriptomic and epigenomic data from spatial multi-ome technologies (e.g., MISAR-seq, Spatial-RNA-ATAC, Slide-tags).

**Input Requirements:**
- A **Seurat object** containing:
  - RNA assay: log-normalized gene expression matrix (genes × spots)
  - ATAC assay: TF-IDF normalized chromatin accessibility matrix (peaks × spots)
  - Spatial coordinates in the `images` slot 
- **Identical spot barcodes** across RNA and ATAC assays
