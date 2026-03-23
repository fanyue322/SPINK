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
### 1. Load data

```
library(SPINK)

# Load example data (included in package)
data(object)
object
```

### 2. Preprocess
```
obj <- spink_preprocess(
  object = object,
  group.by = "domain",
  domain = "R3",           
  refGenome = "hg38",
  distance = 5e+05,          
  num.core = 2
)
```

### 3. Run analysis
```
obj <- spink_analysis(
  object = obj,
  link.gene.thr = 0.1,       
  n.permute = 1000           
)
```

### 4. Get results
```
results <- GetLinkResult(obj, domain = "R3")
head(results)
```
