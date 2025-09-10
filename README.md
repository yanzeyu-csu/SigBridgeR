# **SigBridgeR** <a href="https://wanglabcsu.github.io/SigBridgeR/"><img src="man/figures/logo_white.png" alt="sigbridger website" align="right" height="139"/></a>

[![Repo Status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) [![License: GPL3](https://img.shields.io/badge/license-GPL3-blue.svg)](https://cran.r-project.org/web/licenses/GPL3) [![](https://img.shields.io/badge/devel%20version-2.2.3-blue.svg)](https://github.com/WangLabCSU/SigBridgeR) [![R CMD check](https://github.com/WangLabCSU/SigBridgeR/workflows/R-CMD-check/badge.svg)](https://github.com/WangLabCSU/SigBridgeR/actions) [![registry status badge](https://wanglabcsu.r-universe.dev/badges/:registry)](https://wanglabcsu.r-universe.dev/)

------------------------------------------------------------------------

## 🌐 Overview

SigBridgeR integrates multiple algorithms, using single-cell RNA sequencing data, bulk expression data, and sample-related phenotypic data, to identify the cells most closely associated with the phenotypic data, performing as a bridge to existing tools.

## 🔧 Installation

Usually we recommend installing the latest release from GitHub because of the latest features and bug fixes.

1.  Install the development version from GitHub:

``` r
if(!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("WangLabCSU/SigBridgeR")
```

2.  Install from r-universe:

``` r
install.packages("SigBridgeR", repos = "https://wanglabcsu.r-universe.dev")
```

## 📓 Documentation

Get Started:

-   [Quick Started Guide](vignettes/Quick_Start.md)
-   [Full Tutorial](vignettes/Full_Tutorial.md) for more details
-   Use `?SigBridgeR::function_name` to access the help documents in R.

If you encounter problems, please see:

-   [Troubleshooting Guide](vignettes/Troubleshooting.md)
-   Please use [Github issues](https://github.com/WangLabCSU/SigBridgeR/issues) if you want to file bug reports or feature requests, let us know if you have ideas to make this project better!

Other information:

-   What is *Single Cell Sequencing*?
    -   [Veiw in Wiki](https://en.wikipedia.org/wiki/Single-cell_sequencing)
-   What is *RNA-seq*?
    -   [View in Wiki](https://en.wikipedia.org/wiki/RNA-Seq)