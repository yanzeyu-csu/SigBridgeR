# **SigBridgeR** <a href="https://wanglabcsu.github.io/sigbridger/"><img src="man/figures/logo_white.png" alt="sigbridger website" align="right" height="139"/></a>

[![CRAN Status](https://www.r-pkg.org/badges/version/SigBridgeR)](https://cran.r-project.org/package=SigBridgeR) [![Repo Status](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)[![](https://img.shields.io/badge/devel%20version-0.0.0.9000-blue.svg)](https://github.com/WangLabCSU/SigBridgeR)

------------------------------------------------------------------------

## 🌐 Overview

SigBridgeR integrates multiple algorithms, using single-cell RNA sequencing data, bulk expression data, and sample-related phenotypic data, to identify the cells most closely associated with the phenotypic data, performing as a bridge to existing tools.

## 🔧 Installation

You can install the **SigBridgeR** using the following options:

### Stable release from CRAN

```{r install_from_cran}
install.packages("SigBridgeR")
```

### Development version from GitHub

```{r install_from_github}
if(!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("WangLabCSU/SigBridgeR")
```

## 📓 Documentation

Use `?SigBridgeR` to access the help documents in R.

Get Started:

-   [Quick Started Guide](vignettes/Quick_Start.md)
-   [Full Tutorial](vignettes/Full_Tutorial.md) for more details

If you encounter problems, please see:

-   [Troubleshooting Guide](vignettes/Troubleshooting.md)

Other information:

-   What is *Single Cell Sequencing*?
    -   [Veiw in Wiki](https://en.wikipedia.org/wiki/Single-cell_sequencing)
-   What is *RNA-seq*?
    -   [View in Wiki](https://en.wikipedia.org/wiki/RNA-Seq)

SigbridgeR integrates algorithms from the following repositories:

-   [Github-sunduanchen/Scissor](https://github.com/sunduanchen/Scissor)
-   [Github-Qinran-Zhang/scAB](https://github.com/Qinran-Zhang/scAB/)
-   [Github-WangX-Lab/ScPP](https://github.com/WangX-Lab/ScPP)
-   [Github-aiminXie/scPAS](https://github.com/aiminXie/scPAS)

## 📮 Contact

For support or questions:

Maintainer: Exceret 