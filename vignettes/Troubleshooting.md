# Overview

This vignette provides solutions to common problems you might encounter
when using this package. If you cannot find your problem here, please
file an issue on our GitHub repository.

------------------------------------------------------------------------

> ### Error in normalize.quantiles(dataset0) : ERROR; return code from pthread\_create() is 22

To solve this problem, you need to install `preprocessCore` without
threading support. Try:

    git clone https://github.com/bmbolstad/preprocessCore.git
    cd preprocessCore
    R CMD INSTALL --configure-args="--disable-threading"  .

or

    BiocManager::install(
      "preprocessCore",
      configure.args = "--disable-threading",
      force = TRUE
    )

See
[bioconductor\_docker/issues/22](https://github.com/Bioconductor/bioconductor_docker/issues/22),
[Scissor/issues/15](https://github.com/sunduanchen/Scissor/issues/15)
for more details.

------------------------------------------------------------------------

> ### Error at alpha=0.05:subscript out of bounds
>
> ### Error in `Scissor.v5.optimized()`:
>
> ### ! object ‘fit0’ not found

This may be due to two reasons: first, a mismatch in the dimensions of
the bulk expression data and the phenotype data; second, incorrect
column names in the survival phenotype data leading to a failure to
match.

To check dimension mismatch：

    ncol(your_bulk_data) == nrow(your_phenotype_data) # should be TRUE

    all(unique(colnames(your_bulk_data)) == unique(rownames(your_bulk_data))) # should be TRUE

    all(order(colnames(your_bulk_data)) == order(rownames(your_phenotype_data))) # should be TRUE

Survival phenotype column names should be formatted as `time` and
`status`, ensuring correct capitalization and spelling:

    head(survival_phenotype) # case-sensitive
    #           time status
    # GSM70130 34.80      0
    # GSM70131 35.67      0
    # GSM70136 43.37      0
    # GSM70138 60.77      0
    # GSM70140 33.80      1
    # GSM70144 58.53      0

------------------------------------------------------------------------

> ### Error in `function_name()`:
>
> ### ! lazy-load database ‘/home/user/R/x86\_64-pc-linux-gnu-library/4.4/SigBridgeR/R/SigBridgeR. rdb’ is corrupt

An error occurred during installation, causing the package to be
corrupted. Try reinstalling the package:

    # I reccomend restarting R/RStudio before reinstalling
    remove.packages("SigBridgeR")

    detach("package:SigBridgeR", unload = TRUE)

    if (!requireNamespace("remotes")) {
      install.packages("remotes")
    }
    remotes::install_github("WangLabCSU/SigBridgeR")

------------------------------------------------------------------------

> ### Error:
>
> ### ! Invalid syntax: ‘c(scissor\_umap, scpas\_umap)’

As far as I know, this is due to the R environment being contaminated,
which prevents the use of `%<-%`. You can try restarting the R
environment or clean up the environment.

    # restart R/RStudio
    .rs.api.restartSession()

    rm(list = ls(all.names = TRUE))

------------------------------------------------------------------------

> ### Error:
>
> ### ! Detected n gene(s) with zero variance:
>
> ### ℹ “gene name(s)”

This is due to the presence of genes with zero variance in the bulk
expression data when you are using `scPP` and `binary` phenotype. This
indicates that the expression levels of one (or several) genes are
nearly identical across different samples. You should check your data.
