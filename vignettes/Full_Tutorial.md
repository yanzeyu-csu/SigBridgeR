# **Full Tutorial for SigBridgeR** {#full-tutorial-for-sigbridger}

## 0. Preface

### 0.1 Contents

- [**Full Tutorial for SigBridgeR** {#full-tutorial-for-sigbridger}](#full-tutorial-for-sigbridger-full-tutorial-for-sigbridger)
  - [0. Preface](#0-preface)
    - [0.1 Contents](#01-contents)
    - [0.1 Introduction to SigBridgeR](#01-introduction-to-sigbridger)
  - [1. Installation](#1-installation)
    - [1.1 Stable Release from GitHub](#11-stable-release-from-github)
    - [1.2 Check Dependencies](#12-check-dependencies)
  - [2. Loading and preprocessing data](#2-loading-and-preprocessing-data)
    - [2.1 Single-cell RNA-seq Data](#21-single-cell-rna-seq-data)
      - [2.1.1 (Option A) Start from Raw Matrix](#211-option-a-start-from-raw-matrix)
      - [2.1.2 (Option B) Start from AnnDataR6 Object](#212-option-b-start-from-anndatar6-object)
      - [2.1.3 (Optional) Filter Out Tumor Cells](#213-optional-filter-out-tumor-cells)
    - [2.2 Bulk expression data](#22-bulk-expression-data)
    - [2.3 Phenotype Data](#23-phenotype-data)
  - [3. Screening Cells Associated with Phenotype](#3-screening-cells-associated-with-phenotype)
    - [3.1 (Option A) Scissor Screening](#31-option-a-scissor-screening)
    - [3.2 (Option B) scPAS Screening](#32-option-b-scpas-screening)
    - [3.3 (Option C) scAB Screening](#33-option-c-scab-screening)
    - [3.4 (Option D) scPP Screening](#34-option-d-scpp-screening)
    - [3.5 (Optional) Merge screening results](#35-optional-merge-screening-results)
  - [4. Visualization](#4-visualization)
    - [4.1 UMAP for screening results](#41-umap-for-screening-results)
    - [4.2 Stack bar plot for screening results](#42-stack-bar-plot-for-screening-results)
  - [5. Example {#example}](#5-example-example)
  - [6. Other function details](#6-other-function-details)
  - [7. References](#7-references)

### 0.1 Introduction to SigBridgeR

SigBridgeR (short for **Sig**nificant cell-to-phenotype **Bridge** in **R**) is an R package for screening cells highly associated with phenotype data using single-cell RNA-seq, bulk expression and sample_related phenotype data at a pan-cancer level. It integrates functionality from these R packages: [Github-sunduanchen/Scissor](https://github.com/sunduanchen/Scissor), [Github-Qinran-Zhang/scAB](https://github.com/Qinran-Zhang/scAB/), [Github-WangX-Lab/ScPP](https://github.com/WangX-Lab/ScPP) and [Github-aiminXie/scPAS](https://github.com/aiminXie/scPAS).

------------------------------------------------------------------------

## 1. Installation

Install **SigBridgeR** using one of these methods:

### 1.1 Stable Release from GitHub

```{r install_from_github}
if(!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("WangLabCSU/SigBridgeR")
```

### 1.2 Check Dependencies

You can use this function to quickly verify installed dependencies and their versions:

```{r check_dependencies}
CheckPkgs <- function(packages) {
  CheckSinglePkg <- function(pkg_spec) {
    installed <- requireNamespace(pkg_spec$pkg, quietly = TRUE)
    current_version <- NA_character_
    
    if (installed) {
      current_version <- as.character(utils::packageVersion(pkg_spec$pkg))
      
      if (!is.null(pkg_spec$version)) {
        installed <- current_version >= pkg_spec$version
      }
    }
    
    data.frame(
      Package = pkg_spec$pkg,
      Required_Version = if (is.null(pkg_spec$version)) "Any" else pkg_spec$version,
      Installed = installed,
      Current_Version = if (installed) current_version else NA,
      stringsAsFactors = FALSE
    )
  }
  
  result <- do.call(rbind, lapply(packages, CheckSinglePkg))
  rownames(result) <- NULL  
    
  return(result)
}

CheckPkgs(list(
  list(pkg = "Seurat", version = "5.0.0"),
  list(pkg = "dplyr"),
  list(pkg = "cli"),
  list(pkg = "AUCell", version = "1.20.2"),
  list(pkg = "future"),
  list(pkg = "ggrepel"),
  list(pkg = "grid"),
  list(pkg = "IDConverter"),
  list(pkg = "Matrix"),
  list(pkg = "patchwork"),
  list(pkg = "scAB"),
  list(pkg = "Scissor"),
  list(pkg = "scPAS"),
  list(pkg = "ScPP"),
  list(pkg = "tibble"),
  list(pkg = "tidyr"),
  list(pkg = "data.table", version = "1.14.1")
  )
)

```

------------------------------------------------------------------------

## 2. Loading and preprocessing data

First load the package:

```{r load_package}
library(SigBridgeR)
```

### 2.1 Single-cell RNA-seq Data

You can use function `SCPreProcess` to preprocess your single-cell RNA-seq data. Here are some options:

#### 2.1.1 (Option A) Start from Raw Matrix

When starting from a raw count matrix, `SCPreProcess` will automatically perform Seurat preprocess (including `NormalizeData, FindVariableFeatures, ScaleData, RunPCA, FindNeighbors, FindClusters, RunTSNE, RunUMAP`). `SCPreProcess` returns a fully preprocessed Seurat object for downstream use.

```{r scpreprocessing_raw_matrix}
your_seurat <- SCPreProcess(
  your_matrix,
  project = "Scissor_Single_Cell", # Parameters used in Seurat preprocessing pipeline
  min_cells = 400,
  min_features = 0,
  normalization_method = "LogNormalize",
  scale_factor = 10000,
  selection_method = "vst",
  resolution = 0.6,
  dims = 1:10,
  verbose = TRUE,
  future_global_maxsize = 6 * 1024^3,
  ...
)
```

#### 2.1.2 (Option B) Start from AnnDataR6 Object

`SCPreProcess` also supports AnnData objects. You may reference and use the following code:

```{r scpreprocessing_anndata}
reticulate::use_pythonenv("The_path_to_your_python") 

anndata_obj <- anndata::read_h5ad("path_to_your_file.h5ad")

your_seurat <- SCPreProcess(
  anndata_obj,
  project = "Scissor_Single_Cell", # Parameters used in Seurat preprocessing pipeline
  min_cells = 400,
  min_features = 0,
  normalization_method = "LogNormalize",
  scale_factor = 10000,
  selection_method = "vst",
  resolution = 0.6,
  dims = 1:10,
  verbose = TRUE,
  future_global_maxsize = 6 * 1024^3,
  ...
)
```

The description of data in `anndata_obj$obs` will be add to `your_seurat@meta.data`.

**helpful documentation:**

-   [AnnData for R](https://github.com/dynverse/anndata)

#### 2.1.3 (Optional) Filter Out Tumor Cells

If you aim to filter out phenotype-associated cells in all tumor cells. With a preprocessed Seurat object (containing results from ), `SCPreProcess` will filter out tumor cells using the specified metadata column:

```{r scpreprocessing_seurat}
your_seurat <- SCPreProcess(your_seurat, column2only_tumor = "Tissue")
```

> Note: I don't recommend using columns like `column2only_tumor = "Celltype"` as tumor cell identities vary across tissues. instead:
>
> -   Create a Dedicated Column: Add a new metadata column (e.g., is_tumor) to explicitly label cells:"Tumo(u)r"/"Normal"
>
> -   Code Example:
>
> ```{r }
> # For glioblastoma (GBM)
> seurat_obj@meta.data$is_tumor <- ifelse(
>  grepl("GBM|glioblastoma|astrocytoma_grade_IV", seurat_obj@meta.data$Celltype, ignore.case = TRUE),
>  "Tumor",  # or "Tumour" 
>  "Normal"  # or "Non-Tumor" 
> )
> ```

### 2.2 Bulk expression data

`BulkPreProcess` performs a straightforward task: converting common gene identifiers (e.g., Ensembl IDs, Entrez) to standardized gene symbols by using the [IDConverter](https://github.com/ShixiangWang/IDConverter) package.

```{r bulk_preprocessing}
# genes * samples
your_bulk_data <- read.csv("path_to_your_file.csv", header = TRUE, row.names = 1)

your_bulk_data <- BulkPreProcess(your_bulk_data)
```

You can also use the `org.Hs.eg.db` package for gene symbol matching if you prefer not to use BulkPreProcess's built-in `IDConverter`.

```{r bulk_preprocessing2}
library(org.Hs.eg.db)

your_bulk_data <- read.csv("path_to_your_file.csv", header = TRUE, row.names = 1)

ensembl_ids <- sub("\\..*", "", rownames(your_bulk_data))
gene_symbols <- mapIds(org.Hs.eg.db, 
                      keys = ensembl_ids,
                      column = "SYMBOL",
                      keytype = "ENSEMBL",
                      multiVals = "first")

rownames(your_bulk_data) <- gene_symbols
your_bulk_data <- your_bulk_data[!is.na(rownames(your_bulk_data)), ]

```

### 2.3 Phenotype Data

bascially you can just use your phenotype data directly. If you are confused about the structure `Screen()` requires, please refer to [Section 5](#example).

------------------------------------------------------------------------

## 3. Screening Cells Associated with Phenotype

The function **`Screen`** provide 4 different options for screening cells associated with phenotype, These 4 algorithms come from the repositories mentioned in [Section 0.1](#01-introduction-to-sigbridger), and you can choose one of them to screen your cells.

Key parameters for `Screen`:

-   `matched_bulk`: A data frame of bulk expression data after intersecting samples.
-   `sc_data`: A Seurat object after preprocessing, you can use the output of `Preprocess` function or your own preprocessed Seurat object.
-   `phenotype`: A data frame of phenotype data after intersecting samples.
-   `label_type`: A character value specifying the filtering labels are stored in the `Seurat_object@misc` , default: `NULL`.
-   `phenotype_class`: A character value specifying the phenotype data type, i.e. `"binary"`, `"survival"` or `"continuous"`.
-   `screen_method`: A character value specifying the screening method, i.e. "Scissor", "scPAS", "scAB" or "scPP"
-   `...`: Other parameters for the screening methods.

### 3.1 (Option A) Scissor Screening

Parameters pass to `...` when using `Scissor` method:

-   `path2save_scissor_inputs`: A character value specifying the path to save intermediate data, default: `Scissor_inputs.RData`
-   `path2load_scissor_cahce`: A character value specifying the path to load intermediate data
-   `scissor_alpha`: Parameter used to balance the effect of the l1 norm and the network-based penalties. It can be a number or a searching vector. If alpha = NULL, a default searching vector is used. The range of alpha is in `[0,1]`. A larger alpha lays more emphasis on the l1 norm.
-   `scissor_cutoff`: Cutoff for the percentage of the Scissor selected cells in total cells. This parameter is used to restrict the number of the Scissor selected cells. A cutoff less than 50% (default 20%) is recommended depending on the input data.
-   `reliability_test_alpha`: the same as `scissor_alpha`
-   `reliability_test_n`: Permutation times (default: 10)
-   `nfold`: The fold number in cross-validation (default: 10)

**Usage**:

```{r scissor_screening}
scissor_result = Screen(
  matched_bulk = matched_bulk,
  sc_data = sc_dataset, # A Seurat object after preprocessing
  phenotype = matched_phenotype_data,
  label_type = "TP53", # The filtering labels are stored in the `@misc` 
  phenotype_class = "binary",  
  screen_method = c("Scissor"),
  path2save_scissor_inputs = "Tmp/Scissor_inputs.RData" # Intermediate data
)
```

You can use the intermediate data for repeated runs. This is an inherent feature of the `scissor`.

```{r scissor_screening_cache}
scissor_result = Screen(
  sc_data = sc_dataset, 
  label_type = "TP53", 
  phenotype_class = "binary", 
  screen_method = c("Scissor"),
  path2load_scissor_cahce = "Tmp/Scissor_inputs.RData" # Intermediate data
  path2save_scissor_inputs = "Tmp/Scissor_inputs.RData" 
)
```

If only the parameters `scissor_alpha` and `scissor_cutoff` are adjusted, this method can also be applied.

```{r scissor_screening_param_adjusted}
scissor_result = Screen(
  label_type = "TP53", 
  phenotype_class = "binary", 
  screen_method = c("Scissor"),
  path2load_scissor_cahce = "Tmp/Scissor_inputs.RData", # Intermediate data
  scissor_alpha = 0.05, 
  scissor_cutoff = 0.05 
)

```

**Returning structure**: A list containing:

-   `scRNA_data`: A Seurat object after screening
-   `scissor_result`: The result of Scissor screening
-   `reliability_test`: Reliability test results

**Cell level Evaluation**:

You can use `Sissor::evaluate.cell()` to obtain some supporting information for each Scissor selected cell. First, prepare a benchmark data set.

```{r scissor_screening_cell_level_evaluation}
evaluate_summary <- Scissor::evaluate.cell(
    './data/benchmark_data.RData',
    scissor_result$scissor_result,
    FDR = 0.05,
    bootstrap_n = 100
)
```

helpful documentation:

[Scissor-Cell Level Evaluations](https://sunduanchen.github.io/Scissor/vignettes/Scissor_Tutorial.html#cell-level-evaluations)

### 3.2 (Option B) scPAS Screening

Parameters pass to `...` when using `scPAS` method (basically adapted from the `scPAS`'s documentation):

-   Parameters passed to `scPAS::scPAS()`

    These parameters directly interface with the core `scPAS`() function from the original package:

    -   `assay`: Name of Assay to get.
    -   `imputation`: Logical. imputation or not.
    -   `nfeature`: Numeric. The Number of features to select as top variable features in `sc_data`. Top variable features will be used to intersect with the features of `matched_bulk`. Default is NULL and all features will be used.
    -   `alpha`: Numeric. Parameter used to balance the effect of the l1 norm and the network-based penalties. It can be a number or a searching vector. If `alpha = NULL`, a default searching vector is used. The range of alpha is in `[0,1]`. A larger alpha lays more emphasis on the l1 norm.
    -   `network_class`: The source of feature-feature similarity network. By default this is set to sc and the other one is bulk.

-   Parameters passed to `SigBridgeR::DoscPAS`

    These parameters represent customizable options provided by upstream/downstream processing functions when integrating the `scPAS` workflow

    -   `extra_filter`: A logical value. If `TRUE`, the following options for filtering will be applied. Default is `FALSE`.
    -   `gene_RNAcount_filter`: A numeric value. The threshold for filtering out genes with low RNA counts. Default is `20`.
    -   `bulk_0_filter_thresh`: A numeric value. The threshold for filtering out samples with low RNA counts. Default is `0.05`.

**usage**:

```{r scPAS_screening}
scpas_result = Screen(
  matched_bulk = matched_bulk,
  sc_data = A_Seurat_object,
  phenotype = phenotype,
  label_type = "TP53", # The filtering labels are stored in the `@misc` 
  screen_method = "scpas",
  phenotype_class = "binary"
)
```

**returning structure**: A list containing:

-   `scRNA_data`: A Seurat object after screening

### 3.3 (Option C) scAB Screening

Parameters pass to `...` when using `scAB` method (basically adapted from the `scAB`'s documentation):

-   `alpha`: Coefficient of phenotype regularization, default is `0.005`
-   `alpha_2`: Coefficient of cell-cell similarity regularization, default is `5e-05`
-   `maxiter`: Maximum number of iterations, default is `2000`
-   `tred`: Threshold for early stopping, default is `2`

**usage**:

```{r scAB_screening}
scab_result = Screen(
  matched_bulk = your_matched_bulk,
  sc_data = A_Seurat_object,
  phenotype = your_matched_phenotype,
  label_type = "TP53", # The filtering labels are stored in the `@misc` 
  screen_method = "scAB",
  phenotype_class = "binary",
)
```

**returning structure**: A list containing:

-   `scRNA_data`: A Seurat object after screening
-   `scAB_result`: A list with the submatrix and loss value

### 3.4 (Option D) scPP Screening

Parameters pass to `...` when using `scPP` method :

-   `ref_group`: The reference group for the binary analysis, default is `1`
-   `Log2FC_cutoff`: The cutoff for the log2 fold change of the binary analysis, default is `0.585`
-   `estimate_cutoff`: Effect size threshold for continuous traits, default is `0.2`
-   `probs`: Quantile cutoff for cell classification, default is `0.2`

**usage**:

```{r scPP_screening}
# This will take several hours
scpp_result = Screen(
  matched_bulk = your_matched_bulk,
  sc_data = A_Seurat_object,
  phenotype = your_matched_phenotype,
  label_type = "TP53", # The filtering labels are stored in the `@misc` 
  screen_method = "scpp",
  phenotype_class = "binary",
)
```

**returning structure**: A list containing:

-   `scRNA_data`: A Seurat object after screening

### 3.5 (Optional) Merge screening results

If you have performed multiple screening methods one the same data, you can use the function `MergeResult` to merge the results of these methods. The Seurat object or a results list from `Screen` is accepted.

```{r merge_screening_results}
merged_seurat=MergeResult(
    your_scissor_result, 
    your_scPAS_result, 
    your_scAB_result, 
    your_scPP_result
)

# * mixed input form is alse supported 

merged_seurat=MergeResult(
    your_scissor_result$scRNA_data, 
    your_scPAS_result$scRNA_data, 
    your_scAB_result, 
    your_scPP_result,
)

```

This function performs a simple task of consolidating screening results into the Seurat object's `meta.data` slot and returning the updated object. Please note that the intermediate data (e.g., `scissor_result$reliability.test` and `scab_result$scAB_result`) will not be preserved in this process.

------------------------------------------------------------------------

## 4. Visualization

After screening, you can use **`ScreenFractionPlot`** and **`DimPlot`** from `Seurat` for plotting the results of screening

### 4.1 UMAP for screening results

**example**:

Suppose you have performed `scissor` algorithm screening on your Seurat object and wish to examine the distribution across different celltypes and patient, you may reference and use the following code:

```{r umap_exmaple}
library(patchwork)
library(zeallot)

c(celltype_umap, patient_umap, scissor_umap) %<-%
    purrr::map(
        c("celltype", "patient", "scissor"),
        ~ Seurat::DimPlot(
            your_seurat_obj,
            group.by = .x,
            label.size = 6,
            pt.size = 0.6,
            reduction = "umap"
        ) +
            ggplot2::ggtitle(.x)
    )

# * Show
celltype_umap | patient_umap | scissor_umap
```

This will generate three UMAP plots separately.

Or suppose you have performed `scPAS` screening on your Seurat object and want to visualize the distribution of prediction confidence scores, you may reference and use the following code:

```{r umap_exmaple2}
library(patchwork)
library(zeallot)

c(scPAS_Pvalue_umap, scPAS_NRS_umap) %<-%
    purrr::map(
        c("scPAS_Pvalue", "scPAS_NRS"),
        ~ Seurat::FeaturePlot(
            object = your_seurat_obj,
            features = .x,
            label.size = 6,
            pt.size = 0.6,
            label = TRUE
        ) +
            ggplot2::ggtitle(.x) +
            theme(legend.position = "right")
    )

# * Show
scPAS_Pvalue_umap | scPAS_NRS_umap
```

This will generate two plots, one for each feature specified in `feature`.

**helpful documentation**:

[Browser-Seurat::DimPlot](https://satijalab.org/seurat/reference/DimPlot.html)

[Browser-Seurat::FeaturePlot](https://satijalab.org/seurat/reference/FeaturePlot.html)

### 4.2 Stack bar plot for screening results

Key parameters for `ScreenFractionPlot`:

-   `seurat_obj`: A Seurat object after screening.
-   `group_by`: Used to specify the column of the meta.data in `seurat_obj`. The plot results will be grouped by this parameter.
-   `screen_type`: Screening algorithm used before. (case-sensitive, e.g., "scissor" for Scissor results)
-   `show_null`: Logical whether to show groups with zero cells (default: FALSE).
-   `plot_color` Custom color palette (named vector format):
    -   Required names: "Positive", "Negative", "Neutral"
    -   Default: c("Neutral"="#CECECE", "Positive"="#ff3333", "Negative"="#386c9b")

Suppose you have already performed the `scPAS` algorithm screening on your Seurat object, and you want to check the proportion of positive cells across different patients. You can refer to and use the following code:

```{r stack_bar_plot_example}
# based on `ggplot2`
plot <- ScreenFractionPlot(
   screened_seurat = your_seurat_obj,
   group_by = "patient", 
   screen_type = "scPAS",
   plot_title = "scPAS Screening Results"
)
```

The order of the groups is determined by the proportion of **Positive** cells within each group.

Use `?ScreenFractionPlot` in R to see more details.

------------------------------------------------------------------------

## 5. Example {#example}

Here we use the example data to demonstrate how to use the functions in `SigBridgeR` to screen cells associated with phenotype.

```{r example}
library(SigBridgeR)

```

------------------------------------------------------------------------

## 6. Other function details

-   `AddMisc()` : Add miscellaneous information to the Seurat object. Support for adding multiple attributes to the `SeuratObject@misc` slot simultaneously.

```{r add_misc_example}
# basic usage
seurat_obj <- AddMisc(seurat_obj, "QC_stats" = qc_df)

# Auto-incrementing example when `cover` set to FALSE
seurat_obj <- AddMisc(seurat_obj, markers = markers1)
seurat_obj <- AddMisc(seurat_obj, markers = markers2, cover=FALSE)

# Add multiple attributes to the `SeuratObject@misc` slot simultaneously
seurat_obj <- AddMisc(seurat_obj, markers1 = markers1, markers2 = markers2)
```

Use `?AddMisc` in R to see more details.

------------------------------------------------------------------------

## 7. References

Identifying phenotype-associated subpopulations by integrating bulk and single-cell sequencing data. Nature Biotechnology (2021). https://doi.org/10.1038/s41587-021-01091-3. Duanchen Sun, Xiangnan Guan, Amy E. Moran, Ling-Yun Wu, David Z. Qian, Pepper Schedin, Mu-Shui Dai, Alexey V. Danilov, Joshi J. Alumkal, Andrew C. Adey, Paul T. Spellman and Zheng Xia