# SIGNAL: Spatial varIable Gene detectioN And cLassification
![SIGNAL workflow](man/figures/SIGNAL_Fig1_overview.pdf)

SIGNAL is a statistical method to detect spatially variable genes in single-cell resolution spatial transcriptomics and classify them by the source of their spatial variation. For each gene, SIGNAL fits a multi-kernel variance component model with a cell-type kernel, a spatial kernel and their interaction. An overall test detects any effect. Separate tests then classify each detected gene as having a **cell-type** effect (expression differs between cell types), a **spatial** effect (expression varies over space for all cell types) and/or a **cell-type-specific spatial** effect (ctSVG). A second-stage test finds the cell types that carry each ctSVG pattern. Results are combined across a panel of spatial kernels with the Cauchy combination test.

**Reference**: Hou, T. *et al.* SIGNAL: Spatial varIable Gene detectioN And cLassification. (in preparation)

## Installation
Please run the following codes in R to install SIGNAL package from GitHub.
```
if (!require("devtools", quietly = TRUE)){
  install.packages("devtools")
}
devtools::install_github("TingtHou/SIGNAL", build_vignettes = TRUE)
```

SIGNAL depends only on `CompQuadForm`, `ggplot2` and packages that ship with R (`mgcv`, `methods`, `parallel`, `stats`, `utils`). It does not need compilation.

## Detect and classify spatially variable genes with SIGNAL
The best vignette for getting started with SIGNAL is the [Tutorial](https://tingthou.github.io/SIGNAL/tutorial.html). After installation, you can also open it in R with `vignette("SIGNAL")`.

```
library(SIGNAL)
data(SIGNAL_example)

obj <- createSIGNALobject(counts   = SIGNAL_example$counts,    # genes x cells
                          location = SIGNAL_example$location,  # cells x 2
                          celltype = SIGNAL_example$celltype)  # one label per cell
obj <- SIGNAL_preprocess(obj)            # filter genes, rank-based normal transform
obj <- SIGNAL_kernel(obj)                 # cell-type kernel + 11 spatial kernels
obj <- SIGNAL_test(obj)                   # overall, cell-type, spatial and ctSVG tests
obj <- SIGNAL_celltype(obj)               # which cell types carry each ctSVG

summary(obj)              # number of genes in each category
head(obj@result)          # per-gene p-values and classification
obj@celltype_result       # ctSVG genes x cell types p-values
SIGNAL_plotGene(obj, "gene61", by_celltype = TRUE)
```

`SIGNAL()` runs the first four steps in one call.

The defaults reproduce the original SIGNAL analysis code exactly:

* `transform = "rank"`: each gene is converted to normal scores, `qnorm((rank(y) - 0.5) / n)`, with tied counts sharing their average rank. Results are deterministic. `transform = "impute"` breaks tied zeros at random instead (set `seed`).
* `method = "davies_liu"`: the overall test uses Davies' method and the cell-type, spatial and ctSVG tests use Liu's moment-matching approximation. `method = "davies"` uses Davies for all four tests, and `method = "moment"` uses Liu for all four (fastest).

The overall test needs an eigendecomposition per gene, so for large data sets split the genes into batches, for example one HPC job per batch.

| Column in `obj@result` | Description |
|---|---|
| `p_overall` | test of any effect (cell type, spatial or ctSVG) |
| `p_celltype`, `p_spatial`, `p_ctSVG` | test of each effect on its own |
| `padj_overall` | `p_overall` adjusted for multiple testing (BH by default) |
| `celltype_effect`, `spatial_effect`, `ctSVG_effect` | gene is significant overall and that effect has p < `alpha` |
| `category` | e.g. `"Spatial"`, `"CellType+ctSVG"`, `"Unclassified"`, `"Not significant"` |

The default spatial kernel panel was adapted from CELINA (Shang, Wu & Zhou, *Nat. Commun.* 2025). SIGNAL computes the kernels itself and does not need CELINA.
