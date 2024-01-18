# Gene Trajectory Inference

This is the first version of GeneTrajectory that performs gene trajectory inference in scRNA-seq data. The example notebooks and scripts can be found under Rmd_notebooks, benchmark, and scripts.
More example vignettes will be added soon. Processed example data can be provided upon request.

![Workflow](https://github.com/RihaoQu/IGT/blob/master/images/GT_workflow.png)

## Install

GeneTrajectory can be installed as follows:

```r
install.packages("devtools")
devtools::install_github("RihaoQu/IGT")
```
## Example tutorial

The standard preprocessing can be done by using the Seurat R package which includes: library normalization; finding variable features; scaling; PCA embedding; and UMAP embedding.

```r
require(GeneTrajectory)
require(plot3D)
require(ggplot2)
require(viridis)
require(scales)
require(Seurat)
require(SeuratWrappers)

# Here, we load a preprocessed Seurat object for gene trajectory inference.
data_S <- data_S_WLS_combined_E14.5
DimPlot(data_S, reduction = "umap", label = T, label.size = 5, group.by = "cell_type") & NoAxes()
DimPlot(data_S, group.by = "orig.ident") & NoAxes()
```

