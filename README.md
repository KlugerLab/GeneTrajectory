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

The standard preprocessing can be done by employing the [Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) R package which includes: library normalization; finding variable features; scaling; generating PCA embedding (and UMAP embedding for visualization).

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

Next, we construct the cell-cell kNN graph and calculate cell-cell graph distances.
```r
data_S <- GeneTrajectory::RunDM(data_S)
cell.graph.dist <- GetGraphDistance(data_S, K = 10, dims = 1:10)
```

We then narrow down the gene list for gene-gene distance computation by focusing on the top 2000 variable genes expressed by 1% - 50% of cells. 

```r
assay <- "RNA"
DefaultAssay(data_S) <- assay
data_S <- FindVariableFeatures(data_S, nfeatures = 2000)
all_genes <- data_S@assays[[assay]]@var.features
expr_percent <- apply(as.matrix(data_S[[assay]]@data[all_genes, ]) > 0, 1, sum)/ncol(data_S)
genes <- all_genes[which(expr_percent > 0.01 & expr_percent < 0.5)]
length(genes)
```

The intermediate outputs are written into a local directory which allows for gene-gene Wasserstein distance computation using Python OT package.
```r
cg_output <- CoarseGrain(data_S, cell.graph.dist, genes, N = 1000, dims = 1:10)
dir.path <- "./mouse_dermal/"
dir.create(dir.path, recursive=T)
write.table(cg_output[["features"]], paste0(dir.path, "gene_names.csv"), row.names = F, col.names = F, sep = ",")
write.table(cg_output[["graph.dist"]], paste0(dir.path, "ot_cost.csv"), row.names = F, col.names = F, sep = ",")
Matrix::writeMM(Matrix(cg_output[["gene.expression"]], sparse = T), paste0(dir.path, "gene_expression.mtx"))

system(sprintf("nohup /data/anaconda3/bin/python ./GeneTrajectory/python/gene_distance_cal_parallel.py %s &", dir.path))

```








