# Gene Trajectory Inference

GeneTrajectory is a method for inferring gene trajectories in scRNA-seq data, which faciliates understanding of gene dynamics underlying biological processes. The major workflow of GeneTrajectory comprises the following four main steps:
* Step 1. Build a cell-cell kNN graph in which each cell is connected to its k-nearest neighbors. Find the shortest path connecting each pair of cells in the graph and denote its length as the graph distance between cells.
* Step 2. Compute pairwise graph-based Wasserstein distance between gene distributions, which quantifies the minimum cost of transporting the distribution of a given gene into the distribution of another gene in the cell graph.
* Step 3. Generate a low-dimensional representation of genes (using Diffusion Map by default) based on the gene-gene Wasserstein distance matrix. Identify gene trajectories in a sequential manner.
* Step 4. Determine the order of genes along each gene trajectory.

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
#Load required packages
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
```

The intermediate outputs are written into a local directory which allows for gene-gene Wasserstein distance computation using [Python OT package](https://pythonot.github.io/).
```r
cg_output <- CoarseGrain(data_S, cell.graph.dist, genes, N = 1000, dims = 1:10)
dir.path <- "./mouse_dermal/" #This should be replaced by your local directory path.
dir.create(dir.path, recursive=T)
write.table(cg_output[["features"]], paste0(dir.path, "gene_names.csv"), row.names = F, col.names = F, sep = ",")
write.table(cg_output[["graph.dist"]], paste0(dir.path, "ot_cost.csv"), row.names = F, col.names = F, sep = ",")
Matrix::writeMM(Matrix(cg_output[["gene.expression"]], sparse = T), paste0(dir.path, "gene_expression.mtx"))

#This can also be run in the terminal. Please make sure to install the latest version of POT module (python), using the following:
#pip install -U https://github.com/PythonOT/POT/archive/master.zip
system(sprintf("nohup /data/anaconda3/bin/python ./GeneTrajectory/python/gene_distance_cal_parallel_Iter50000.py %s &", dir.path)) #Python paths should be adjusted accordingly.
```

When the computation is finished, a file named `emd.csv` is generated under the same directory. Gene trajectory extraction and visualization can be done using the following code.

```r
gene.dist.mat <- LoadGeneDistMat(dir.path, file_name = "emd.csv")
gene_embedding <- GetGeneEmbedding(gene.dist.mat, K = 5)$diffu.emb
gene_trajectory <- ExtractGeneTrajectory(gene_embedding, gene.dist.mat, N = 3, t.list = c(9,16,5), K = 5)

#Visualizing gene trajectories using the leading three non-trivial eigenvectors.
par(mar = c(1.5,1.5,1.5,1.5))
scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3],
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 1, theta = 45, phi = 0,
          col = ramp.col(c(hue_pal()(3))))
```

To examine how each given gene trajectory is reflected over the cell graph, we can track how these genes are expressed across different regions in the cell embedding. Here, we would recommend users to apply [ALRA](https://github.com/KlugerLab/ALRA/blob/master/README.md) imputation to smooth the expression values for generating gene bin plots.
```r
N.bin = 7
data_S <- RunALRA(data_S)
data_S <- AddGeneBinScore(data_S, gene_trajectory, N.bin = N.bin, trajectories = 1:3, assay = "alra")
FeaturePlot(data_S, pt.size = 0.05, features = paste0("Trajectory",1,"_genes", 1:N.bin), ncol = N.bin, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes()
FeaturePlot(data_S, pt.size = 0.05, features = paste0("Trajectory",2,"_genes", 1:N.bin), ncol = N.bin, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes()
FeaturePlot(data_S, pt.size = 0.05, features = paste0("Trajectory",3,"_genes", 1:N.bin), ncol = N.bin, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes()
```





