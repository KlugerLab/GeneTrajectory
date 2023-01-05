require(GeneTrajectory)
require(plot3D)
require(ggplot2)
require(viridis)
require(scales)
require(SeuratWrappers)


################################################
##### Standard preprocessing
################################################
data_S <- readRDS("./data/pbmc_10k_myeloid_seurat.rds")
DefaultAssay(data_S) <- "RNA"
data_S <- NormalizeData(data_S)
data_S <- FindVariableFeatures(data_S)
data_S <- ScaleData(data_S)
data_S <- RunPCA(data_S, npcs = 30, verbose = F)
data_S <- RunUMAP(data_S, dims = 1:30)
data_S <- FindNeighbors(data_S, dims = 1:30)
data_S <- FindClusters(data_S, resolution = 0.3)
data_S$cluster <- data_S$RNA_snn_res.0.3
DimPlot(data_S, reduction = "umap", label = T, label.size = 5, group.by = "cluster")

################################################
##### Annotate cell types
################################################
DefaultAssay(data_S) <- "RNA"
Idents(data_S) <- "cluster"
cluster_markers <- FindAllMarkers(
  data_S, only.pos = T, min.diff.pct = 0.1
)
cluster_markers$pct.diff <- cluster_markers$pct.1 - cluster_markers$pct.2
cluster_markers <- cluster_markers[which(cluster_markers$p_val_adj <= 0.05), ]
top8 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 8, wt = avg_log2FC)
DoHeatmap(data_S, features = top8$gene, slot = "scale.data") + theme(axis.text.y = element_text(size = 10)) + scale_fill_viridis()
DotPlot(data_S, features = unique(top8$gene), group.by = "cluster") & RotatedAxis()
cluster_relabel <- c("0" = "CD14+CD16- monocytes",
                     "1" = "CD14+HLA-DR(high) monocytes",
                     "2" = "CD14-CD16+ monocytes",
                     "3" = "CD14- dendritic cells (type-2)")
data_S$celltype_new <- cluster_relabel[as.character(data_S$cluster)]
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
DimPlot(data_S, group.by = "celltype_new", cols = cbPalette) & NoAxes()
FeaturePlot(data_S, features = c("CD14", "HLA-DMA", "FCGR3A", "CD1C"), ncol = 4) & NoAxes() & scale_color_viridis()

Idents(data_S) <- "celltype_new"
DoHeatmap(data_S, features = top8$gene, group.colors = cbPalette[c(3,4,2,1)], slot = "scale.data", size = 2.5) + theme(axis.text.y = element_text(size = 10)) + scale_fill_viridis()
DotPlot(data_S, features = unique(top8$gene), group.by = "celltype_new") & RotatedAxis()


################################################
##### Prepare for gene-gene distance computation
################################################
##### Select genes
assay <- "RNA"
DefaultAssay(data_S) <- assay
data_S <- FindVariableFeatures(data_S, nfeatures = 2000)
all_genes <- data_S@assays[[assay]]@var.features
expr_percent <- apply(as.matrix(data_S[[assay]]@data[all_genes, ]) > 0, 1, sum)/ncol(data_S)
genes <- all_genes[which(expr_percent > 0.005 & expr_percent < 0.75)]
length(genes)

##### Prepare the input for gene-gene Wasserstein distance computation
data_S <- GeneTrajectory::RunDM(data_S)
cell.graph.dist <- GetGraphDistance(data_S, K = 10)
cg_output <- CoarseGrain(data_S, cell.graph.dist, genes, N = 1000)

#####Save the files for computation in Python
dir.path <- "./data/human_myeloid/N1000/"
dir.create(dir.path)
write.table(cg_output[["features"]], paste0(dir.path, "gene_names.csv"), row.names = F, col.names = F, sep = ",")
write.table(cg_output[["graph.dist"]], paste0(dir.path, "ot_cost.csv"), row.names = F, col.names = F, sep = ",")
Matrix::writeMM(Matrix(cg_output[["gene.expression"]], sparse = T), paste0(dir.path, "gene_expression.mtx"))



################################################
##### Wasserstein distance computation in Python
################################################
# Make sure to install the latest version of POT module (python), using the following:
# pip install -U https://github.com/PythonOT/POT/archive/master.zip
# Run the following command to get the gene-gene Wasserstein distance matrix
# python gene_distance_cal.py ./data/human_myeloid/N1000/



################################################
##### Gene trajectory inference and visualization
################################################
gene.dist.mat <- LoadGeneDistMat(dir.path, file_name = "emd.csv")
gene_embedding <- GetGeneEmbedding(gene.dist.mat, K = 5)$diffu.emb
gene_trajectory <- ExtractGeneTrajectory(gene_embedding, gene.dist.mat, N = 3, t.list = c(11,21,8), K = 5)
table(gene_trajectory$selected)


gene_labels <- paste("-----", rownames(gene_embedding))
names(gene_labels) <- rownames(gene_embedding)
genes_selected <- c("CCR2", "ICAM2", "FCGR3A",  "SELL", "C1QA", "C1QB",
                    "CD1C", "CLEC10A", "CD2", "CD72", "CCR5",
                    "PKIB", "RETN", "CLEC5A", "CSF1R")
gene_labels[names(gene_labels) %notin% genes_selected] <- ""
par(mar = c(1.5,1.5,1.5,1.5))
scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3],
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 1, theta = 45, phi = 0,
          col = ramp.col(c(hue_pal()(3))))
text3D(gene_embedding[,1],
       gene_embedding[,2],
       gene_embedding[,3],  labels = gene_labels,
       add = T, colkey = FALSE, cex = 0.5)

scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3], alpha = 1,
          bty = "b2", colvar = gene_trajectory$Pseudoorder1,
          main = "", pch = 19, cex = 1, theta = 45, phi = 0,
          col = ramp.col(rev(viridis(12)[2:11])))
scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3], alpha = 1,
          bty = "b2", colvar = gene_trajectory$Pseudoorder2,
          main = "", pch = 19, cex = 1, theta = 45, phi = 0,
          col = ramp.col(rev(viridis(12)[2:11])))
scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3], alpha = 1,
          bty = "b2", colvar = gene_trajectory$Pseudoorder3,
          main = "", pch = 19, cex = 1, theta = 45, phi = 0,
          col = ramp.col(rev(viridis(12)[2:11])))

gene_list <- list()
for (i in 1:3){
  message(paste0("Trajectory-", i))
  gene_trajectory_sub <- gene_trajectory[which(gene_trajectory$selected == paste0("Trajectory-", i)),]
  genes <- rownames(gene_trajectory_sub)[order(gene_trajectory_sub[, paste0("Pseudoorder", i)])]
  message(paste(rev(genes), collapse = ", "))
  gene_list[[i]] <- genes
}


data_S <- RunALRA(data_S)
data_S <- AddGeneBinScore(data_S, gene_trajectory, N.bin = 5, trajectories = 1:3, assay = "alra", reverse = c(F, F, T))
FeaturePlot(data_S, pt.size = 0.05, features = paste0("Trajectory",1,"_genes", 1:5), ncol = 5, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes()
FeaturePlot(data_S, pt.size = 0.05, features = paste0("Trajectory",2,"_genes", 1:5), ncol = 5, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes()
FeaturePlot(data_S, pt.size = 0.05, features = paste0("Trajectory",3,"_genes", 1:5), ncol = 5, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes()


