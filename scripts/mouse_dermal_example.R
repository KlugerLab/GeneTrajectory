require(GeneTrajectory)
require(plot3D)
require(ggplot2)
require(viridis)
require(scales)
require(SeuratWrappers)


################################################
##### Standard preprocessing
################################################
data_S <- readRDS("./data/E14.5_mouse_dermal_seurat.rds")
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
data_S$cell_type <- "UD"
data_S$cell_type[which(data_S$cluster %in% c(0,4))] <- "LD"
data_S$cell_type[which(data_S$cluster %in% c(7))] <- "DC"
DimPlot(data_S, reduction = "umap", label = T, label.size = 5, group.by = "cell_type")
cbPalette <- c("#E69F00", "#56B4E9", "#009E73")
DimPlot(data_S, group.by = "cell_type", cols = cbPalette, split.by = "orig.ident") & NoAxes()
FeaturePlot(data_S, features = c("Dkk1", "Dkk2", "Sox2", "Lef1"), ncol = 4) & NoAxes() & scale_color_viridis()
DimPlot(data_S, group.by = "orig.ident") & NoAxes()


################################################
##### Prepare for gene-gene distance computation
################################################
##### Select genes
assay <- "RNA"
DefaultAssay(data_S) <- assay
data_S <- FindVariableFeatures(data_S, nfeatures = 2000)
all_genes <- data_S@assays[[assay]]@var.features
expr_percent <- apply(as.matrix(data_S[[assay]]@data[all_genes, ]) > 0, 1, sum)/ncol(data_S)
genes <- all_genes[which(expr_percent > 0.01 & expr_percent < 0.5)]
length(genes)

##### Prepare the input for gene-gene Wasserstein distance computation
data_S <- GeneTrajectory::RunDM(data_S)
cell.graph.dist <- GetGraphDistance(data_S, K = 10, dims = 1:10)
cg_output <- CoarseGrain(data_S, cell.graph.dist, genes, N = 1000, dims = 1:10)

#####Save the files for computation in Python
dir.path <- "./data/mouse_dermal/N1000/"
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
# python gene_distance_cal.py ./data/mouse_dermal/N1000/



################################################
##### Gene trajectory inference and visualization
################################################
gene.dist.mat <- LoadGeneDistMat(dir.path, file_name = "emd.csv")
gene_embedding <- GetGeneEmbedding(gene.dist.mat, K = 5)$diffu.emb
gene_trajectory <- ExtractGeneTrajectory(gene_embedding, gene.dist.mat, N = 3, t.list = c(9,16,4), K = 5)
table(gene_trajectory$selected)


par(mar = c(1.5,1.5,1.5,1.5))
scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3],
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 1, theta = 45, phi = 0,
          col = ramp.col(c(hue_pal()(3))))

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
data_S <- AddGeneBinScore(data_S, gene_trajectory, N.bin = 7, trajectories = 1:3, assay = "alra")
FeaturePlot(data_S, pt.size = 0.05, features = paste0("Trajectory",1,"_genes", 1:5), ncol = 5, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes()
FeaturePlot(data_S, pt.size = 0.05, features = paste0("Trajectory",2,"_genes", 1:5), ncol = 5, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes()
FeaturePlot(data_S, pt.size = 0.05, features = paste0("Trajectory",3,"_genes", 1:5), ncol = 5, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes()


