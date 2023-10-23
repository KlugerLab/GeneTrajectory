require(Seurat)
require(GeneTrajectory)
require(ggplot2)
require(mgcv)

##### Load Seurat object
data_S <- readRDS("./data/E14.5_mouse_dermal_seurat.rds")

##### Find the root cell for the CTI-based methods
sub_data_S <- subset(data_S, cell_type %in% c("DC", "UD"))
data.matrix <- sub_data_S@reductions$pca@cell.embeddings
dist.mat <- sapply(rownames(data.matrix),function(cell_name) sqrt(sum((data.matrix[cell_name,] - data.matrix[root_cell,])**2)))
sub_data_S@meta.data$eu_dist <- dist.mat
sub_data_S@meta.data$pseudo_orders <- rank(-sub_data_S@meta.data$eu_dist)
FeaturePlot(sub_data_S,reduction = "umap",features = "pseudo_orders") & scale_color_viridis_c()
sub_data_S$root <- 0
sub_data_S$root[which(sub_data_S$cluster == 1)][which.min(sub_data_S@meta.data$pseudo_orders[which(sub_data_S$cluster == 1)])] <- 1
FeaturePlot(sub_data_S,reduction = "umap",features = "root", order = T) 

##### Key gene markers
genes_to_show <- c("Dkk1","Lef1","Sox2","Gli1","Ptch1","Grem1","Cdkn1a","Sox18","Foxd1","Bmp4")
pathways <- c("Wnt","Wnt","DC","SHH","SHH","Wnt","Cdkn1a","DC","DC","Wnt")
names(pathways) <- genes_to_show


##### GeneTrajectory
dir.path <- paste0("./data/mouse_dermal/N1000/")
gene.dist.mat <- LoadGeneDistMat(dir.path, file_name = "emd.csv")
gene_embedding <- GetGeneEmbedding(gene.dist.mat, K = 5)$diffu.emb
gene_trajectory <- ExtractGeneTrajectory(gene_embedding, gene.dist.mat, N = 3, t.list = c(9,16,5), K = 5)
par(mar = c(1.5,1.5,1.5,1.5))
scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3],
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 1, theta = 45, phi = 0,
          col = ramp.col(c(hue_pal()(3))))
gene_list <- list()
for (i in 1:3){
  message(paste0("Trajectory-", i))
  gene_trajectory_sub <- gene_trajectory[which(gene_trajectory$selected == paste0("Trajectory-", i)),]
  genes <- rownames(gene_trajectory_sub)[order(gene_trajectory_sub[, paste0("Pseudoorder", i)])]
  message(paste(genes, collapse = ", "))
  gene_list[[i]] <- genes
}
E14.5_DC_genes <- gene_list[[1]]

dir.path <- paste0("./data/mouse_dermal_noCC/N1000/")
gene.dist.mat <- LoadGeneDistMat(dir.path, file_name = "emd.csv")
gene_embedding <- GetGeneEmbedding(gene.dist.mat, K = 5)$diffu.emb
gene_trajectory <- ExtractGeneTrajectory(gene_embedding, gene.dist.mat, N = 3, t.list = c(9,16,5), K = 5)
par(mar = c(1.5,1.5,1.5,1.5))
scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3],
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 1, theta = 45, phi = 0,
          col = ramp.col(c(hue_pal()(3))))
gene_list <- list()
for (i in 1:3){
  message(paste0("Trajectory-", i))
  gene_trajectory_sub <- gene_trajectory[which(gene_trajectory$selected == paste0("Trajectory-", i)),]
  genes <- rownames(gene_trajectory_sub)[order(gene_trajectory_sub[, paste0("Pseudoorder", i)])]
  message(paste(genes, collapse = ", "))
  gene_list[[i]] <- genes
}
E14.5_DC_genes_noCC <- gene_list[[1]]

##### Find common genes
shared_genes <- intersect(E14.5_DC_genes_noCC, E14.5_DC_genes)

E14.5_DC_genes_idx <- 1:length(shared_genes)
names(E14.5_DC_genes_idx) <- E14.5_DC_genes
E14.5_DC_genes_noCC_idx <- 1:length(shared_genes)
names(E14.5_DC_genes_noCC_idx) <- E14.5_DC_genes_noCC[which(E14.5_DC_genes_noCC %in% shared_genes)]


gene_order_vec <- 1:length(gene_list[[1]])
names(gene_order_vec) <- gene_list[[1]]
GT_gene_order_vec <- gene_order_vec
GT_gene_order_vec_noCC <- E14.5_DC_genes_noCC_idx[names(GT_gene_order_vec)]
GT_gene_order_vec_noCC 



##### CellRank
cellrank_gene_order <- read.csv("./data/mouse_dermal_E14_WLS_pairs_noLD/cellrank_gene_order.csv", header = F)
cellrank_gene_order$V1[which(cellrank_gene_order$V1 == 'Pakap')] <- 'Pakap.1'
cellrank_gene_order_vec <- 1:length(cellrank_gene_order$V1)
names(cellrank_gene_order_vec) <- cellrank_gene_order$V1

cellrank_gene_order_noCC <- read.csv("./data/mouse_dermal_E14_WLS_pairs_noLD/cellrank_gene_order_noCC.csv", header = F)
cellrank_gene_order_noCC$V1[which(cellrank_gene_order_noCC$V1 == 'Pakap')] <- 'Pakap.1'
cellrank_gene_order_vec_noCC <- 1:length(cellrank_gene_order_noCC$V1)
names(cellrank_gene_order_vec_noCC) <- cellrank_gene_order_noCC$V1


##### Other methods
setwd("./data/mouse_dermal_E14_WLS_pairs_noLD/")
gene_expression_mat <- data_S[["RNA"]]@data[E14.5_DC_genes,]
methods <- c("PAGA","Monocle2","Monocle3","Slingshot")

gene_order_vec_list <- list()
for (method in methods){
  cell_pt <- read.csv(sprintf("./data/mouse_dermal_E14_WLS_pairs_noLD/mouse_dermal_E14_WLS_pairs_%s.csv", method), header = F)
  cells <- which(!is.na(cell_pt[,1]))
  length(cells)
  peak_time <- c()
  for (i in 1:nrow(gene_expression_mat)){
    message(i)
    dat <- data.frame(x=cell_pt[cells,1], y = gene_expression_mat[i,cells])
    mod_gam = gam(y ~ s(x), data = dat)
    newd <- data.frame(x = c(0:(1000*max(cell_pt[cells,1])))/1000)
    peak_time <- c(peak_time, newd$x[which.max(predict.gam(mod_gam,newd))])
  }
  peak_time
  gene_order_vec <- rank(peak_time)
  names(gene_order_vec) <- rownames(gene_expression_mat)
  gene_order_vec_list[[method]] <- gene_order_vec
}
gene_order_vec_list[["GeneTrajectory"]] <- GT_gene_order_vec
gene_order_vec_list[["CellRank"]] <- cellrank_gene_order_vec


gene_order_vec_noCC_list <- list()
for (method in methods){
  cell_pt <- read.csv(sprintf("./data/mouse_dermal_E14_WLS_pairs_noLD/mouse_dermal_E14_WLS_pairs_%s_noCC.csv", method), header = F)
  cells <- which(!is.na(cell_pt[,1]))
  length(cells)
  peak_time <- c()
  for (i in 1:nrow(gene_expression_mat)){
    message(i)
    dat <- data.frame(x=cell_pt[cells,1], y = gene_expression_mat[i,cells])
    mod_gam = gam(y ~ s(x), data = dat)
    newd <- data.frame(x = c(0:(1000*max(cell_pt[cells,1])))/1000)
    peak_time <- c(peak_time, newd$x[which.max(predict.gam(mod_gam,newd))])
  }
  peak_time
  gene_order_vec <- rank(peak_time)
  names(gene_order_vec) <- rownames(gene_expression_mat)
  gene_order_vec_noCC_list[[method]] <- gene_order_vec
}
gene_order_vec_noCC_list[["GeneTrajectory"]] <- GT_gene_order_vec_noCC
gene_order_vec_noCC_list[["CellRank"]] <- cellrank_gene_order_vec_noCC



##### Visualization
all_methods <- names(gene_order_vec_list)

pdf("gene_order_plot.pdf",height = 2, width = 9.5)
for (method in all_methods){
  gene_order_vec <- gene_order_vec_list[[method]]
  genes = genes_to_show[which(genes_to_show %in% names(gene_order_vec))]
  df_plot <- data.frame(x = gene_order_vec[genes],
                        y = .5, 
                        gene = genes,
                        pathways = pathways[genes])
  
  options(ggrepel.max.overlaps = Inf)
  print(
    ggplot(data = df_plot, aes(x=x, y=y, label = gene)) +
      geom_point()+
      ggrepel::geom_label_repel(data = df_plot, aes(color = pathways)) +
      scale_color_manual(values = c("Wnt" = "forestgreen",
                                    "DC" = "blue",
                                    "SHH" = "darkorange2",
                                    "Cdkn1a" = "darkred")) +
      geom_segment(aes(xend=x), yend=-20, linetype=1) +
      theme(axis.line.x = element_line(size=1.5, arrow = grid::arrow(length = unit(0.3, "cm")))) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
      xlab("Gene ordering") + ylab("") + xlim(0,91) + ggtitle(method)
  )
}
dev.off()


pdf("gene_order_plot_noCC.pdf",height = 2, width = 9.5)
for (method in all_methods){
  gene_order_vec <- gene_order_vec_noCC_list[[method]]
  genes = genes_to_show[which(genes_to_show %in% names(gene_order_vec))]
  df_plot <- data.frame(x = gene_order_vec[genes],
                        y = .5, 
                        gene = genes,
                        pathways = pathways[genes])
  
  options(ggrepel.max.overlaps = Inf)
  print(
    ggplot(data = df_plot, aes(x=x, y=y, label = gene)) +
      geom_point()+
      ggrepel::geom_label_repel(data = df_plot, aes(color = pathways)) +
      scale_color_manual(values = c("Wnt" = "forestgreen",
                                    "DC" = "blue",
                                    "SHH" = "darkorange2",
                                    "Cdkn1a" = "darkred")) +
      geom_segment(aes(xend=x), yend=-20, linetype=1) +
      theme(axis.line.x = element_line(size=1.5, arrow = grid::arrow(length = unit(0.3, "cm")))) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
      xlab("Gene ordering") + ylab("") + xlim(0,91) + ggtitle(method)
  )
}
dev.off()

