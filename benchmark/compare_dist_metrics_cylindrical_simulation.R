require(Seurat)
require(lsa)
require(philentropy)
require(GeneTrajectory)
require(ggplot2)

data_S_list <- readRDS("./data/data_S_list_cylinder.rds")

compare_distance_metrics_cylinder <- function(param){
  N_cells = param[1]
  sparsity = param[2]
  i  = param[3]
  dir.path <- sprintf("./data/robustness_evaluation/cylinder_simulation/NC_%s_SPARSITY_%s/data/replicate%s/", N_cells, sparsity, i)
  if (! dir.exists(dir.path)) stop(c(NA,NA,NA))
  setwd(dir.path)
  data_S <- data_S_list[[sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i)]]
  message(sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i))
  expr_N <- apply(data_S@assays$RNA@data > 0, 1, sum)
  
  ground_truth <- as.numeric(gsub(".1$","",rownames(data_S)))
  
  gene.dist.mat.all.list <- list()
  gene.dist.mat.all.list[["GT"]] <- LoadGeneDistMat(paste0(dir.path, "/DIM_10_K_10/"), file_name = "emd.csv")
  gene.dist.mat.all.list[["Euclidean"]] <- as.matrix(dist(data_S@assays$RNA@data, method = "euclidean"))
  gene.dist.mat.all.list[["Manhattan"]] <- as.matrix(dist(data_S@assays$RNA@data, method = "manhattan"))
  gene.dist.mat.all.list[["Pearson"]] <- 1-cor(t(as.matrix(data_S@assays$RNA@data)), method = "pearson")
  gene.dist.mat.all.list[["Spearman"]] <- 1-cor(t(as.matrix(data_S@assays$RNA@data)), method = "spearman")
  gene.dist.mat.all.list[["Cosine"]] <- 1-lsa::cosine(t(as.matrix(data_S@assays$RNA@data)))
  gene.dist.mat.all.list[["JS"]] <- philentropy::JSD(as.matrix(data_S@assays$RNA@data), test.na = FALSE, est.prob = "empirical")
  
  
  cor.all.list <- list()
  for (d in names(gene.dist.mat.all.list)){
    gene.dist.mat.all <- gene.dist.mat.all.list[[d]]
    
    genes <- c(1:500)[which(expr_N[1:500] >= 1)]
    gene.dist.mat <- gene.dist.mat.all[genes, genes]
    dim(gene.dist.mat)
    gene_embedding <- GetGeneEmbedding(gene.dist.mat, K = 5)$diffu.emb
    gene_trajectory <- ExtractGeneTrajectory(gene_embedding, gene.dist.mat, N = 1, t.list = c(1000), K = 5)
    table(gene_trajectory$selected)
    cor_linear <- abs(cor(gene_trajectory$Pseudoorder1, ground_truth[genes], method = "spearman"))
    message(cor_linear)
    
    genes <- c(501:1000)[which(expr_N[501:1000] >= 1)]
    gene.dist.mat <- gene.dist.mat.all[genes, genes]
    dim(gene.dist.mat)
    gene_embedding <- GetGeneEmbedding(gene.dist.mat, K = 5)$diffu.emb
    return_angle <- function(x){
      if (x[2]>=0) return(acos(x[1]/sqrt(x[1]^2+x[2]^2)))
      if (x[2]<0) return(acos(-x[1]/sqrt(x[1]^2+x[2]^2))+pi)
    }
    angle_vec <- apply(gene_embedding[,1:2], 1, return_angle)
    cor_vec <- c(abs(cor(x = angle_vec, y = ground_truth[genes], method = "spearman")))
    for (idx in 2:length(genes)){
      cor_vec <- c(cor_vec, abs(cor(x = angle_vec[c(idx:length(genes),1:(idx-1))], y = ground_truth[genes], method = "spearman")))
    }
    idx = which.max(cor_vec)
    cor_circular <- max(cor_vec)
    message(cor_circular)
    
    cor_all = c(cor_linear, cor_circular)
    
    cor.all.list[[d]] <- cor_all
  }
  
  cor.all.df <- as.data.frame(do.call(rbind, cor.all.list))
  colnames(cor.all.df) <- c("linear", "circular")
  cor.all.df <- cbind(metric = names(gene.dist.mat.all.list),
                      cor.all.df)
  cor.all.df <- cbind(sample = sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i),
                      cor.all.df)
  cor.all.df
}


param_list0 <- list()
for (N_cells in c(1000)){
  for (sparsity in c(0.025, 0.05, 0.1)){
    for (i in 1:10){
      param_list0[[sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i)]] <- c(N_cells, sparsity, i)
    }
  }
}

distance_metrics_comparison_output <- lapply(param_list0, compare_distance_metrics_cylinder)
saveRDS(distance_metrics_comparison_output, "distance_metrics_comparison_output_cylinder.rds")

distance_metrics_comparison_output <- distance_metrics_comparison_output_cylinder
summary_df <- do.call(rbind, distance_metrics_comparison_output)
summary_df$metric[which(summary_df$metric == "GT")] <- "Wasserstein"
summary_df$sparsity <- unlist(strsplit(summary_df$sample, "_"))[c(1:nrow(summary_df))*6-2]
head(summary_df)

ggplot(data = summary_df, aes(x = metric, y = linear, fill = metric)) + 
  geom_boxplot() + 
  facet_grid(~ sparsity, margins = FALSE) & RotatedAxis() & ylim(c(0,1)) & ylab("Spearman correlation")
ggplot(data = summary_df, aes(x = metric, y = circular, fill = metric)) + 
  geom_boxplot() + 
  facet_grid(~ sparsity, margins = FALSE) & RotatedAxis() & ylim(c(0.5,1)) & ylab("Spearman correlation")

