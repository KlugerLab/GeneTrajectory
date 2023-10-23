
compare_GT_Monocle3_cyclic <- function(param){
  N_cells <- param[1]
  sparsity <- param[2]
  i  <- param[3]
  message(sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i))
  data_S <- data_S_list[[sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i)]]
  ground_truth <- as.numeric(gsub(".1$","",rownames(data_S)))
  dir.path <- sprintf("./data/benchmark/cylinder_simulation/NC_%s_SPARSITY_%s/Monocle3/replicate%s", N_cells, sparsity, i)
  setwd(dir.path)
  
  N_g = 100
  set.seed(1024)
  gene_ind1 <- sample(c(1:500)[which(apply(data_S@assays$RNA@data[1:500,], 1, sum) > 0)], N_g)
  set.seed(4201)
  gene_ind2 <- sort(sample(c(501:1000)[which(apply(data_S@assays$RNA@data[501:1000,], 1, sum) > 0)], N_g))
  
  cell_pt <- read.csv("CT_cyclic.csv", header = F)[,1]
  
  tic()
  peak_time <- c()
  for (i in gene_ind2){
    dat <- data.frame(x=cell_pt, y = data_S[["RNA"]]@data[i,])
    mod_gam = gam(y ~ s(x), data = dat)
    
    newd <- data.frame(x = c(1:(100*max(cell_pt)))/100)
    peak_time <- c(peak_time, newd$x[which.max(predict.gam(mod_gam,newd))])
  }
  peak_time
  toc()
  
  
  cor_vec <- c(abs(cor(x = peak_time, y = ground_truth[gene_ind2], method = "spearman")))
  for (idx in 2:length(gene_ind2)){
    cor_vec <- c(cor_vec, abs(cor(x = peak_time[c(idx:length(gene_ind2),1:(idx-1))], y = ground_truth[gene_ind2], method = "spearman")))
  }
  idx = which.max(cor_vec)
  cor_cyclic <- max(cor_vec)
  message(cor_cyclic)
  
  cor_cyclic
}


compare_GT_Monocle2_cyclic <- function(param){
  N_cells <- param[1]
  sparsity <- param[2]
  i  <- param[3]
  message(sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i))
  data_S <- data_S_list[[sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i)]]
  ground_truth <- as.numeric(gsub(".1$","",rownames(data_S)))
  dir.path <- sprintf("./data/benchmark/cylinder_simulation/NC_%s_SPARSITY_%s/Monocle2/replicate%s", N_cells, sparsity, i)
  setwd(dir.path)
  
  N_g = 100
  set.seed(1024)
  gene_ind1 <- sample(c(1:500)[which(apply(data_S@assays$RNA@data[1:500,], 1, sum) > 0)], N_g)
  set.seed(4201)
  gene_ind2 <- sort(sample(c(501:1000)[which(apply(data_S@assays$RNA@data[501:1000,], 1, sum) > 0)], N_g))
  
  cell_pt <- read.csv("CT_cyclic.csv", header = F)[,1]
  
  tic()
  peak_time <- c()
  for (i in gene_ind2){
    dat <- data.frame(x=cell_pt, y = data_S[["RNA"]]@data[i,])
    mod_gam = gam(y ~ s(x), data = dat)
    
    newd <- data.frame(x = c(1:(100*max(cell_pt)))/100)
    peak_time <- c(peak_time, newd$x[which.max(predict.gam(mod_gam,newd))])
  }
  peak_time
  toc()
  
  
  cor_vec <- c(abs(cor(x = peak_time, y = ground_truth[gene_ind2], method = "spearman")))
  for (idx in 2:length(gene_ind2)){
    cor_vec <- c(cor_vec, abs(cor(x = peak_time[c(idx:length(gene_ind2),1:(idx-1))], y = ground_truth[gene_ind2], method = "spearman")))
  }
  idx = which.max(cor_vec)
  cor_cyclic <- max(cor_vec)
  message(cor_cyclic)
  
  cor_cyclic
}



compare_GT_PAGA_cyclic <- function(param){
  N_cells <- param[1]
  sparsity <- param[2]
  i  <- param[3]
  message(sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i))
  data_S <- data_S_list[[sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i)]]
  ground_truth <- as.numeric(gsub(".1$","",rownames(data_S)))
  dir.path <- sprintf("./data/benchmark/cylinder_simulation/NC_%s_SPARSITY_%s/PAGA_DPT/replicate%s", N_cells, sparsity, i)
  setwd(dir.path)
  
  N_g = 100
  set.seed(1024)
  gene_ind1 <- sample(c(1:500)[which(apply(data_S@assays$RNA@data[1:500,], 1, sum) > 0)], N_g)
  set.seed(4201)
  gene_ind2 <- sort(sample(c(501:1000)[which(apply(data_S@assays$RNA@data[501:1000,], 1, sum) > 0)], N_g))
  
  cell_pt <- read.csv("CT_cyclic.csv", header = F)[,1]
  
  tic()
  peak_time <- c()
  for (i in gene_ind2){
    dat <- data.frame(x=cell_pt, y = data_S[["RNA"]]@data[i,])
    mod_gam = gam(y ~ s(x), data = dat)
    
    newd <- data.frame(x = c(1:(100*max(cell_pt)))/100)
    peak_time <- c(peak_time, newd$x[which.max(predict.gam(mod_gam,newd))])
  }
  peak_time
  toc()
  
  cor_vec <- c(abs(cor(x = peak_time, y = ground_truth[gene_ind2], method = "spearman")))
  for (idx in 2:length(gene_ind2)){
    cor_vec <- c(cor_vec, abs(cor(x = peak_time[c(idx:length(gene_ind2),1:(idx-1))], y = ground_truth[gene_ind2], method = "spearman")))
  }
  idx = which.max(cor_vec)
  cor_cyclic <- max(cor_vec)
  message(cor_cyclic)
  
  cor_cyclic
}


compare_GT_GeneTrajectory_cyclic <- function(param){
  N_cells <- param[1]
  sparsity <- param[2]
  i  <- param[3]
  message(sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i))
  data_S <- data_S_list[[sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i)]]
  ground_truth <- as.numeric(gsub(".1$","",rownames(data_S)))
  dir.path <- sprintf("./data/benchmark/cylinder_simulation/NC_%s_SPARSITY_%s/GeneTrajectory/replicate%s/DIM_5_K_10/", N_cells, sparsity, i)
  
  N_g = 100
  
  gene.dist.mat <- LoadGeneDistMat(dir.path, file_name = "emd.csv")
  genes <- which(rownames(data_S)[501:1000] %in% rownames(gene.dist.mat))+500
  set.seed(4201)
  gene_ind2 <- sort(sample(genes, N_g))
  gene.dist.mat <- gene.dist.mat[gene_ind2, gene_ind2]
  dim(gene.dist.mat)
  gene_embedding <- GetGeneEmbedding(gene.dist.mat, K = 5)$diffu.emb
  
  return_angle <- function(x){
    if (x[2]>=0) return(acos(x[1]/sqrt(x[1]^2+x[2]^2)))
    if (x[2]<0) return(acos(-x[1]/sqrt(x[1]^2+x[2]^2))+pi)
  }
  
  
  angle_vec <- apply(gene_embedding[,1:2], 1, return_angle)
  
  cor_vec <- c(abs(cor(x = angle_vec, y = ground_truth[gene_ind2], method = "spearman")))
  for (idx in 2:length(gene_ind2)){
    cor_vec <- c(cor_vec, abs(cor(x = angle_vec[c(idx:length(gene_ind2),1:(idx-1))], y = ground_truth[gene_ind2], method = "spearman")))
  }
  idx = which.max(cor_vec)
  cor_cyclic <- max(cor_vec)
  message(cor_cyclic)
  
  cor_cyclic
}



runSlingShot_cyclic <- function(param){
  N_cells <- param[1]
  sparsity <- param[2]
  i  <- param[3]
  message(sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i))
  data_S <- data_S_list[[sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i)]]
  sce <- as.SingleCellExperiment(data_S)
  reducedDims(sce)
  sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP')
  
  DimPlot(data_S, reduction = "umap", label=T, group.by = "seurat_clusters", repel=T)
  
  start.clus = sce$seurat_clusters[which.min(sce$cell_pseudotime_cyclic)]
  
  lin1 <- getLineages(reducedDims(sce)$UMAP, sce$seurat_clusters, start.clus = start.clus)#, end.clus = '0')
  crv1 <- getCurves(lin1)
  
  spt <- as.data.frame(slingPseudotime(crv1))
  head(spt)
  
  if (ncol(spt) == 1){
    CT <- spt[,1]
  }
  if (ncol(spt) > 1){
    CT <- spt[,which.max(unlist(lapply(1:ncol(spt), function(x){
      sum(!is.na(spt[,x]))
    })))]
  }
  dir.path <- sprintf("./data/benchmark/cylinder_simulation/NC_%s_SPARSITY_%s/SlingShot/replicate%s", N_cells, sparsity, i)
  dir.create(dir.path, recursive = TRUE)
  write.table(CT, paste0(dir.path, "/CT_cyclic.csv"), row.names = F, col.names = F, sep = ",")
}




