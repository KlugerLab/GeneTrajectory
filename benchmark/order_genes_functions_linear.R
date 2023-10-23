
compare_GT_Monocle3_linear <- function(param){
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
  gene_ind1 <- sample(c(1:500)[which(apply(data_S@assays$RNA@data[1:500,] > 0, 1, sum) > 0)], N_g)
  set.seed(4201)
  gene_ind2 <- sort(sample(c(501:1000)[which(apply(data_S@assays$RNA@data[501:1000,] > 0, 1, sum) > 0)], N_g))
  
  cell_pt <- read.csv("CT_linear.csv", header = F)[,1]
  
  tic()
  peak_time <- c()
  for (i in gene_ind1){
    dat <- data.frame(x=cell_pt, y = data_S[["RNA"]]@data[i,])
    mod_gam = gam(y ~ s(x), data = dat)
    
    newd <- data.frame(x = c(1:(100*max(cell_pt)))/100)
    peak_time <- c(peak_time, newd$x[which.max(predict.gam(mod_gam,newd))])
  }
  peak_time
  toc()
  
  cor_linear <- abs(cor(peak_time, ground_truth[gene_ind1], method = "spearman"))
  
  cor_linear
}


compare_GT_Monocle2_linear <- function(param){
  N_cells = 500
  sparsity = 0.025
  i = 1
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
  gene_ind1 <- sample(c(1:500)[which(apply(data_S@assays$RNA@data[1:500,] > 0, 1, sum) > 0)], N_g)
  set.seed(4201)
  gene_ind2 <- sort(sample(c(501:1000)[which(apply(data_S@assays$RNA@data[501:1000,] > 0, 1, sum) > 0)], N_g))
  
  cell_pt <- read.csv("CT_linear.csv", header = F)[,1]
  
  tic()
  peak_time <- c()
  for (i in gene_ind1){
    dat <- data.frame(x=cell_pt, y = data_S[["RNA"]]@data[i,])
    mod_gam = gam(y ~ s(x), data = dat)
    
    newd <- data.frame(x = c(1:(100*max(cell_pt)))/100)
    peak_time <- c(peak_time, newd$x[which.max(predict.gam(mod_gam,newd))])
  }
  peak_time
  toc()
  
  cor_linear <- abs(cor(peak_time, ground_truth[gene_ind1], method = "spearman"))
  message(cor_linear)
  
  cor_linear
}


compare_GT_PAGA_linear <- function(param){
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
  gene_ind1 <- sample(c(1:500)[which(apply(data_S@assays$RNA@data[1:500,] > 0, 1, sum) > 0)], N_g)
  set.seed(4201)
  gene_ind2 <- sort(sample(c(501:1000)[which(apply(data_S@assays$RNA@data[501:1000,] > 0, 1, sum) > 0)], N_g))
  
  cell_pt <- read.csv("CT_linear.csv", header = F)[,1]
  
  tic()
  peak_time <- c()
  for (i in gene_ind1){
    dat <- data.frame(x=cell_pt, y = data_S[["RNA"]]@data[i,])
    mod_gam = gam(y ~ s(x), data = dat)
    
    newd <- data.frame(x = c(1:(100*max(cell_pt)))/100)
    peak_time <- c(peak_time, newd$x[which.max(predict.gam(mod_gam,newd))])
  }
  peak_time
  toc()
  
  cor_linear <- abs(cor(peak_time, ground_truth[gene_ind1], method = "spearman"))
  
  message(cor_linear)
  
  cor_linear
}


compare_GT_GeneTrajectory_linear <- function(param){
  N_cells <- param[1]
  sparsity <- param[2]
  i  <- param[3]
  message(sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i))
  data_S <- data_S_list[[sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i)]]
  ground_truth <- as.numeric(gsub(".1$","",rownames(data_S)))
  dir.path <- sprintf("./data/benchmark/cylinder_simulation/NC_%s_SPARSITY_%s/GeneTrajectory/replicate%s/DIM_5_K_10/", N_cells, sparsity, i)
  
  N_g = 100
  if (! file.exists(paste0(dir.path, "emd.csv"))) dir.path <- sprintf("./data/benchmark/cylinder_simulation/NC_%s_SPARSITY_%s/GeneTrajectory/replicate%s/DIM_5_K_15/", N_cells, sparsity, i)
  
  gene.dist.mat <- LoadGeneDistMat(dir.path, file_name = "emd.csv")
  genes <- which(rownames(data_S)[1:500] %in% rownames(gene.dist.mat))
  set.seed(1024)
  gene_ind1 <- sample(genes, N_g)
  gene.dist.mat <- gene.dist.mat[gene_ind1, gene_ind1]
  dim(gene.dist.mat)
  gene_embedding <- GetGeneEmbedding(gene.dist.mat, K = 5)$diffu.emb
  gene_trajectory <- ExtractGeneTrajectory(gene_embedding, gene.dist.mat, N = 1, t.list = c(1000), K = 5)
  table(gene_trajectory$selected)
  cor_linear <- abs(cor(gene_trajectory$Pseudoorder1, ground_truth[gene_ind1], method = "spearman"))
  
  cor_linear
}



runSlingShot_linear <- function(param){
  N_cells <- param[1]
  sparsity <- param[2]
  i  <- param[3]
  message(sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i))
  data_S <- data_S_list[[sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i)]]
  sce <- as.SingleCellExperiment(data_S)
  reducedDims(sce)
  sce <- slingshot(sce, clusterLabels = 'seurat_clusters', reducedDim = 'UMAP')
  
  DimPlot(data_S, reduction = "umap", label=T, group.by = "seurat_clusters", repel=T)
  
  start.clus = sce$seurat_clusters[which.min(sce$cell_pseudotime_linear)]
  
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
  write.table(CT, paste0(dir.path, "/CT_linear.csv"), row.names = F, col.names = F, sep = ",")
}


compare_GT_SlingShot_linear <- function(param){
  N_cells <- param[1]
  sparsity <- param[2]
  i  <- param[3]
  message(sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i))
  data_S <- data_S_list[[sprintf("NC_%s_SPARSITY_%s_replicate_%s", N_cells, sparsity, i)]]
  ground_truth <- as.numeric(gsub(".1$","",rownames(data_S)))
  dir.path <- sprintf("./data/benchmark/cylinder_simulation/NC_%s_SPARSITY_%s/SlingShot/replicate%s", N_cells, sparsity, i)
  setwd(dir.path)
  
  cell_pt <- read.csv("CT_linear.csv", header = F)[,1]
  data_S <- data_S[,which(!is.na(cell_pt))]
  gene_expr_N <- apply(as.matrix(data_S[["RNA"]]@data) > 0, 1, sum)#/ncol(data_S)
  genes <- which(gene_expr_N[1:500] > 0)
  
  N_g = 100
  set.seed(1024)
  gene_ind1 <- sample(genes, N_g)
  
  tic()
  peak_time <- c()
  for (i in gene_ind1){
    dat <- data.frame(x=cell_pt, y = data_S[["RNA"]]@data[i,])
    mod_gam = gam(y ~ s(x), data = dat)
    
    newd <- data.frame(x = c(1:(100*max(cell_pt)))/100)
    peak_time <- c(peak_time, newd$x[which.max(predict.gam(mod_gam,newd))])
  }
  peak_time
  toc()
  
  cor_linear <- abs(cor(peak_time, ground_truth[gene_ind1], method = "spearman"))
  
  cor_linear
}




