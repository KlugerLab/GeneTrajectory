AddGeneBinScore <- function(object, 
                               gene_trajectory, 
                               N_bin = 5, 
                               trajectories = 1:2,
                               assay = "RNA",
                               prefix = "Trajectory",
                               use.module_score = F,
                               binarize = F,
                               reverse = NULL,
                               echo = F){
  data_S <- object
  prefix0 <- prefix
  for (trajectory in trajectories){
    prefix <- paste0(prefix0, trajectory, "_genes")
    gene_trajectory_reordered <- gene_trajectory[order(gene_trajectory[, paste0("Pseudoorder", trajectory)]),]
    genes <- rownames(gene_trajectory_reordered)[which(gene_trajectory_reordered[, paste0("Pseudoorder", trajectory)] > 0)]
    if (! is.null(reverse)){
      if (isTRUE(reverse[trajectory])) genes <- rev(genes)
    }
    length(genes)
    
    step <- length(genes)/N_bin
    metadata <- data.frame(rep(0, ncol(data_S)))
    for (i in 1:N_bin){
      start <- ceiling((i-1)*step+1)
      end <- min(ceiling(i*step), length(genes))
      #print(start:end)
      genes_subset <- genes[start:end]
      if (echo) message(paste(genes_subset, collapse = ", "))
      if (!use.module_score){
        normalized_gc <- data_S[[assay]]@data[genes_subset,]#/apply(data_S[[assay]]@data[genes_subset,],1,sum)
        #print(apply(normalized_gc,1,sum))
        if (binarize) normalized_gc <- (normalized_gc > 0)
        metadata[,i] <- apply(normalized_gc,2,sum)
        metadata[,i] <- metadata[,i]/length(genes_subset)#sum(metadata[,i])
      }
      if (use.module_score){
        feature_list <- list()
        feature_list[["tmp"]] <- genes_subset
        data_S <- AddModuleScore(data_S, features = feature_list, assay = assay, name = "tmp")
        metadata[,i] <- data_S$tmp1
      }
    }
    rownames(metadata) <- colnames(data_S)
    colnames(metadata) <- paste0(prefix, 1:N_bin)
    #print(head(metadata))
    data_S <- AddMetaData(data_S, metadata)
  }
  
  data_S
}