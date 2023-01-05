#' Add gene bin score
#'
#' @param object Seurat object
#' @param gene.trajectory Gene trajectory data frame
#' @param N.bin integer; How many gene bins
#' @param trajectories vector; Which gene trajectories to define gene bin score
#' @param assay string; Which assay to use
#' @param prefix string; String added to the names in the metadata
#' @param reverse vector; Whether to reverse the order of genes along each trajectory
#' @param echo boolean; Whether to print messages
#'
#'
#' @return
#' Seurat object with gene bin score added in the metadata
#'
#' @export
#'
#' @examples
#'
#'
AddGeneBinScore <- function(object,
                            gene.trajectory,
                            N.bin = 5,
                            trajectories = 1:2,
                            assay = "RNA",
                            prefix = "Trajectory",
                            reverse = NULL,
                            echo = F){

  data.S <- object
  prefix0 <- prefix
  for (trajectory in trajectories){
    prefix <- paste0(prefix0, trajectory, "_genes")
    gene.trajectory.reordered <- gene.trajectory[order(gene.trajectory[, paste0("Pseudoorder", trajectory)]),]
    genes <- rownames(gene.trajectory.reordered)[which(gene.trajectory.reordered[, paste0("Pseudoorder", trajectory)] > 0)]
    if (! is.null(reverse)){
      if (isTRUE(reverse[trajectory])) genes <- rev(genes)
    }
    length(genes)

    step <- length(genes)/N.bin
    metadata <- data.frame(rep(0, ncol(data.S)))
    for (i in 1:N.bin){
      start <- ceiling((i-1)*step+1)
      end <- min(ceiling(i*step), length(genes))

      genes.subset <- genes[start:end]
      if (echo) message(paste(genes.subset, collapse = ", "))
      normalized.gc <- data.S[[assay]]@data[genes.subset,]
      normalized.gc <- (normalized.gc > 0)
      metadata[,i] <- apply(normalized.gc,2,sum)
      metadata[,i] <- metadata[,i]/length(genes.subset)
    }

    rownames(metadata) <- colnames(data.S)
    colnames(metadata) <- paste0(prefix, 1:N.bin)
    data.S <- AddMetaData(data.S, metadata)
  }

  data.S
}
