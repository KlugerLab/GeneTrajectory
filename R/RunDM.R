#' Run Diffusion Map on a Seurat object
#'
#' @param object Seurat object
#' @param reduction string; Dimensionality reduction to use, default: 'pca'
#' @param dims integer; Which dimensions to use as input features
#' @param K integer; Adaptive kernel bandwidth for each point set to be the distance to its `K`-th nearest neighbor.
#' @param sigma numeric; Fixed kernel bandwidth, `sigma` will be ignored if `K` is specified.
#' @param n.components integer; Number of leading nontrivial eigenvectors to export
#' @param t integer; Number of diffusion times
#' @param dist.mat matrix; Precomputed distance matrix (optional)
#' @param reduction.key string; Dimensional reduction key, specifies the string before the number for the dimension names. 'DM_' by default
#'
#' @return
#' Returns a Seurat object containing a Diffusion Map representation
#'
#'
#' @export
#'
#' @examples
#'
#'

RunDM <- function(object,
                  reduction = "pca",
                  dims = 1:30,
                  K = 10,
                  sigma = NULL,
                  n.components = 30,
                  t = 1,
                  dist.mat = NULL,
                  reduction.key = "DM_"
                  ){
  data.S <- object
  if (is.null(dist.mat)) dist.mat <- as.matrix(dist(data.S@reductions[[reduction]]@cell.embeddings[, dims]))
  dm.emb <- diffusion.map(dist.mat,
                          sigma = sigma,
                          K = K,
                          nEV = n.components+1,
                          t = t
                          )[["diffu.emb"]][, 2:(n.components+1)] #remove the first eigenvector (trivial)
  colnames(dm.emb) <- paste0(reduction.key, 1:ncol(dm.emb))

  data.S[["dm"]] <- CreateDimReducObject(embeddings = dm.emb, key = reduction.key, assay = DefaultAssay(data.S))

  data.S
}
