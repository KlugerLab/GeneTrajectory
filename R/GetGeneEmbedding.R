#' Run Diffusion Map on the gene-gene Wasserstein distance matrix
#'
#' @param dist.mat matrix; Gene-gene Wasserstein distance matrix (symmetric)
#' @param K integer; Adaptive kernel bandwidth for each point set to be the distance to its `K`-th nearest neighbor.
#' @param sigma numeric; Fixed kernel bandwidth, `sigma` will be ignored if `K` is specified.
#' @param nEV integer; Number of leading eigenvectors to export
#' @param t integer; Number of diffusion times
#'
#' @return 
#' List with the following elements:
#' \tabular{ll}{
#'    \code{diffu.emb} \tab Diffusion Map embedding with the leading `K` eigenvectors \cr
#'    \tab \cr
#'    \code{eigen.vals} \tab Eigenvalues associated with corresponding eigenvectors \cr
#' }
#' 
#' 
#' @export
#'
#' @examples
#' 
#' 

GetGeneEmbedding <- function(dist.mat, 
                               K = 10, 
                               sigma = NULL, 
                               nEV = 30,
                               t = 1){
  dm.res <- diffusion.map(dist.mat, 
                          K = K, 
                          sigma = sigma, 
                          nEV = nEV+1,
                          t = t)
  dm.res[["diffu.emb"]] <- dm.res[["diffu.emb"]][, 2:(nEV+1)]
  colnames(dm.res[["diffu.emb"]]) <- paste0("DM_", 1:ncol(dm.res[["diffu.emb"]]))
  dm.res[["eigen.vals"]] <- dm.res[["eigen.vals"]][2:(nEV+1)]
  dm.res
}