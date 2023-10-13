#' Run Diffusion Map on a precomputed distance matrix
#'
#' @param dist.mat matrix; Precomputed distance matrix (symmetric)
#' @param K integer; Adaptive kernel bandwidth for each point set to be the distance to its `K`-th nearest neighbor.
#' @param sigma numeric; Fixed kernel bandwidth, `sigma` will be ignored if `K` is specified.
#' @param nEV integer; Number of leading eigenvectors to export
#' @param t integer; Number of diffusion times
#' 
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

diffusion.map <- function(dist.mat, 
                          K = 10, 
                          sigma = NULL,
                          nEV = 30,
                          t = 1
                          ){
  K <- min(K, nrow(dist.mat))
  dists <- as.matrix(dist.mat)
  sigma.list <- c()
  if (is.null(sigma)){
    for (i in 1:nrow(dists)){
      sigma.list[i] <- sort(dists[i,])[K]
    }
  }
  if (!is.null(sigma)){
    sigma.list <- rep(sigma, nrow(dists))
  }
  
  affinity.matrix <- matrix(0, nrow=nrow(dists), ncol=ncol(dists))
  
  for (i in 1:nrow(affinity.matrix)){
    dist.vec <- dists[i, ]
    dist.vec[is.na(dist.vec)] <- 10^6 #default value to fit in the missing values in the distance matrix
    affinity.matrix[i, ] <- exp(-dist.vec^2/(sigma.list[i]^2))
  }
  
  affinity.matrix.2 <- (affinity.matrix + t(affinity.matrix))/2
  normalized.vec <- sqrt(1/apply(affinity.matrix.2, 1, sum))
  affinity.matrix.3 <- t(affinity.matrix.2 * normalized.vec) * normalized.vec
  
  N.EV <- min(nrow(affinity.matrix.3), nEV)
  E.list <- rARPACK::eigs_sym(affinity.matrix.3, k = N.EV)
  
  diffu.emb <- matrix(0, nrow = nrow(E.list$vector), ncol = N.EV)
  for(i in 1:N.EV){
    diffu.emb[,i] <- E.list$vector[,i]*normalized.vec*(E.list$values[i]^t)
  }
  colnames(diffu.emb) <- paste0("EV.", 1:ncol(diffu.emb))
  rownames(diffu.emb) <- colnames(dist.mat)
  
  res <- list()
  res[["diffu.emb"]] <- diffu.emb
  res[["eigen.vals"]] <- E.list$values[1:N.EV]
  
  class(res) <- c("diffusion.map", "list")
  return(res)
}