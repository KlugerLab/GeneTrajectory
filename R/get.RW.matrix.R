#' Convert a distance matrix into a random-walk matrix based on adaptive Gaussian kernel
#'
#' @param dist.mat matrix; Precomputed distance matrix (symmetric)
#' @param K integer; Adaptive kernel bandwidth for each point set to be the distance to its `K`-th nearest neighbor.
#'
#'
#' @return
#' Random-walk matrix
#'
#'
#' @export
#'
#' @examples
#'
#'

get.RW.matrix <- function(dist.mat,
                          K = 10){
  dists <- as.matrix(dist.mat)
  sigma.list <- c()
  for (i in 1:nrow(dists)){
    sigma.list[i] <- sort(dists[i,])[K]
  }

  affinity.matrix <- matrix(0, nrow=nrow(dists), ncol=ncol(dists))
  for (i in 1:nrow(affinity.matrix)){
    dist.vec <- dists[i, ]
    dist.vec[is.na(dist.vec)] <- 10^6
    affinity.matrix[i, ] <- exp(-dist.vec^2/(sigma.list[i]^2))
  }
  affinity.matrix.2 <- (affinity.matrix + t(affinity.matrix))/2
  normalized.vec <- 1/apply(affinity.matrix.2, 1, sum)
  diffusion.mat <- affinity.matrix.2 * normalized.vec
  diffusion.mat
}
