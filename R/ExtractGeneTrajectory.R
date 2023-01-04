

ExtractGeneTrajectory <- function(gene.embedding,
                                  dist.mat,
                                  dims = 1:5,
                                  N,
                                  t.list,
                                  K = 10,
                                  quantile = 0.02
                                  ){
  dist.to.origin <- sqrt(apply(gene.embedding[, dims]^2, 1, sum))
  summary.df <- as.data.frame(gene.embedding[, dims])
  summary.df$selected <- "Other"
  nBranch <- N
  genes <- rownames(gene.embedding)
  diffusion.mat <- get.RW.matrix(dist.mat, K = K)

  for (i in 1:nBranch){
    if (length(which(summary.df$selected == "Other")) == 0) stop("Wrong!")

    if (length(which(summary.df$selected != "Other")) != 0) dist.to.origin[which(summary.df$selected != "Other")] <- -Inf

    message(genes[which.max(dist.to.origin)])

    seed <- rep(0, nrow(gene.embedding))
    seed[which.max(dist.to.origin)] <- 1

    seed.diffused <- seed
    max.T <- t.list[i]
    for (ii in 1:max.T){
      seed.diffused <- diffusion.mat %*% matrix(seed.diffused, ncol = 1)
    }

    cutoff <- max(seed.diffused) * quantile

    summary.df$selected[which(seed.diffused > cutoff & summary.df$selected == "Other")]  <- paste0("Trajectory-", i)
    summary.df$seed.diffused <- -log(seed.diffused)


    summary.df <- cbind(summary.df,
                       construct.pseudoorder(dist.mat,
                                             which(summary.df$selected  == paste0("Trajectory-", i)),
                                             max.id = genes[which.max(dist.to.origin)]))

    colnames(summary.df)[ncol(summary.df)] <- paste0("Pseudoorder", i)
  }

  summary.df
}
