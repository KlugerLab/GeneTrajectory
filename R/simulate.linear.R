#' Simulate linear gene process
#'
#' @param N_cells integer; Number of cells
#' @param N_genes integer; Number of genes
#' @param model character; Count model ("poisson" or "negbin")
#' @param meanlog numeric; Mean of log normal distribution
#' @param sdlog numeric; Standard deviation of log normal distribution
#' @param scale numeric; Scale of UMI counts
#' @param seed integer; Random seed
#' @param maxT numeric; Maximum cell pseudotime
#' @param sort boolean; Whether to sort genes based on their peak times
#' @param sparsity numeric; Sparsity of count matrix
#' @param theta numeric; Dipersion parameter for negative binomial model
#'
#' @return
#' Returns a gene-by-cell count matrix
#'
#' @export simulate.linear
#'
#' @examples
#'
#'

simulate.linear <- function(N_cells = 1000,
                            N_genes = 200,
                            model = "poisson",
                            meanlog = 0,
                            sdlog = 0.25,
                            scale = 25,
                            seed = 1,
                            maxT = 15,
                            sort = TRUE,
                            sparsity = NULL,
                            theta = 10){
  if (!is.null(sparsity)) {
    message(sprintf("Simulating the linear process, %s cells, %s genes, sparsity: %s, random.seed: %s, using %s model", N_cells, N_genes, sparsity, seed, model))
  }
  if (is.null(sparsity)) {
    message(sprintf("Simulating the linear process, %s cells, %s genes, no sparsification, random.seed: %s, using %s model", N_cells, N_genes, seed, model))
  }

  cell_pt_list <- list()
  gene_pt_list <- list()
  Np <- length(N_cells)
  for (i in 1:Np){
    set.seed(seed+10*i+1)
    cell_pt_list[[i]] <- sort(runif(n = N_cells[i], min = 0, max = maxT))
    set.seed(seed+10*i+2)
    gene_pt_list[[i]] <- sort(runif(n = N_genes[i], min = 0, max = maxT))
  }

  cell_pt <- cell_pt_list[[1]]
  gene_pt <- gene_pt_list[[1]]

  gc_mat <- matrix(0, nrow = sum(N_genes), ncol = sum(N_cells))

  N <- sum(N_genes)
  set.seed(seed)
  sigma <- 2*rlnorm(n = N, meanlog = meanlog, sdlog = sdlog)
  set.seed(seed+1)
  alpha <- scale * rlnorm(n = N, meanlog = meanlog, sdlog = sdlog)

  for (i in 1:N){
    mean <- gene_pt[i]
    sd <- sigma[i]
    scale <- alpha[i]
    dist <- abs(cell_pt - mean)
    expectation <- floor(scale * exp(-dist^2/(2*sd^2)))

    set.seed(seed + i)
    if (model == "poisson") gc_mat[i,] <- rpois(n = length(expectation), lambda = expectation)
    set.seed(seed + i)
    if (model == "negbin") gc_mat[i,] <- rnegbin(n = length(expectation), mu = expectation, theta = theta)
  }

  if (!is.null(sparsity)) {
    flattened_mat <- as.numeric(gc_mat)
    sparsity <- min(sparsity, sum(flattened_mat>0)/length(flattened_mat))
    set.seed(seed)
    idx <- sample(1:length(flattened_mat), size = floor(sparsity*length(flattened_mat)), prob = flattened_mat/sum(flattened_mat))
    binarized_vec <- rep(0, length(flattened_mat))
    binarized_vec[idx] <- 1
    gc_mat <- matrix(binarized_vec * flattened_mat, nrow = nrow(gc_mat))
  }

  rownames(gc_mat) <- unlist(gene_pt_list)
  colnames(gc_mat) <- unlist(cell_pt_list)
  gc_mat
}


