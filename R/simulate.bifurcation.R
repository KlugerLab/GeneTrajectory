#' Simulate bifurcating gene process
#'
#' @param N_cells integer vector; Number of cells in the parent process and two daughter processes
#' @param N_genes integer vector; Number of genes associated with the parent process and two daughter processes
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
#' @export simulate.bifurcation
#'
#' @examples
#'

simulate.bifurcation <- function(N_cells = 5*c(100, 50, 50),
                                 N_genes = 5*c(100, 50, 50),
                                 model = "poisson",
                                 meanlog = 0,
                                 sdlog = 0.25,
                                 scale = 25,
                                 seed = 1,
                                 maxT = 15,
                                 sort = TRUE,
                                 sparsity = 0.1,
                                 theta = 10){

  if (!is.null(sparsity)) {
    message(sprintf("Simulating the bifurcation process, %s cells, %s genes, sparsity: %s, random.seed: %s, using %s model", sum(N_cells), sum(N_genes), sparsity, seed, model))
  }
  if (is.null(sparsity)) {
    message(sprintf("Simulating the bifurcation process, %s cells, %s genes, no sparsification, random.seed: %s, using %s model", sum(N_cells), sum(N_genes), seed, model))
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

  cell_pt_mat <- matrix(0, nrow = Np, ncol = sum(N_cells))
  gene_pt_mat <- matrix(0, nrow = Np, ncol = sum(N_genes))

  cell_pt_mat[1,] <- maxT
  cell_pt_mat[1,1:N_cells[1]] <- cell_pt_list[[1]]
  gene_pt_mat[1,] <- maxT
  gene_pt_mat[1,1:N_genes[1]] <- gene_pt_list[[1]]

  N_cells_summed <- unlist(lapply(1:Np, FUN = function(i){
    sum(N_cells[1:i])
  }))
  N_genes_summed <- unlist(lapply(1:Np, FUN = function(i){
    sum(N_genes[1:i])
  }))

  for (i in 2:Np){
    cell_pt_mat[i,(N_cells_summed[i-1]+1):N_cells_summed[i]] <- cell_pt_list[[i]]
    gene_pt_mat[i,(N_genes_summed[i-1]+1):N_genes_summed[i]] <- gene_pt_list[[i]]
  }

  gc_mat <- matrix(0, nrow = sum(N_genes), ncol = sum(N_cells))

  N <- sum(N_genes)
  set.seed(seed)
  sigma <- 2*rlnorm(n = N, meanlog = meanlog, sdlog = sdlog)
  set.seed(seed+1)
  alpha <- scale * rlnorm(n = N, meanlog = meanlog, sdlog = sdlog)

  for (i in 1:N){
    mean <- gene_pt_mat[, i]
    sd <- sigma[i]
    scale <- alpha[i]

    dist <- unlist(lapply(1:ncol(cell_pt_mat), function(i){
      sum(abs(cell_pt_mat[, i] - mean))
    }))

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
  gc_mat_tree <- gc_mat

  gc_mat
}


