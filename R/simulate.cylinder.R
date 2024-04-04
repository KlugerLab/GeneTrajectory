#' Simulate concurrent independent linear and cyclic gene processes
#'
#' @param N_cells integer; Number of cells
#' @param N_genes integer vector; Number of genes associated with the linear and cyclic processes
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
#' @export simulate.cylinder
#'
#' @examples
#'

simulate.cylinder <- function(N_cells = 5000,
                              N_genes = rep(500,2),
                              model = "poisson",
                              meanlog = 0,
                              sdlog = 0.25,
                              scale = 25,
                              seed = 1,
                              maxT = 15,
                              sort = TRUE,
                              sparsity = 0.1,
                              theta = 10){
  gc_mat_linear <- simulate.linear(N_cells = N_cells,
                                   N_genes = N_genes[1],
                                   model = model,
                                   meanlog = meanlog,
                                   sdlog = sdlog,
                                   scale = scale,
                                   seed = seed,
                                   maxT = maxT,
                                   sort = sort,
                                   sparsity = sparsity,
                                   theta = theta)
  gc_mat_cyclic <- simulate.cyclic(N_cells = N_cells,
                                     N_genes = N_genes[2],
                                     model = model,
                                     meanlog = meanlog,
                                     sdlog = sdlog,
                                     scale = scale,
                                     seed = seed,
                                     maxT = maxT,
                                     sparsity = sparsity,
                                     theta = theta)

  set.seed(1)
  idx = sample(1:ncol(gc_mat_cyclic))
  gc_mat_cylinder <- rbind(gc_mat_linear, gc_mat_cyclic[,idx])
  colnames(gc_mat_cylinder) <- paste0(colnames(gc_mat_linear), "|", colnames(gc_mat_cyclic[,idx]))

  gc_mat_cylinder
}

