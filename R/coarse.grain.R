#' Base function called by CoarseGrain
#'
#' @param cell.embedding matrix; Cell embedding for kNN clustering
#' @param gene.expression matrix; Cell-by-gene count matrix (library normalized)
#' @param graph.dist matrix; Cell-cell graph distance matrix
#' @param N integer; Number of meta-cells after coarse-graining
#' @param random.seed; Random seed for kNN clustering
#'
#' @return
#' Returns a list
#'
#'
#' @export
#'
#' @examples
#'
#'

coarse.grain <- function(cell.embedding, 
                         gene.expression, 
                         graph.dist,
                         N = 1000, 
                         random.seed = 1){
  
  message("Run k-means clustering")
  set.seed(random.seed)
  km.res <- stats::kmeans(cell.embedding, N)
  
  message("Coarse-grain matrices")
  
  #####KNN graph
  KNN.membership.mat <- matrix(0, nrow = N, ncol = nrow(cell.embedding))
  for (i in 1:ncol(KNN.membership.mat)){
    KNN.membership.mat[km.res$cluster[i], i] <- 1
  }
  KNN.membership.mat <- KNN.membership.mat/apply(KNN.membership.mat, 1, sum)
  
  #####Coarse-grain the gene expression matrix
  gene.expression.updated <- biclust::binarize(KNN.membership.mat, 0) %*% gene.expression
  
  #####Coarse-grain the cost matrix for the calculation of EMD
  graph.dist.updated <- KNN.membership.mat %*% graph.dist %*% t(KNN.membership.mat)
  
  #####Collect the output 
  output <- list()
  output[["gene.expression"]] <- gene.expression.updated
  output[["graph.dist"]] <- graph.dist.updated
  output[["features"]] <- colnames(gene.expression.updated)
  
  output
}