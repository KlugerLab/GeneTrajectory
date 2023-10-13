#' Run Diffusion Map on a Seurat object
#'
#' @param object Seurat object
#' @param reduction string; Dimensionality reduction to use, default: 'dm'
#' @param dims integer; Which dimensions to use as input features for kNN graph construction
#' @param K integer; Adaptive kernel bandwidth for each point set to be the distance to its `K`-th nearest neighbor.
#'
#' @return
#' Returns a cell-cell graph distance matrix
#'
#'
#' @export
#'
#' @examples
#'
#'

GetGraphDistance <- function(object,
                             reduction = "dm",
                             dims = 1:5,
                             K = 10){
  data.S <- object
  cell.embedding <- data.S@reductions[[reduction]]@cell.embeddings[, dims]
  
  message("Constructing kNN graph")
  knn.result <- FNN::get.knn(cell.embedding, k = K)
  KNN.adj.mat <- matrix(0, nrow = nrow(cell.embedding), ncol = nrow(cell.embedding))
  for (i in 1:nrow(KNN.adj.mat)){
    KNN.adj.mat[i, knn.result$nn.index[i,]] <- 1
  }
  
  message("Constructing graph distance matrix")
  knn.graph <- igraph::graph_from_adjacency_matrix(KNN.adj.mat, mode = "undirected")
  graph.dist.mat <- igraph::distances(knn.graph)
  if (is.infinite(max(graph.dist.mat))) stop("The cell-cell kNN graph has disconnected components. Please increase K.")
  message(sprintf("The largest graph distance is %s", max(graph.dist.mat)))
  
  graph.dist.mat
}