#' Coarse-graining
#'
#' @param object Seurat object
#' @param features vector;
#' @param N integer;
#' @param assay string;
#' @param reduction string;
#' @param dims vector;
#' @param slot string;
#' @param random.seed integer;
#'
#' @return
#' Returns a list
#'
#' @export
#'
#' @examples
#'
#'

CoarseGrain <- function(object,
                        graph.dist,
                        features,
                        N = 1000,
                        assay = "RNA",
                        reduction = "dm",
                        dims = 1:5,
                        slot = "data",
                        random.seed = 1){
  
  data.S <- object
  cell.embedding <- data.S@reductions[[reduction]]@cell.embeddings[, dims]
  gene.expression <- t(as.matrix(GetAssayData(data.S, slot = slot, assay = assay)[features, ]))
  cg.output <- coarse.grain(cell.embedding = cell.embedding,
                            gene.expression = gene.expression,
                            graph.dist = graph.dist,
                            N = N,
                            random.seed = random.seed)
  
  cg.output
}


