#' Load gene-gene Wasserstein matrix
#'
#' @param dir.path Path to gene-gene Wasserstein matrix
#' @param file.name Name of the matrix file
#' 
#'
#' @return
#' Gene-gene Wasserstein matrix
#'
#'
#' @export
#'
#' @examples
#'

LoadGeneDistMat <- function(dir.path, 
                         file_name = "emd.csv"){
  emd_mat <- as.matrix(read.csv(paste0(dir.path, file_name), header = F))
  genes <- read.csv(paste0(dir.path, "gene_names.csv"), header = F)[,1]
  rownames(emd_mat) <- genes
  colnames(emd_mat) <- genes
  emd_mat
}