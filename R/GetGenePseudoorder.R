#' Order genes along a given trajectory
#'
#' @param dist.mat matrix; Gene-gene Wasserstein distance matrix (symmetric)
#' @param subset vector; Genes in a given trajectory
#' @param max.id integer; Index of the terminal gene
#'
#'
#' @return
#' Pseudoorder of genes along a given trajectory
#'
#'
#' @export
#'
#' @examples
#'
#'


GetGenePseudoorder <- function(dist.mat,
                                  subset,
                                  max.id = NULL){
  emd <- dist.mat[subset, subset]
  #dm.emb <- RunDM_emd_v2(emd, PLOT = F)
  dm.emb <- diffusion.map(emd)$diffu.emb

  pseudoorder <- rank(dm.emb[,2])
  names(pseudoorder) <- rownames(dm.emb)
  if (!is.null(max.id)){
    if (pseudoorder[max.id] <= max(pseudoorder)/2){
      pseudoorder <- rank(-dm.emb[,2])
      names(pseudoorder) <- rownames(dm.emb)
    }
  }

  pseudoorder.all <- rep(0, nrow(dist.mat))
  names(pseudoorder.all) <- rownames(dist.mat)
  pseudoorder.all[names(pseudoorder)] <- pseudoorder

  pseudoorder.all
}
