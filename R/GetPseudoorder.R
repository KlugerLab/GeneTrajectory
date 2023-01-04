GetPseudoorder <- function(dist_mat, 
                                  subset,
                                  max.id = NULL){
  emd <- dist_mat[subset, subset]
  dm_emb <- RunDM_emd_v2(emd, PLOT = F)
  
  pseudoorder <- rank(dm_emb[,2])
  names(pseudoorder) <- rownames(dm_emb)
  if (!is.null(max.id)){
    #if(names(pseudoorder)[which.max(pseudoorder)] != max.id){
    if (pseudoorder[max.id] <= max(pseudoorder)/2){
      pseudoorder <- rank(-dm_emb[,2])
      names(pseudoorder) <- rownames(dm_emb)
    }
  }
  
  pseudoorder_all <- rep(0, nrow(dist_mat))
  names(pseudoorder_all) <- rownames(dist_mat)
  pseudoorder_all[names(pseudoorder)] <- pseudoorder
  
  pseudoorder_all
}