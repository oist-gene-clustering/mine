####################################################################
#                                                                  #
#               MINE ALGORITHM CLUSTERING PART                     #
#                                                                  #
####################################################################

###################################
# Pairwise similarity computation #
###################################

#' Get similarity value for two cells
#'
#' Gets the similarity value for two cells, see report.
#'
#' @param k first cell number
#' @param l second cell number
#' @param copulas copula list for all cells in data
#' @param P gene expression pattern matrix
#' @param Q pattern-oriented gene expression threshold matrix
#' @return similarity value for cells @k and @l
#' 
#' @importFrom stats mean
#' 
#' @export
slm <- function(k, l, copulas, P, Q, dpattern) {
  computeProb <- function(j) return(pBCopula(Q[, j], copulas[[j]]))
  ## Identity ##
  if (k == l) return(1)
  ## Complete dissimilarity ##
  if (mean(P[, k] == P[, l]) <= dpattern) return(0)
  
  pl <- computeProb(l)
  pk <- computeProb(k)
  print(paste0(k, " ", l, " ", pl, " ", pk, " ", pl*pk))
  
  ## Cell mutual independence is assumed ##
  return(pl*pk)
}

###################################
# Clustering: merging agreement  #
###################################

#' Agreement for same-cluster neighbors
#'
#' Returns boolean for agreement of same-cluster neighbors.
#'
#' @param m number of conditions
#' @param uf Union Find structure that stores the clusters
#' @param j the index of the considered condition
#' @param root_j root of the cluster to which @j belongs
#' @param d the similarity threshold for neighbors
#' @param Kj threshold number of neighbors needed for merging for the cluster associated with @j
#' @param distMatrix distance matrix between conditions
#' @return boolean for -dis-agreement from one cluster
#' 
#' @export
agree <- function(m, uf, j, root_j, d, Kj, distMatrix) {
  ## list of all neighbors in cluster ##
  neighbors_j <- (1:m)[findUF((1:m), uf) == root_j]
  neighbors_j <- neighbors_j[neighbors_j != j]
  ## list of agreeing -to the merging- neighbors in the cluster ##
  agreeing_nj <- neighbors_j[1-distMatrix[j, neighbors_j] > d]
  return((length(agreeing_nj) >= Kj))
}

#' Agreement for mergeing two clusters
#'
#' Returns boolean for agreement from both clusters.
#'
#' @param uf union find structure containing the clusters
#' @param d threshold probability for similarity between two cells
#' @param Kk threshold number of neighbors needed for merging first cluster
#' @param Kl threshold number of neighbors needed for merging second cluster
#' @param k first cell number
#' @param l second cell number
#' @param distMatrix distance matrix for cells = 1-similarity matrix
#' @param m number of cells
#' @return boolean for -dis-agreement from both clusters
#' 
#' @export
agreeMerge <- function(uf, d, Kk, Kl, k, l, distMatrix, m) {
  tmp <- sapply(1:m, function(j) findUF(j, uf))
  root_k <- findUF(k, uf)
  root_l <- findUF(l, uf)
  if (root_k == root_l) return(FALSE)
  return(agree(m, uf, k, root_k, d, Kk, distMatrix) 
         & agree(m, uf, l, root_l, d, Kl, distMatrix))
}