########################################################################
#                                                                      #
#               UNION FIND STRUCTURE IMPLEMENTATION                    #
#                                                                      #
########################################################################

source("~/Documents/algorithm/utils.R")

#########################################
# Initialization of Union Find & utils. #
#########################################

#' Create Union Find
#'
#' Creates a Union Find structure from an array in the global environment --side effect.
#'
#' @param array array of elements
#' @param uf union find name
#' 
#' @export
createUF <- function(array, uf) {
  assignit(paste(uf, "_parent", sep = ""), array)
  assignit(paste(uf, "_rank", sep = ""), rep(1, length(array)))
  return(NULL)
}

#' Get rank array of Union Find
#'
#' Gets the ranks of a Union Find structure.
#'
#' @param uf union find name
#' @return ranks of @uf
#' 
#' @export
getRanks <- function(uf) 
  return(getit(paste(uf, "_rank", sep = "")))

#' Get parent array of Union Find
#'
#' Gets the parent/root array of a Union Find structure.
#'
#' @param uf union find name
#' @return parents of @uf
#' 
#' @export
getParents <- function(uf) 
  return(getit(paste(uf, "_parent", sep = "")))

#' Get rank of element
#'
#' Gets the rank of an element in the Union Find structure.
#'
#' @param i number of element in UF
#' @param uf union find name
#' @return rank of @i in @uf
getRank <- function(i, uf) 
  return(getRanks(uf)[i])

#' Get parent of element
#'
#' Gets the parent/root of an element in the Union Find structure.
#'
#' @param i number of element in UF
#' @param uf union find name
#' @return parent of @i in @uf
#' 
#' @export
getParent <- function(i, uf) 
  return(getParents(uf)[i])

#' Assign parent to element
#'
#' Assigns a parent to an element in the Union Find structure --side effect.
#'
#' @param i number of element in UF
#' @param x parent to assign
#' @param uf union find name
#' 
#' @export
assignParent <- function(i, x, uf) {
  parent <- getParents(uf)
  parent[i] <- x
  assignit(paste(uf, "_parent", sep = ""), parent)
  return(NULL)
}

#' Assign rank to element
#'
#' Assigns a rank to an element in the Union Find structure --side effect.
#'
#' @param i number of element in UF
#' @param x rank to assign
#' @param uf union find name
#' 
#' @export
assignRank <- function(i, x, uf) {
  rank <- getRanks(uf)
  rank[i] <- x
  assignit(paste(uf, "_rank", sep = ""), rank)
  return(NULL)
}

#' Check change between union find structures
#'
#' Checks if the two union find structures are equal.
#'
#' @param uf_old first UF structure
#' @param uf_new second UF structure
#' @param k size of both UF structures
#' @return TRUE iff the two structures are the same
#' 
#' @export
noChangeUF <- function(uf_old, uf_new, k) {
  parentUFOLD <- getParents(uf_old)
  parentUFNEW <- getParents(uf_new)
  if (!(length(parentUFNEW) == length(parentUFOLD))) return(F)
  for (i in 1:k) {
    tmp <- findUF(i, uf_old)
    tmp <- findUF(i, uf_new)
  }
  return(all(parentUFOLD == parentUFNEW))
}

#' Copy UF structure
#'
#' Copies an UF structure --side effect.
#'
#' @param uf_old name of the UF structure to copy
#' @param uf_new sname of the UF structure copy
#' 
#' @export
copyUF <- function(uf_old, uf_new) {
  assignit(paste(uf_new, "_parent", sep = ""), getParents(uf_old))
  assignit(paste(uf_new, "_rank", sep = ""), getParents(uf_old))
  return(NULL)
}

#' Destroy UF structure
#'
#' Destroys an UF structure and does garbage collection --side effect.
#'
#' @param uf name of the UF structure to destroy
#' 
#' @export
destroyUF <- function(uf) {
  rmit(as.character(sprintf("%s_parent", uf)))
  rmit(as.character(sprintf("%s_rank", uf)))
  gc()
  return(NULL)
}

###################
# Main operations #
###################

#' Find parent of element in UF structure
#'
#' Finds the parent of an element in the UF structure.
#'
#' @param i number of element in the UF structure
#' @param uf name of the UF structure
#' @return parent of element @i
#' 
#' @export
findUF <- function(i, uf) {
  if (!(getParent(i, uf) == i)) assignParent(i, findUF(getParent(i, uf), uf), uf)
  return(getParent(i, uf))
}

#' Merge two groups in UF structure
#'
#' Merges two groups of elements using path compression in the UF structure --side effect.
#'
#' @param i number of first element in the UF structure
#' @param j number of second element in the UF structure
#' @param uf name of the UF structure
#' 
#' @export
unionUF <- function(i, j, uf) {
  root_i <- findUF(i, uf)
  root_j <- findUF(j, uf)
  if (root_i == root_j) return(NULL)
  if (getRank(root_i, uf) < getRank(root_j, uf)) assignParent(root_i, root_j, uf)
  if (getRank(root_j, uf) < getRank(root_i, uf)) assignParent(root_j, root_i, uf)
  else {
    assignParent(root_j, root_i, uf) 
    assignRank(root_i, getRank(root_i, uf) + 1, uf) 
  }
  return(NULL)
}
  