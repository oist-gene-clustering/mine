################################################################
#                                                              #
#               MINE ALGORITHM CORE CODE                       #
#                                                              #
################################################################

library(MAST)
library(doParallel)
library(foreach)
library(caret)
library(Matrix)
library(corpcor)
library(fastcluster)
library(mixtools)
library(mnormt)

path="~/Documents/algorithm/"
sourceIt <- function(file) source(paste(path, file, sep=""))
loadIt <- function(file) load(paste(path, file, sep=""))

sourceIt("heap.R")
sourceIt("unionFind.R")
sourceIt("hashTable.R")

sourceIt("utils.R")
sourceIt("utils_pattern.R")
sourceIt("utils_computation.R")
sourceIt("utils_clustering.R")

###########
# Options #
###########

options(expressions = 5e5, warn=-1)
ncores <- detectCores()-1
toclean <- NULL
lambda_0=0.1

#############
# Main code #
#############

#' Cell clustering of gene expression matrix
#'
#' Clusters cell according to gene expression matrix.
#'
#' @param M un- or normalized gene expression matrix
#' @param d threshold >= 0.5 for cell similarity: default 0.5: user-selected
#' @param d_neighbors threshold for pattern merging: default d
#' @param dprime threshold for feature selection: default 0.9
#' @param dpattern threshold for pattern resemblance: default 0
#' @param f frequency-based gene trimming criterion: default 0.7
#' @param plotDendrogram should a cell hierarchy be plotted? default TRUE
#' @param normalized is matrix M normalized? default FALSE
#' @param returnMatrix returns the trimmed matrix
#' @param returnCopulas returns the fitted copulas
#' @param returnDistMatrix returns the computed distance matrix
#' @return clustering: vector of labels for each cell/condition
#'
#' @importFrom MAST thresholdSCRNACountMatrix
#' @importFrom parallel makeCluster
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel clusterExport
#' @importFrom SummarizedExperiment assay
#' @importFrom stats rnorm
#' @importFrom stats dendrogram
#' @importFrom stats hclust
#' @importFrom stats dist
#' @importFrom corpcor cor.shrink
#' 
#' @export
mine <- function(M, d=0.5, d_neighbors=d, dprime=0.9, dpattern=0.9, f=0.7,  
                 plotDendrogram=T, normalized=F, returnMatrix=T, returnCopulas=T, returnDistMatrix=T) {
  
    #________________________________________________#
    #       Initialization of parameters             #
    #________________________________________________#
  
    ## Number of conditions ##
    m <- dim(M)[2]
    maxit=10
    
    #____________________________________________________________________________________________________#
    #                           FIND MOST PROBABLE PATTERNS                                              #
    # <=> Classification for each cell of every gene into dropout, mildly-expressed and highly-expressed #
    # and filtering of non-informative genes (freq > f or variance = 0 + dimension reduction)            #
    #____________________________________________________________________________________________________# 
    
    start <- clock("Start")
    p <- nrow(M)
    
    print(p)
    
    ## Trim and filter the gene expression matrix ##
    sca <- filterConvertData(M, f, normalized, dprime)
    p <- nrow(assay(sca))
    
    print(p)
    
    ## Compute the gene expression thresholds ##
    res <- thresholdSCRNACountMatrix(assay(sca), nbins=floor(p/100), min_per_bin=2)
    assay(sca) <- res$counts_threshold
    p <- nrow(assay(sca))
    tmp <- getPatterns(sca, res, p, m)
    if (is.null(tmp)) return(NULL)
    modelDropout <- tmp[[1]] 
    modelMild <- tmp[[2]] 
    assay(sca) <- assay(sca)
    P <- tmp[[3]]
    Q <- tmp[[4]]
    is_mild <- tmp[[5]]
    p <- tmp[[6]]
    peaks <- tmp[[7]]

    print(p)
    
    start <- clock("End Part (I): patterns", start)
    
    #__________________________________________________________________________________________________________________________#
    #                       COMPUTATION OF THE PROBABILITY DISTRIBUTION                                                        #
    # <=> Computation of parameters for marginal distributions (Continuous approximation of Poisson-Negative Binomial (CPNB)   #
    #     mixture Model) for each (gene, cell) and for bimodal Gaussian copula for each cell                                   #
    #__________________________________________________________________________________________________________________________# 
    
    ## Compute some of the margin parameters ##
    tmp <- computeMarginParam(p, m, P, modelDropout, sca, is_mild)
    mu <- tmp[[1]]
    delta <- tmp[[2]]
    pdrArray <- tmp[[3]]
    indices_high <- tmp[[4]]
    indices_mild <- tmp[[5]]
    
    start <- clock("Part (II): margin parameters estimation", start)
    
    ## Data with small gaussian noise                                       ##
    ## in order not to apply the margin parameters on data I fitted them on ##
    tdata <- assay(sca) + matrix(rnorm(length(assay(sca))), ncol=m, nrow=p)
    for (j in 1:m) {
        for (i in 1:p) {
          tdata[i,j] <- pmixture(tdata[i, j], P[i, j], pdrArray[i], lambda_0, mu[i], delta[i])
          if (is.nan(tdata[i, j])) tdata[i, j] <- NA
          if (!is.na(tdata[i, j])) tdata[i, j] <- bisectionInv(pnorm, tdata[i, j])
        }
    }
    tdata[is.infinite(tdata) & tdata < 0] <- -5
    tdata[is.infinite(tdata) & tdata > 0] <- 5
    
    start <- clock("Part (II): data transformation", start)
    
    missingLines <- apply(tdata, 1, function(x) any(is.na(x)))
    lines <- 1:p
    lines[missingLines] <- NA
    lines[!missingLines] <- 1:(p-length(missingLines[missingLines]))
    
    getIndList <- function(idxList) return(lapply(idxList, function(x) {xx <- sapply(x, function(e) lines[e]); xx[!is.na(xx)]}))
    indices_high <- getIndList(indices_high)
    indices_mild <- getIndList(indices_mild)
    tdata <- tdata[!missingLines, ]
    p <- dim(tdata)[1]
    
    print(p)
    
    cvm <- cor.shrink(t(tdata))[, ]
    
    start <- clock("Part (II): generating training data", start)
    
    cl <- makeCluster(ncores, type="PSOCK", outfile="")
    clusterEvalQ(cl=cl, library(SummarizedExperiment))
    clusterEvalQ(cl=cl, library(stats))
    clusterEvalQ(cl=cl, library(mnormt))
    clusterExport(cl=cl, varlist=ls(.GlobalEnv), envir=.GlobalEnv)
    copulas <- parLapply(cl, 1:m, function(j) {
      cat(sprintf("Copula %d", j), "\n", file="")
      tmp <- initParam(m, p, j, tdata, modelMild[[j]], list(indices_high[[j]], indices_mild[[j]]), P, mu, delta, pdrArray, peaks)
      bimodalNormalCopula(p, j, assay(sca), c(tmp[[2]], 1-tmp[[2]]), tmp[[3]],
                          rep(list(cvm), length(tmp[[3]])), P, pdrArray, lambda_0=lambda_0, mu, delta)
    })
    stopCluster(cl)
    
    start <- clock("End Part (II): bimodal copulas", start)
    
    #______________________________________________________________________________________________#
    #                       PATTERN MERGING                                                        #
    # <=> Compute the similarity between each pair of cells                                        #
    #______________________________________________________________________________________________# 

    ## Create a priority queue of the mergeable pairs, ordered by decreasing similarity ##
    mergeablePairs <- createHeap()
    distMatrix <- diag(rep(0, m))
    for (k in 1:(m-1)) {
      for (l in (k+1):m) {
        ## Compute similarity between cell k and cell l ##
        sim <- slm(k, l, copulas, P, Q, dpattern)
        distMatrix[l][k] <- 1-sim
        distMatrix[k][l] <- distMatrix[l][k]
        ## If cells k and l are mergeable, then add them to the queue ##
        if (sim > d) pushHeap(list(k, l), sim, mergeablePairs)
      }
    }
    
    start <- clock("End Part (III): pattern merging", start)
    
    print(mergeablePairs)
    
    return(list(copulas, simValues, assay(sca), tmp[[3]]))
    
    #_________________________________________________________________________________________#
    #                       CLUSTERING                                                        #
    #_________________________________________________________________________________________# 
    
    ## Initialize the Union Find structure that contains the clustering      ##
    createUF(1:m, "uf")
    ## Number of neighbors that need to agree on the merging of two clusters    ##
    ## At each iteration i, the number of neighbors needed for a cluster C_i is ##
    ## floor(sqrt(#C_i-1)) (my definition of neighbor does not include the      ##
    ## considered node itself)                                                  ##
    ## At first, since there is only one element in each cluster, no neighbor   ##
    ## is required                                                              ##
    K <- rep(0, m)
    clusterSizes <- rep(1, m)
    ## Iterate the following steps until convergence of the clusters         ##
    while (T) {
      copyUF("uf_old", "uf")
      queue <- mergeablePairs
      while (!isEmpty(queue)) {
        pair <- extractMax(queue)
        k <- pair[[1]][[1]]
        l <- pair[[1]][[2]]
        ## Check the agreement to the merging of K neighbors for each cluster ##
        if (agreeMerge(uf, d_neighbors, K[k], K[l], k, l, distMatrix, m)) {
          unionUF(k, l, "uf")
          clusterSizes[k] <- clusterSizes[k]+clusterSizes[l]
          clusterSizes[l] <- clusterSizes[k]
          K[k] <- floor(sqrt(clusterSizes[k]-1))
          K[l] <- K[k]
        }
      }
      if (noChangeUF("uf", "uf_old", m)) break
    }
    if (plotDendrogram) plot(hclust(as.dist(scale(distMatrix)), method="average"))
    
    start <- clock("End Part (IV): clustering", start)
    
    ## Cleaning ##
    destroyUF("uf_old")
    houseKeeping()
    
    res <- list(get("uf", envir = sys.frame(0)), assay(sca))
    destroyUF("uf")
    if (returnMatrix) res <- append(res, list(M))
    if (returnCopulas) res <- append(res, list(copulas))
    if (returnDistMatrix) res <- append(res, list(distMatrix))
    return(res)
}
  