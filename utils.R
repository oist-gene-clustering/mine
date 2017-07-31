##################################################################
#                                                                #
#               UTILS FOR MINE ALGORITHM                         #
#                                                                #
##################################################################

###################
# Type conversion #
###################

#' Convert type matrix
#'
#' Converts a matrix to a certain type-valued matrix.
#'
#' @param obj matrix
#' @param t character sting definying type
#' @return t valued matrix
#' 
#' @export
convertTypeMatrix <- function(obj, t) 
  return(apply(obj, 2, function(x) {storage.mode(x) <- t; x}))

#' Convert type vector
#'
#' Converts a vector to a certain type-valued vector.
#'
#' @param obj vector
#' @param t character sting definying type
#' @return t valued vector
#' 
#' @export
convertTypeVector <- function(obj, t) 
  return(sapply(obj, function(x) {storage.mode(x) <- t; x}))

#################
# Normalization #
#################

#' log2 quasi-TPM+1 normalization of matrix
#'
#' Converts a matrix into a log2 quasi-TPM+1 normalized matrix. It does NOT normalize for gene length,
#' and it accounts for sequencing depth.
#'
#' @param M raw counts matrix
#' @return normalized matrix
normalizeData <- function(M) 
  return(as.matrix(apply(M, 2, function(x) log(x/sum(x)*10**6+1)/log(2))))

#' log2+1 transformation of matrix
#'
#' Converts a matrix into a log2 +1 transformed matrix.
#'
#' @param M normalized matrix
#' @return regularized matrix
#' 
#' @export
normalizeLogData <- function(M) 
  return(as.matrix(log(M+1)/log(2)))

######################################################
# Filtering and conversion into SingleCellAssay type #
######################################################

#' feature selection for features
#'
#' Applies custom feature selection to remove genes that are too positively correlated to others.
#'
#' @param M normalized nonnegative matrix
#' @param dprime threshold for similarity between features
#' @return list of indices of kept lines
#' 
#' @importFrom stats cor
#' @importFrom stats cutree
#' @importFrom stats mean
#' @importFrom stats var
#' @importFrom stats dist
#' @importFrom caret findCorrelation
#' @importFrom fastcluster hclust
#' 
#' @export
featureSelection <- function(M, dprime, useCaret=T) {
  ## Delete features of null variance ##
  MM <- M[apply(M, 2, function(x) !(var(x) == 0 | sum(x) == 0)), ]
  MM <- scale(MM)
  MM <- na.omit(MM)[,]
  corM <- cor(t(MM), method = "spearman")
  p <- nrow(MM)
  if (useCaret) return((1:p)[!(1:p %in% caret::findCorrelation(corM, cutoff=dprime))])
  means <- apply(M, 1, mean)
  stdevs <- apply(M, 1, function(x) sqrt(var(x)))
  msm <- means-stdevs
  msp <- means+stdevs
  ## Genes are considered truly redundant iff. they are highly correlated ##
  ## and their mean expression value is roughly similar                   ##
  ms <- matrix(sapply(means, function(x) x >= msm & x <= msp), byrow = T, ncol=p)
  ## maximum distance: 1 ##
  cordist <- (1-abs(corM))/2
  cordist[!ms] = 2
  ## Average linkage method ##
  hc <- fastcluster::hclust(as.dist(cordist), method="average")
  hc <- cutree(hc, h=1-dprime)
  return(match(1:max(hc), hc))
}

#' Convert a matrix into a correct SingleCellAssay
#'
#' Converts a matrix into a correct SingleCellAssay.
#'
#' @param M un- or normalized nonnegative matrix
#' @param f frequency-based filtering criterion 0 < f < 1
#' @param normalized boolean = true iff M is normalized
#' @param dprime the threshold correlation for genes for feature selection
#' @return SingleCellAssay sca where assay@sca = log2-regularized normalized trimmed M
#' 
#' @importFrom MAST FromMatrix
#' @importFrom MAST freq
#' @importFrom stats na.fail
#' @importFrom SummarizedExperiment assay
#' 
#' @export
filterConvertData <- function(M, f, normalized, dprime) {
  ## Delete features where there is at least one missing value ##
  M <- na.omit(M)[,]
  M <- if (normalized) normalizeLogData(M) else normalizeData(M)
  fData = data.frame(primerid=row.names(M), row.names=row.names(M))
  cData = data.frame(wellKey=colnames(M), row.names=colnames(M))
  sca <- FromMatrix(M, cData, fData)
  sca <- sca[which(freq(sca)<f), ]
  sca <- sca[featureSelection(assay(sca), dprime), ]
  return(sca)
}

#################
# GENERAL UTILS #
#################

#' Create value in global environment
#'
#' Creates a value in the global environment --side effect.
#'
#' @param name value name
#' @param value value to assign
#' 
#' @export
assignit <- function(name, value) {
  assign(name, value, envir=sys.frame(0))
  toclean <- c(toclean, get(name, envir=sys.frame(0)))
  return(NULL)
}

#' Get value in global environment
#'
#' Gets the value of a variable in the global environment.
#'
#' @param name value name
#' @return value of variable name in global environment
#' 
#' @export
getit <- function(name) 
  return(get(name, envir=sys.frame(0)))

#' Remove value in global environment
#'
#' Removes a variable in the global environment --side effect.
#'
#' @param name variable name
#' 
#' @export
rmit <- function(name) {
  rm(name, envir=sys.frame(0))
  return(NULL)
}

#' Print current runtime since start
#'
#' Prints the runtime since start.
#'
#' @param msg message to print with the runtime
#' @param start starting time point
#' @return new starting point
#' 
#' @importFrom stats time
#' 
#' @export
clock <- function(msg, start=NULL) {
  if (!is.null(start)) cat(sprintf("%s %2.2f", msg, Sys.time()-start), "\n")
  return(Sys.time())
}

#' Log conversion wrapper
#'
#' Returns the log-value or the value of the result.
#'
#' @param value result
#' @param log boolean indicating if the value should be log-converted
#' @return log or not converted value
#' 
#' @export
retLog <- function(value, log) 
  if (log) return(log(value)) else return(value)

#' Global values cleaner and garbage collector
#'
#' Deletes the global values in argument and performs garbage collecting --side effect.
#' 
#' @export
houseKeeping <- function() {
  for(e in toclean) rm(e)
  gc()
  return(NULL)
}

#' Finds point of inverse function by bisection
#'
#' Finds -quickly- a such as:
#' - f of a belongs to [x-eps, x+eps], 
#' - if y0 is the inf of y such as f of y belongs to [x-eps, x+eps], then |a-y0| < limInf
#'
#' @param f MONOTONOUS NONDECREASING, continuous and univariate function
#' @param x value of the codomain of f
#' @param eps threshold for image value
#' @param step length of initial search interval
#' @param lim stop criterion for bisection
#' @param xmin lower bound of the codomain of f
#' @param xmax upper bound of the codomain of f
#' @return fiber of x by f, or -1 if x is not in the codomain of f, or if it looped
#' 
#' @export
bisectionInv <- function(f, x, eps=1e-7, step=100, lim=1e-7, xmin=0, xmax=1, maxiter=1e5, limInf=1e-4) {
  if (x > xmax | x < xmin) return(-1)
  if (x == xmax) return(Inf)
  if (x == xmin) return(-Inf)
  l <- 0
  r <- step
  while(f(r) < x & maxiter > 0) {r <- r+step; maxiter <- maxiter-1}
  if (maxiter == 0) return(-2)
  while(r-l > lim) {
    m <- (r+l)/2
    if (abs(f(m)-x) < eps) {
      while(abs(f(m)-x) < eps) m <- m-limInf
      return(m+limInf)
    }
    if (f(m) > x) r <- m else l <- m
  }
  return(m)
}

#####################################
# MV normal distribution estimation #
#####################################

#' Initialization of parameters 
#'
#' Initializes the parameters of each of the 2 components.
#'
#' @param m number of conditions
#' @param p number of features
#' @param j index of the considered cell
#' @param data data matrix
#' @param modelMild_j logistic regression model for copula j
#' @param indices_high_mild_j indices for the two components of the mixture: high, mild
#' @param P pattern matrix
#' @param mu mean parameter values for margin distributions
#' @param delta dispersion vector parameter values for margin distributions
#' @param pdrArray mixture rate parameter values for margin distributions
#' @param peaks list of the two highest peaks
#' @return list of initialization of transformed data, mixture rate parameter, mean parameters
#' 
#' @importFrom stats density
#' @importFrom stats rnorm
#' @importFrom stats predict
#' @importFrom stats mean
#' 
#' @export
initParam <- function(m, p, j, tdata, modelMild_j, indices_high_mild_j, P, mu, delta, pdrArray, peaks) {
  if (length(indices_high_mild_j) == 0) return(NULL)
  if (length(indices_high_mild_j) > 1) {
    if (length(indices_high_mild_j[[1]]) == 0 & length(indices_high_mild_j[[2]]) == 0) return(NULL)
    if (length(indices_high_mild_j[[2]]) == 0) indices_high_mild_j <- list(indices_high_mild_j[[1]])
    if (length(indices_high_mild_j[[1]]) == 0) indices_high_mild_j <- list(indices_high_mild_j[[2]])
  }
  kk <- length(indices_high_mild_j)
  ## Empirical generation of data ##
  dens <- lapply(1:kk, function(k) density(tdata[indices_high_mild_j[[k]], j], na.rm=T))
  ## Having number of samples > number of features ##
  N=2*p+1
  q <- lapply(1:kk, function(k) sample(tdata[indices_high_mild_j[[k]], j], N, replace=T) + rnorm(N, 0, dens[[k]]$bw))
  if (kk < 2) q <- list(q[[1]])
  else q <- list(q[[1]], q[[2]])
  q <- lapply(1:kk, function(k) base::t(matrix(q[[k]], p)))
  if (kk < 2) q <- rbind(q[[1]])
  else q <- rbind(q[[1]], q[[2]])
  comp <- lapply(1:kk, function(k) tdata[indices_high_mild_j[[k]], j])
  meancomp <- lapply(1:kk, function(k) mean(tdata[indices_high_mild_j[[k]], ], na.rm=T))
  muCop <- lapply(1:kk, function(k) sapply(1:p, function(i) peaks[[i]][k]))
  for (i in 1:kk) muCop[[i]] <- sapply(muCop[[i]], function(x) if (is.na(x[1])) meancomp[[i]] else x)
  lambda <- predict(modelMild_j, data.frame(Value=mean(tdata[, j])), type="response")[1]
  return(list(q, lambda, muCop))
}

#' Compute -generalized- inverse of matrix
#'
#' Computes -generalized- inverse of matrix.
#'
#' @param rho matrix
#' @param epsilon tolerance to detect zero values
#' @return if @rho is singular, returns Moore-Penrose generalized inverse M: rho*M*rho = rho
#' else returns exact inverse -uses Cholesky factorization if rho is symmetric PD
#' 
#' @export
invCustom <- function(rho, epsilon=sqrt(.Machine$double.eps), evalues=F) 
  return(chol2inv(chol(rho)))

#' Multivariate normal distribution density -cf mvtnorm::dmvnorm
#'
#' Computes the density for one point for normal distribution.
#'
#' @param q vector: already transformed: qnorm @ margin
#' @param mu mean vector parameter
#' @param rho covariance matrix parameter
#' @param log logical: should density be returned in log?
#' @return density of X = x, X having a multivariate normal distribution
#' 
#' @export
dmvnormCustom <- function (q, mu, rho, invrho = NULL) {
  if (is.null(invrho)) invrho <- invCustom(rho)
  x <- x-mu
  if (!is.matrix(x)) x <- matrix(x, ncol=p)
  distval <- rowSums(x %*% invrho * x)
  logdet <- sum(log(eigen(invrho, symmetric=T, only.values=T)$values))
  logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
  return(exp(logretval))
}

# ##########################################################
# # Tests and plots for gene expression distribution assay #
# ##########################################################
# #GOAL: Hypothesis on the gene expression distribution at cell level and at gene level
# 
# # library(cionagedata)
# # library(celegansdata)
# biase <- readRDS("~/Documents/benchmark/biase.rds")
# library(scater)
# library(psych)
# #sceset <- ciona_sceset
# #sceset <- celegans_sceset
# #ccs <- counts(sceset)
# ccs <- fpkm(biase)
# # Accounting for dropouts
# ccs_dropout <- ccs[rowSums(ccs) == 0,]
# ccs_nodropout <- ccs[rowSums(ccs) > 0,]
# 
# # R help for function psych::pairs.panels
# # Adapted from the help page for pairs, pairs.panels shows a scatter
# # plot of matrices (SPLOM), with bivariate scatter plots below the
# # diagonal, histograms on the diagonal, and the Pearson correlation
# # above the diagonal. Useful for descriptive statistics of small
# # data sets.  If lm=TRUE, linear regression fits are shown for both
# # y by x and x by y.  Correlation ellipses are also shown. Points
# # may be given different colors depending upon some grouping
# # variable.
# # 1 < n, m < 11
# 
# # This function plots the dependence between n randomly selected genes for 1 cell
# # => gives a idea of the marginal laws and gene dependence relationship
# depGene1Cell <- function(mat, n=5) {
#   m <- t(mat)
#   return(pairs.panels(m[sample(1:dim(m)[1], 2), sample(1:dim(m)[2], n)]))
# }
# 
# # This function plots the dependence between m randomly selected cells for n randomly selected genes
# # => gives an idea of the dependence relationships between cell copulas, and gene dependence for 1-cell copula
# # (on the diagonal part)
# depCellnGene <- function(mat, n=5, m=5) return(pairs.panels(mat[sample(1:dim(mat)[1], n), sample(1:dim(mat)[2], m)]))
# 
# getList <- function(field, param, interval) {
#   tmp <- unlist(lapply(interval, function(i) if (pData(sceset)[i, field] == param) i))
#   return(tmp[!(is.null(tmp))])
# }
# 
# # This function plots the dependence between m randomly selected cells (should come from the same embryo or cell or stage
# # to be relevant) for n randomly selected genes
# # => gives an idea about dependence between copulas, and gene dependence for the same 1-cell copula for
# #one parameter (cell stage, embryo, etc.)
# depCellnGeneSelect <- function(mat, n=5, m=5, em=NULL, ct=NULL, stage=NULL) {
#   tmp <- 1:dim(mat)[2]
#   if (!is.null(em)) tmp <- getList("Embryo", em, tmp)
#   if (!is.null(ct)) tmp <- getList("Cell.ID", ct, tmp)
#   if (!is.null(stage)) tmp <- getList("EmbryoStage", stage, tmp)
#   return(depCellnGene(mat[, tmp], n, m))
# }
# 
# # This function plots the cell dependency for 1 randomly selected gene (for m randomly selected cells)
# # => gives a idea about the 1-gene marginal distribution
# depCellGene <- function(mat, type="log2(TPM+1)", m=5) {
#   gene <- sample(1:dim(mat)[1], 1)
#   plot(1:dim(mat)[2], mat[gene, ], main=sprintf("Gene expression in the dataset for gene %s", rownames(mat)[gene]),
#          xlab = "Cells", ylab = sprintf("%s counts", type), xaxt='n')
#   axis(side=1, at=1:dim(mat)[2], labels=colnames(mat))
#   return(NULL)
# }
# 
# ##################################################
# # Tests and plots for binary logistic regression #
# ##################################################
# #GOAL: Determine the relationship between dropout events and gene expression values
# #GOAL: Determine the relationship between mild-expressiveness events and gene expression values
# 
# #Tests linearity between dependent variables (dropout or mild-expressiveness) and logits
# #   levels = 1/(1+exp(a*centered_values)) (levels being 0 or 1 vector)
# testLinearity <- function(lvls, values, event) {
#   lvls <- ifelse(lvls == lvls[1], 1, 0)
#   values <- values-mean(values)
#   linearModel <- glm(lvls ~ values, family="binomial")
#   nullModel <- glm(lvls ~ 1, family="binomial")
#   #McFadden R^2 index ~ R^2 coefficient for linear regression
#   #1-logLik(linearModel)/logLik(nullModel)
#   #= library(pscl) pR2(linearModel)[["McFadden"]]
#   res <- list(linearModel$coefficients, 1-logLik(linearModel)/logLik(nullModel), anova(linearModel, test="Chisq"))
#   plot(values, lvls, col="green", main=sprintf("Plotting centered gene expression values against %s events", event),
#        xlab="centered gene expression values", ylab=sprintf("%s events", event))
#   abline(res[[1]][1], res[[1]][2], col="red")
#   return(res)
# }
# 
# #Tests independence between residuals
# testIndependence <- function(lvls, values) {
#   x <- values-mean(values)
#   y <- ifelse(lvls == lvls[1], 1, 0)
#   linearModel <- glm(y ~ x, family="binomial")
#   #plot(linearModel)
#   library(lmtest)
#   yy <- linearModel$residuals[-length(linearModel$residuals)]
#   xx <- linearModel$residuals[-1]
#   model2 = lm(yy ~ xx)
#   residuals <- linearModel$residuals
#   return(list(dwtest(y ~ x), linearModel, chisq.test(abs(linearModel$residuals)), model2, acf(residuals)))
# }

##############################
# Tests and plots for copula #
##############################
#GOAL: See if the fitted copula distribution model is consistent with the data

testCopula <- function(copulas=NULL, dataset=NULL, j=NULL, copulaTest=F, allTest = F) {
  if (allTest) {
    library(psych)
    cl <- makeCluster(ncores, type="FORK")
    x <- as.matrix(parSapply(cl, 1:length(copulas), function(k) rBCopula(1, copulas[k][[1]])))
    stopCluster(cl)
    return(pairs.panels(x))
  }
  if (copulaTest) {
    cl <- makeCluster(ncores, type="FORK")
    x <- parSapply(cl, 1:length(copulas), function(k) rBCopula(1, copulas[k][[1]])[j])
    stopCluster(cl)
    return(hist(x, xlab=sprintf("gene expression values for gene %s", row.names(dataset)[j]), col="grey",
                main="Histogram of gene expression value generated by copula"))
  }
  return(hist(dataset[, j], xlab=sprintf("gene expression values for gene %s", row.names(dataset)[j]), col="grey",
              main=paste("Histogram of gene expression value in", deparse(substitute(dataset)))))
}

retPlot <- function(i,k) {
  for (j in i:k) {
    plot(seq(-2, -1.5, 0.001), sapply(seq(-2,-1.8,0.001), 
                                      function(i) dBCopula(rep(i, dim(MM)[2]), copulas[[j]])), col=palette(rainbow(32))[j])
    par(new=T)
  }
  return(NULL)
}

retHist <- function(i,k) {
  block <- NULL
  for (j in i:k) {
    block <- c(block, sapply(seq(-2,-1.8,0.001), function(i) dBCopula(rep(i, dim(MM)[2]), copulas[[j]])))
  }
  return(hist(block, fill="grey"))
}