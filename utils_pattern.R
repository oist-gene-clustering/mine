##################################################################
#                                                                #
#      MINE ALGORITHM: COMPUTATION OF CELL PATTERNS              #
#                                                                #
##################################################################

##############################
# Binary logistic regression #
##############################

#' Perform binary logistic regression 
#'
#' Performs a binary logistic regression on the training data. Also, if specified, performs validation tests.
#'
#' @param p number of features
#' @param levelList list of the two levels
#' @param trainLevel level list on which the regression should be performed
#' @param trainValue numeric values on which the regression should be performed
#' @param validLevel validation set for levels
#' @param validValue validation set for values
#' @param validation should validation tests be performed?
#' @return the model obtained by the binary logistic regression
#' 
#' @importFrom stats glm
#' @importFrom stats predict
#' @importFrom stats mean
#' 
#' @export
trainLevelLR <- function(p, levelList, trainLevel, trainValue, validLevel=NULL, validValue=NULL, validation = T) {
  ## To lessen the effects of multicollinearity                   ##
  ## Scale "independent variables"                                ##
  tsc <- scale(trainValue)
  missingLines <- sapply(tsc, function(x) all(is.na(x)))
  trainValue[!missingLines] <- tsc[!missingLines]
  df <- data.frame(Level=trainLevel, Value=trainValue)
  model <- glm(Level ~ Value, family=binomial(link='logit'), data=df)
  if (validation & !is.null(validLevel) & !is.null(validValue)) {
    nd <- data.frame(Level=validLevel, Value=scale(validValue))
    validationResults <- ifelse(predict(model, newdata=nd, type="response") > 0.5, levelList[1], levelList[2])
    misClassificError <- mean(validationResults != validLevel)
    acc <- 1-misClassificError
    return(list(model, acc))
  }
  return(model)
}

############################################
# Cross-validation: leave-one-out strategy #
############################################

#' Perfom cross validation
#'
#' Performs cross validation on the data to determine the best model with binary logistic regression.
#'
#' @param p number of features
#' @param k integer: cutting the dataset in k parts
#' @param is_level vector with categorical values for each element of the dataset
#' @param values vector with numeric values for each element of the dataset
#' @return the training set that gives the best model: highest accuracy value -that is the coefficient of correctly classified data, 
#' the associated accuracy value and model
#' 
#' @importFrom caret createFolds
#' 
#' @export
crossValidation <- function(p, k, levelList, is_level, values, it=0) {
  sets <- createFolds(1:p, k=k, list=T, returnTrain = F)
  validset <- (1:p)[!(1:p %in% sets[[1]])]
  tmp <- trainLevelLR(p, levelList, is_level[sets[[1]]], values[sets[[1]]], is_level[validset], values[validset])
  model <- tmp[[1]]
  acc <- tmp[[2]]
  set <- sets[[1]]
  tryCatch(
  for (i in 1:length(sets)) {
    validset <- (1:p)[!(1:p %in% sets[[i]])]
    trainLevel <- is_level[sets[[i]]]
    trainValue <- values[sets[[i]]]
    validLevel <- is_level[validset]
    validSet <- values[validset]
    if (length(unique(trainLevel)) == 1) {
      missingLevel <- levelList[levelList != unique(trainLevel)][1]
      trainLevel <- c(trainLevel, missingLevel)
      idx <- match(missingLevel, validSet)[1]
      trainValue <- c(trainValue, validSet[idx])
      y <- rep(T, length(validSet))
      y[idx] <- F
      validSet <- validSet[y]
      validLevel <- validLevel[y]
    }
    tmp <- trainLevelLR(p, levelList, trainLevel, trainValue, validLevel, validSet)
    nmodel <- tmp[[1]]
    nacc <- tmp[[2]]
    nset <- sets[[1]]
    if (nacc > acc) {
      acc <- nacc
      model <- nmodel
      acc <- nacc
      set <- sets[[i]]
    }
    }, error=function(e) {
    if (it < maxit) {cat("Restart CV...", "\n"); return(crossValidation(p, k, levelList, is_level, values, it=it+1))}
    cat("Too many tries...", "\n")
    return(NULL)
    })
  return(list(set, acc, model))
}

#################
# Main function #
#################

#' Get pattern 
#'
#' Gets cell patterns for a given dataset, and models for dropout and mild-expressiveness for genes.
#'
#' @param sca SingleCellAssay that contains the trimmed log-normalized matrix gene expression matrix
#' @param res MAST thresholdCountMatrix object associated with @sca
#' @param p number of features/genes
#' @param m number of conditions/cells
#' @return both models, the pattern matrix, the pattern-oriented threshold matrix and a vector 
#' indicating the mild-expressiveness of a gene in each cell
#' 
#' @importFrom SummarizedExperiment assay
#' 
#' @export
getPatterns <- function(sca, res, p, m) {
  ## Dropout events depend on the considered cell and the expression value of the gene across the cells (cf. report) ##
  value_dropout <- apply(assay(sca), 2, function(x) if (sum(x) == 0) rep(0, p) else x/sum(x))
  
  is_dropout <- rep(0, p)
  idx <- rowSums(value_dropout) == 0
  is_dropout[idx] <- "dropout"
  is_dropout[!idx] <- "not_dropout"
  
  ## Get the gene expression thresholds ##
  threshold <- sapply(as.character(res$bin), function(bin) res$cutpoint[bin])
  value_mild <- assay(sca)
  is_mild <- apply(value_mild, 2, function(x) ifelse(x > threshold, "high", "mild"))
  
  ## Matrix with column "Dropout" and column "Median" (log2(TPM+1) value across all cells) ##
  ## Initialize the training and validation sets                                           ##
  ## Perform binary logistic regression to fit dropout distribution                        ##
  modelDropout <- crossValidation(p, k=8, unique(is_dropout), is_dropout, rowSums(value_dropout))[[3]]
  if (is.null(modelDropout)) return(NULL)
  
  start <- clock("Part (I): Dropout model fitted")
    
  ## Matrix with column "Mild" and column "Value" (log2(TPM+1) value in each cell) ##
  ## Initialize the training and validation sets                                   ##
  ## Perform binary logistic regression to fit mild-expressiveness distribution    ##
  modelMild <- lapply(1:m, function(j) 
    crossValidation(p, k=3, unique(c(is_mild)), is_mild[, j], unlist(value_mild[, j]))[[3]])
  if (any(sapply(modelMild, is.null))) return(NULL)
  
  start <- clock("Part (I): Expressiveness model fitted")
    
  ## Resulting pattern matrix ##
  P <- is_mild == "high"
  P[P] <- 1
  
  ## Resulting pattern-oriented threshold matrix (see report) ##
  Q <- matrix(sapply(1:m, function(j) ifelse(P[, j] == 1, -threshold[j], threshold[j])), byrow=F, nrow=p)
  
  start <- clock("Part (I): Pattern matrices")
  
  ## Get the two highest peaks, if they exist ##
  peaks <- lapply(as.character(res$bin), function(bin) 
    if (length(res$peaks[bin][[1]]) == 1) rep(res$peaks[bin][[1]][1, 1], 2)
    else {
      if (length(res$peaks[bin][[1]]) == 0) 
        rep(NA, 2) 
      else 
        res$peaks[bin][[1]][1:2, 1]
      }
    )
  
  return(list(modelDropout, modelMild, P, Q, is_mild, p, peaks))
}