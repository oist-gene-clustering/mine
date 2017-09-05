################################################################
#                                                              #
#               MINE ALGORITHM COPULA PART                     #
#                                                              #
################################################################

#############################
# Continuous approximations #
#############################

#_______________#
#   Poisson     #
#_______________# 

#' Continuous approximation of Poisson -CAP- density function
#'
#' Density function for CAP.
#'
#' @param x value
#' @param lambda mean and variance parameter
#' @return for X following this distribution, returns Prob of X = x
#' 
#' @importFrom stats dgamma
#' 
#' @export
dcaP <- function(x, lambda, log=FALSE) 
  return(dgamma(x, shape=lambda, scale=1, log=log))

#' Continuous approximation of Poisson -CAP- cumulative probability function
#'
#' Cumulative probability function for CAP.
#'
#' @param q quantile
#' @param lambda mean and variance parameter
#' @param lower.tail if set to FALSE, returns Prob of X > q
#' @param log.p if set to TRUE, returns the log-probability
#' @return for X following this distribution, returns Prob of X <= q
#' 
#' @importFrom stats pgamma
#' 
#' @export
pcaP <- function(q, lambda, lower.tail = TRUE, log.p = FALSE) 
  return(pgamma(q, shape=lambda, scale=1, lower.tail = lower.tail, log.p = log.p))

#' Continuous approximation of Poisson -CAP- "inverse" probability function
#'
#' "Inverse" probability function for CAP.
#'
#' @param p probability
#' @param lambda mean and variance parameter
#' @param lower.tail if set to FALSE, returns inf of x such as Prob of X > x is > p
#' @param log.p if set to TRUE, returns the log-probability
#' @return for X following this distribution, returns inf Prob of X <= q
#' 
#' @export
qcaP <- function(p, lambda, lower.tail = TRUE, log.p = FALSE) 
  return(bisectionInv(function(q) pcaP(q, lambda, lower.tail = lower.tail, log.p = log.p), p))
         
#' Continuous approximation of Poisson -CAP- random generation function
#'
#' Random generation function for CAP.
#'
#' @param n number of samples
#' @param lambda mean and variance parameter
#' @return returns an array with n random samples with the CAP distribution
#' 
#' @importFrom stats rgamma
#' 
#' @export
rcaP <- function(n, lambda) return(rgamma(n, shape=lambda, scale=1))

#_________________________#
#   Negative-Binomial     #
#_________________________# 

#' Continuous approximation of Negative-Binomial -CANB- density function
#'
#' Density function for CANB.
#'
#' @param x value
#' @param mu mean parameter
#' @param delta dispersion parameter: var = mean+mean^2*delta
#' @return for X following this distribution, returns Prob of X = x
#' 
#' @importFrom stats dgamma
#' 
#' @export
dcaNB <- function(x, mu, delta, log=FALSE) 
  return(dgamma(x, shape=mu/(1+delta*mu), scale=1+mu*delta, log=log))

#' Continuous approximation of Negative-Binomial -CANB- cumulative probability function
#'
#' Cumulative probability function for CANB.
#'
#' @param q quantile
#' @param mu mean parameter
#' @param delta dispersion parameter: var = mean+mean^2*delta
#' @param lower.tail if set to FALSE, returns Prob of X > q
#' @param log.p if set to TRUE, returns the log-probability
#' @return for X following this distribution, returns Prob of X <= q
#' 
#' @importFrom stats GammaDist
#' 
#' @export
pcaNB <- function(q, mu, delta, lower.tail=TRUE, log.p=FALSE) {
  ## Returns poor NaN results due to implementation? https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=8528
  return(pgamma(q, shape=mu/(1+delta*mu), scale=1+mu*delta, lower.tail=lower.tail, log.p=log.p))
}

#' Continuous approximation of Negative-Binomial -CANB- "inverse" probability function
#'
#' "Inverse" probability function for CANB.
#'
#' @param p probability
#' @param mu mean parameter
#' @param delta dispersion parameter: var = mean+mean^2*delta
#' @param lower.tail if set to FALSE, returns inf of x such as Prob of X > x is > p
#' @param log.p if set to TRUE, returns the log-probability
#' @return for X following this distribution, returns inf Prob of X <= q
#' 
#' @export
qcaNB <- function(p, mu, delta, lower.tail=TRUE, log.p=FALSE) 
  return(bisectionInv(function(q) pcaNB(q, mu, delta, lower.tail=lower.tail, log.p=log.p), p))

#' Continuous approximation of Negative-Binomial -CANB- random generation function
#'
#' Random generation function for CANB.
#'
#' @param n number of samples
#' @param mu mean parameter
#' @param delta dispersion parameter: var = mean+mean^2*delta
#' @return returns an array with n random samples with the CANB distribution
#' 
#' @importFrom stats rgamma
#' 
#' @export
rcaNB <- function(n, mu, delta) 
  return(rgamma(n, shape=mu/(1+delta*mu), scale=1+mu*delta))

#___________#
#   CPNB    #
#___________# 

# Type: density probability distribution function (DDF) // Proof: linearity of the integral operator ##
#' CPNB density function
#'
#' Density function for CPNB.
#'
#' @param x value
#' @param pdr dropout probability
#' @param lambda_0 mean and variance for Poisson-like component
#' @param mu_i mean for Negative Binomial-like component
#' @param delta overdispersion for Negative Binomial-like component
#' @return for X following this distribution, returns Prob of X = x
#' 
#' @export
dCPNB <- function(x, pdr, lambda_0, mu_i, delta, log=FALSE)
  return(retLog(pdr*dcaP(x, lambda=lambda_0)+(1-pdr)*dcaNB(x, mu=mu_i, delta=delta), log))

# Type: cumulative probability distribution function (CDF) // Proof: definition of the mixture       ##
#' CPNB cumulative probability function
#'
#' Cumulative probability function for CPNB.
#'
#' @param q quantile
#' @param pdr dropout probability
#' @param lambda_0 mean and variance for Poisson-like component
#' @param mu_i mean for Negative Binomial-like component
#' @param delta overdispersion for Negative Binomial-like component
#' @param lower.tail if set to FALSE, returns Prob of X > q
#' @param log.p if set to TRUE, returns the log-probability
#' @return for X following this distribution, returns Prob of X <= q
#' 
#' @export
pCPNB <- function(q, pdr, lambda_0, mu_i, delta, lower.tail=TRUE, log.p=FALSE) {
  return(retLog(pdr*pcaP(q, lambda=lambda_0, lower.tail=lower.tail)+(1-pdr)*pcaNB(q, mu=mu_i, delta=delta, 
                                                                                  lower.tail=lower.tail), log.p))
}

# Type: quantile function // Returns x in CDF(mixture)^(-1)(p) = {x, CDF(mixture)(x) = p} computed by bisection ##
# (with an error of eps=1e-7, see docs for @bisectionInv)
#' CPNB "inverse" probability function
#'
#' "Inverse" probability function for CPNB.
#'
#' @param p probability
#' @param pdr dropout probability
#' @param lambda_0 mean and variance for Poisson-like component
#' @param mu_i mean for Negative Binomial-like component
#' @param delta overdispersion for Negative Binomial-like component
#' @param lower.tail if set to FALSE, returns inf of x such as Prob of X > x is > p
#' @param log.p if set to TRUE, returns the log-probability
#' @return for X following this distribution, returns inf Prob of X <= q
#' 
#' @export
qCPNB <- function(p, pdr, lambda_0, mu_i, delta, lower.tail=TRUE, log.p=FALSE)
  return(bisectionInv(function(q) return(pCPNB(ifelse(lower.tail, q, 1-q), pdr, lambda_0, mu_i, delta)), 
                      ifelse(log.p, exp(p), p)))

# Type: random generator // Samples a value u in [0,1] uniformly and returns x in CDF(mixture)^(-1)(u) ##
#' CPNB random generation function
#'
#' Random generation function for CPNB.
#'
#' @param n number of samples
#' @param pdr dropout probability
#' @param lambda_0 mean and variance for Poisson-like component
#' @param mu_i mean for Negative Binomial-like component
#' @param delta overdispersion for Negative Binomial-like component
#' @return returns an array with n random samples with the CPNB distribution
#' 
#' @export
rCPNB <- function(n, pdr, lambda_0, mu_i, delta)
  return(sapply(runif(n, 0, 1), function(u) qCPNB(u, pdr, lambda_0, mu_i, delta)))

##################################
# Marginal distribution creation #
##################################

#___________#
#   Margins #
#___________# 

## for x < 0, Prob(-X = x) = Prob(X = -x) ##
#' Mixture density function
#'
#' Density function for mixture.
#'
#' @param x value
#' @param pattern_ji pattern for gene i, cell j -0 or 1
#' @param pdr dropout probability
#' @param lambda_0 mean and variance for Poisson-like component
#' @param mu_i mean for Negative Binomial-like component
#' @param delta overdispersion for Negative Binomial-like component
#' @return for X following this distribution, returns Prob of X = x
#' 
#' @export
dmixture <- function(x, pattern_ji, pdr, lambda_0, mu_i, delta, log=FALSE) {
  return(dCPNB(x, pdr, lambda_0, mu_i, delta, log))
}

## for x < 0, Prob(-X <= x) = Prob(X >= -x) = 1 - Prob(X < -x) + symmetry ##
#' Mixture cumulative probability function
#'
#' Cumulative probability function for mixture.
#'
#' @param q quantile
#' @param pattern_ji pattern for gene i, cell j -0 or 1
#' @param pdr dropout probability
#' @param lambda_0 mean and variance for Poisson-like component
#' @param mu_i mean for Negative Binomial-like component
#' @param delta overdispersion for Negative Binomial-like component
#' @param lower.tail if set to FALSE, returns Prob of X > q
#' @param log.p if set to TRUE, returns the log-probability
#' @return for X following this distribution, returns Prob of X <= q
#' 
#' @export
pmixture <- function(q, pattern_ji, pdr, lambda_0, mu_i, delta, lower.tail=TRUE, log.p=FALSE) {
  if (pattern_ji == 0) return(pCPNB(q, pdr, lambda_0, mu_i, delta, lower.tail=lower.tail, log.p=log.p))
  return(1-pCPNB(q, pdr, lambda_0, mu_i, delta, lower.tail=lower.tail, log.p=log.p))
}

## inf {x in R+ | Prob(-X <= -x) >= p} = inf {x in R+ | Prob(X >= x) >= p} ##
#' Mixture "inverse" probability function
#'
#' "Inverse" probability function for mixture.
#'
#' @param p probability
#' @param pattern_ji pattern for gene i, cell j -0 or 1
#' @param pdr dropout probability
#' @param lambda_0 mean and variance for Poisson-like component
#' @param mu_i mean for Negative Binomial-like component
#' @param delta overdispersion for Negative Binomial-like component
#' @param lower.tail if set to FALSE, returns inf of x such as Prob of X > x is > p
#' @param log.p if set to TRUE, returns the log-probability
#' @return for X following this distribution, returns inf Prob of X <= q
#' 
#' @export
qmixture <- function(p, pattern_ji, pdr, lambda_0, mu_i, delta, lower.tail=TRUE, log.p=FALSE) {
  if (pattern_ji == 1) return(-qCPNB(p, pdr, lambda_0, mu_i, delta, lower.tail=!lower.tail, log.p=log.p))
  return(qCPNB(p, pdr, lambda_0, mu_i, delta, lower.tail=lower.tail, log.p=log.p))
}

#' Mixture random generation function
#'
#' Random generation function for mixture.
#'
#' @param n number of samples
#' @param pattern_ji pattern for gene i, cell j -0 or 1
#' @param pdr dropout probability
#' @param lambda_0 mean and variance for Poisson-like component
#' @param mu_i mean for Negative Binomial-like component
#' @param delta overdispersion for Negative Binomial-like component
#' @return returns an array with n random samples with the CPNB distribution
#' 
#' @export
rmixture <- function(n, pattern_ji, pdr, lambda_0, mu_i, delta) {
  if (pattern_ji == 1) return(-rCPNB(n, pdr, lambda_0, mu_i, delta))
  return(rCPNB(n, pdr, lambda_0, mu_i, delta))
}

######################################
# Distribution parameters estimation #
######################################

## Helper ##
#' Return X^T * X -vector
#'
#' Returns X^T * X -vector.
#'
#' @param x vector
#' @return x^t * x
#' 
#' @export
sumSquares <- function(x) return(x %*% x)

## IFM method: estimate the best values for margin parameters, and then fit the copula parameters using these estimators ##

#' Estimate margin parameters
#'
#' Estimates margin distribution parameters.
#'
#' @param p number of features in assay in object @sca
#' @param m number of conditions in assai in object @sca
#' @param P binary matrix of gene expression patterns
#' @param sca SingleCellAssay containing the data
#' @param is_mild vector indicating for each cell for each gene if this gene is mildly-expressed 
#' in the given cell
#' @return mean, mixture rate and dispersion parameters for each marginal distribution
#' 
#' @importFrom SummarizedExperiment assay
#' @importFrom stats median
#' @importFrom stats optim
#' @importFrom stats predict
#' 
#' @export
computeMarginParam <- function(p, m, P, modelDropout, sca, is_mild) {
  
  ## De-log the gene expression matrix ##
  M <- exp(assay(sca))
  
  ##_______________________________________________________________________##
  ## Size factor estimation (adapted from Anders et al. 2013)              ##
  ##_______________________________________________________________________##
  
  ## Size factors                                                          ##
  ## Pseudo count of 1 to avoid null issues                                ##
  ## Computation of the geometric mean of the rows                         ##
  pp <- apply(M, 1, function(x) prod((x+1)**(1/m)))
  sizeFact <- apply(M, 2, function(x) median(x/pp))
  
  ## Continuous approximation of Poisson-Negative Binomial (CAPNB) mixture parameters = ##
  ## {meanNB_cell,gene; overdispersionNB_cell; probability of dropout_gene}             ##
  
  ##_______________________________________________________________________##
  ## Mean estimation (adapted from Yu et al. 2013)                         ##
  ##_______________________________________________________________________##
  
  ## access to mu_ik = mean with mu[i], i gene                             ##
  mu <- log(apply(M, 1, function(x) sum(x/(m*sizeFact))))
  
  ##_______________________________________________________________________##
  ## Delta (dispersion parameter) estimation (adapted from Yu et al. 2013) ##
  ##_______________________________________________________________________##
  
  ##     var = mean + mean^2 * delta                                       ## 
  ## The difference with Yu et al. is that we have access only to the      ##
  ## normalized gene expression values                                     ##
  ## Naive estimation of variance                                          ##
  v <- sapply(1:p, function(i) sumSquares(M[i, ]/sizeFact-mu[i])/(m-1))
  ## Unskrunk dispersion parameter                                         ##
  est_delta <- sapply(1:p, function(i) max(0, (m*v[i]-mu[i]*sum(1/sizeFact))/(mu[i]**2*sum(1/sizeFact))))
  weight <- function(zeta) return((sumSquares(est_delta-mean(est_delta))/sumSquares(est_delta-zeta))*(p-2)/(p-1))
  ## slope of function zeta -> 1/sum_i(est_delta_i-zeta)^2 > 0               ##
  slope <- function(zeta) return(-2*sum(est_delta-zeta)/(sumSquares(est_delta-zeta))**2)
  ## Gradient of slope respect to zeta                                       ##
  #gslope <- function(zeta) return(-2*p*(-1)/(sumSquares(est_delta-zeta))**2
  #                                -(-2)*sum(est_delta-zeta)*2*sumSquares(est_delta-zeta)
  #                                 *(-2)*sum(est_delta-zeta)/(sumSquares(est_delta-zeta))**4)
  gslope <- function(zeta) return(2*p/(sumSquares(est_delta-zeta))**2
                                  -8*(sum(est_delta-zeta)**2)/(sumSquares(est_delta-zeta))**3)
  epsilon <- 0.05
  ## Get an unconstrained minimization optimization problem with log-barrier ##
  ##         0 > slope(zeta) > -epsilon and zeta > 0                         ##
  ## and a penalty on the complexity of the model (L2 norm)                  ##
  f <- function(zeta) return(zeta-log(epsilon+slope(zeta))-log(-slope(zeta)*zeta)+(zeta)**2)
  ## Gradient of f (respect to the sole variable zeta)                       ##
  g <- function(zeta) return(1-gslope(zeta)/(epsilon+slope(zeta))+(slope(zeta)+zeta*gslope(zeta))/(slope(zeta)*zeta)+2*zeta)
  zeta <- min(est_delta)
  fit <- optim(par=c(zeta), fn=f, method="BFGS", gr=g)
  zeta <- fit$par
  ## Shrinked dispersion parameter                                          ##
  delta <- (1-weight(zeta))*(est_delta-zeta)+zeta
  pdrArray <- sapply(1:p, function(i) predict(modelDropout, 
                                              data.frame(Value=sum(sapply(1:m, 
                                                                          function(j) assay(sca)[i, j]/sum(assay(sca)[, j])))), 
                                                         type="response"))
  indices_high <- apply(is_mild, 2, function(x) which(x == "high"))
  indices_mild <- lapply(indices_high, function(a) (1:p)[!(1:p %in% a)])
  
  return(list(mu, delta, pdrArray, indices_high, indices_mild))
}

################################
# Copula distribution creation #
################################

#' Create custom bimodal Gaussian copula
#'
#' Creates a bimodal Gaussian copula for cell j.
#'
#' @param p the number of features
#' @param j index of cell
#' @param M trimmed log-normalized gene expression matrix
#' @param mixParam list of mixture rates for bimodal copula
#' @param muCop list of mean parameters for bimodal copula
#' @param rhoCop list of covariance parameters for bimodal copula
#' @param P pattern matrix
#' @param pdrArray array of dropout event probabilities for each gene
#' @param lambda_0 parameter for dropout component in gene expression
#' @param mu mean for regular expressiveness component in gene expression
#' @param delta dispersion parameter for gene expression
#' @return a 2-Gaussian copula mixture
#' 
#' @importFrom stats qnorm
#' @imporFrom stats pnorm
#' @importFrom mixtools rmvnorm 
#' @importFrom mnormt pmnorm
#' @importFrom mnormt dmnorm 
#' @importFrom mnormt rmnorm
#' 
#' @export
bimodalNormalCopula <- function(p, j, M, mixParam, muCop, rhoCop, P, pdrArray, lambda_0, mu, delta) {
  pattern_j = P[, j]
  kk <- length(muCop)
  invCop <- lapply(1:kk, function(k) invCustom(rhoCop[[k]]))
  rnormmv <- function(rho, mu) return(pnorm(rmnorm(n=n, mean=mu, varcov=rho)))

  ##______________________________________________##
  ## Proof for Gaussian copulas in Song (2000)    ##
  ##______________________________________________##
  ## x is a NON-TRANSFORMED VECTOR                ##
  dCopula <- function(x, log=F) {
    q <- qnorm(sapply(1:p, function(i) pmixture(x[i], pattern_j[i], pdrArray[i], lambda_0, mu[i], delta[i])))
    q[is.infinite(q) & q < 0] <- -5
    q[is.infinite(q) & q > 0] <- 5
    #retvals <- sum(sapply(1:kk, function(k) mixParam[k]*dmnorm(q, mean=muCop[[k]], varcov=rhoCop[[k]])))
    retvals <- sum(sapply(1:kk, function(k) mixParam[k]*dmvnormCustom(q, muCop[[k]], rhoCop[[k]], invrho = invCop[[k]])))
    return(ifelse(log, log(retvals), retvals))
  }
  
  ##_________________________________________##
  ## Copula definition                       ##
  ##_________________________________________##
  ## x is a NON-TRANSFORMED VECTOR           ##
  pCopula <- function(x, log.p=F, lower.tail=T) {
    q <- sapply(1:p, function(i) pmixture(x[i], pattern_j[i], pdrArray[i], lambda_0, mu[i], delta[i]))
    print(q)
    q <- qnorm(q)
    print(q)
    q[is.infinite(q) & q < 0] <- -5
    q[is.infinite(q) & q > 0] <- 5
    retvals <- sapply(1:kk, function(k) mixParam[k]*pmnorm(x=as.numeric(q), mean=muCop[[k]], varcov=rhoCop[[k]]))
    print(retvals)
    retvals <- sum(retvals)
    print(retvals)
    if (log.p) retvals <- log(retvals)
    if (!lower.tail) retvals <- 1-retvals
    return(retvals)
  }
  
  ##___________________##
  ## Random generation ##
  ##___________________##
  rCopula <- function(n) {
    m <- sample(1:kk, prob=mixParam, size=n, replace=T)
    rvectors <- base::t(matrix(unlist(lapply(m, function(k) rmnorm(n=1, mean=muCop[[k]], varcov=rhoCop[[k]]))), nrow=n))
    pvectors <- pnorm(rvectors)
    p <- nrow(pvectors)
    uvectors <- NULL
    for (j in 1:n) {
        uvectors <- cbind(uvectors, sapply(1:p, function(i) {print(i); qmixture(pvectors[i, j], pattern_j[i], pdrArray[i], lambda_0, 
                                     mu[i], delta[i])}))
    }
    #uvectors <- lapply(1:n, function(j) sapply(1:p, 
    #                                           function(i) qmixture(pvectors[i, j], pattern_j[i], pdrArray[i], lambda_0, mu[i], delta[i])))
    #uvectors <- matrix(unlist(uvectors), ncol=n)
    return(uvectors)
  }
  
  return(list(fct=list(dCopula=dCopula, 
              pCopula=pCopula, 
              rCopula=rCopula),
              ## For debugging purposes ##
              param=list(mixParam=mixParam, rhoCop=rhoCop, muCop=muCop, pattern_j=pattern_j, 
                         pdrArray=pdrArray, lambda_0=lambda_0, mu=mu, delta=delta)))
  }

dBCopula <- function(x, copula, log=F) return(copula$fct$dCopula(x, log))
pBCopula <- function(q, copula, log.p=F, lower.tail=T) return(copula$fct$pCopula(q, log.p, lower.tail))
rBCopula <- function(n, copula) return(copula$fct$rCopula(n))
