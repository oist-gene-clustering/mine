load("tintori.Rdata")
source("mine.R")
library(scater)

generateCopulas <- function(cellitv=1:30, geneitv=1:3000) {
     s <- system.time(cl <- mine(counts(tintori)[geneitv, cellitv]))
     copulas <- cl[[1]]
     matr <- cl[[2]]
     save(s, copulas, matr, file="tintori-utils.Rdata")
}

plotGEvalues <- function(cellidx) {
      load("tintori-utils.Rdata")
      hist(matr[, cellidx], main=paste0("Histogram of gene expression values in cell ", strsplit(colnames(matr)[cellidx], "_")[[1]][2], " embryo ", strsplit(colnames(matr)[cellidx], "_")[[1]][1]), col="grey", xlab = "Gene expression values")
} 

plotCopula <- function(cellidx) {
      load("tintori-utils.Rdata")
      tt <- system.time(rr <- rBCopula(1, copulas[[cellidx]][1]))
      save(rr, tt, s, copulas, matr, file="tintori-copula.Rdata")
      hist(rr, main=paste0("Histogram of gene expression values in cell ", strsplit(colnames(matr)[cellidx], "_")[[1]][2], " embryo ", strsplit(colnames(matr)[cellidx], "_")[[1]][1], " from copula"), col="grey", xlab = "Gene expression values")
}

plotCDF <- function(itv, cellidx) {
      load("tintori-utils.Rdata")
      pp <- NULL
      yy <- system.time(
      for (j in itv) {
           pp <- cbind(pp, pBCopula(sapply(copulas[[cellidx]][2]$param$pattern_j, function(i) ifelse(i == 0, j, -j)), copulas[[cellidx]][1]))
      })
      plot(itv, pp)
      save(pp, yy, file="tintori-function.Rdata")
}

plotPDF <- function(itv, cellidx) {
      load("tintori-utils.Rdata")
      pp <- NULL
      yy <- system.time(
      for (j in itv) {
           pp <- cbind(pp, dBCopula(sapply(copulas[[cellidx]][2]$param$pattern_j, function(i) ifelse(i == 0, j, -j)), copulas[[cellidx]][1]))
      })
      plot(itv, pp)
      save(pp, yy, file="tintori-funct.Rdata")
}
