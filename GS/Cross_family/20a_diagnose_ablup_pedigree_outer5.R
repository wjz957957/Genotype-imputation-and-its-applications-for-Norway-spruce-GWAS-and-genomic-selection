#!/usr/bin/env Rscript
# Diagnostic only: compute the pedigree A and its 5 held-out cross blocks.
# No phenotype is loaded and no model is fitted.
args <- commandArgs(trailingOnly=TRUE)
work <- if (length(args)>=1L) args[1L] else "."
work <- normalizePath(work, mustWork=TRUE)
pedfile <- file.path(work, "pedigree.csv")
assignfile <- file.path(work,"Cross_family/config/progeny_parent_connected_outer_fold_assignment.tsv")
stopifnot(file.exists(pedfile), file.exists(assignfile))
p <- read.csv(pedfile, colClasses="character", check.names=FALSE, fileEncoding="UTF-8-BOM")
x <- read.delim(assignfile, check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(all(c("ID","Sire","Dam") %in% names(p)),
          all(c("Sample.Filename","Outer_test_fold","Parent_connected_group") %in% names(x)))
stopifnot(!anyDuplicated(p$ID), !anyDuplicated(x$Sample.Filename), all(x$Sample.Filename %in% p$ID))
n <- nrow(p)
parent_index <- function(v) {
  z <- is.na(v) | v %in% c("0", "")
  i <- match(v,p$ID); i[z] <- 0L
  if (anyNA(i)) stop("Pedigree contains parent IDs not found in pedigree")
  i
}
s <- parent_index(p$Sire); d <- parent_index(p$Dam)
if (any(s>=seq_len(n) | d>=seq_len(n))) stop("Pedigree parent order invalid")
A <- matrix(0, n, n)
for (i in seq_len(n)) {
  if (i>1L) {
    j <- seq_len(i-1L)
    v <- numeric(i-1L)
    if (s[i]>0L) v <- v + A[s[i],j]/2
    if (d[i]>0L) v <- v + A[d[i],j]/2
    A[i,j] <- v; A[j,i] <- v
  }
  A[i,i] <- 1 + if(s[i]>0L && d[i]>0L) A[s[i],d[i]]/2 else 0
}
rownames(A) <- colnames(A) <- p$ID
fold <- as.integer(x$Outer_test_fold)
for (o in 1:5) {
  test <- x$Sample.Filename[fold==o]
  train <- x$Sample.Filename[fold!=o]
  C <- A[test, train, drop=FALSE]
  cat(sprintf("outer_%02d test=%d train=%d max_abs_A_test_train=%.15g n_nonzero_entries=%d\n",
    o, length(test), length(train), max(abs(C)), sum(abs(C)>1e-12)))
}
cat("If every max_abs_A_test_train is zero, ABLUP is intercept-only in all five folds.\n")
