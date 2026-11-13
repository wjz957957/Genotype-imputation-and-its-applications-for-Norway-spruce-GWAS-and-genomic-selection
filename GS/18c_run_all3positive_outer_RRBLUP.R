#!/usr/bin/env Rscript
# Fold-paired RR-BLUP on OLD Dardel server; uses Step18b inputs.
# A baseline model is fitted on the identical CHIP matrix for every fold.
# Only the selected QTN count differs from Step17; no test leakage in fitting.
suppressPackageStartupMessages(library(rrBLUP))
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L) stop('Usage: Rscript 18c_run_all3positive_outer_RRBLUP.R WORK OUTER TRAIT')
WORK <- normalizePath(args[1], mustWork=TRUE)
OUTER <- args[2]; TRAIT <- args[3]
if (!(OUTER %in% sprintf('%02d',1:10)) || !(TRAIT %in% c('BB','H6'))) stop('Bad fold/trait')
BASE <- file.path(WORK,'10_outer_test','08_final_genotypes_crosschip_r2_0p8_all3positive',paste0('outer_',OUTER),TRAIT)
OUT <- file.path(WORK,'10_outer_test','09_outer_prediction_crosschip_r2_0p8_all3positive',paste0('outer_',OUTER),TRAIT)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

read_ids <- function(path) {
  if (!file.exists(path)) stop('Missing CV IDs: ',path)
  ids <- scan(path,what='',quiet=TRUE)
  if (anyDuplicated(ids) || length(ids)<1L) stop('Invalid CV IDs')
  ids
}
read_raw <- function(path) {
  if (!file.exists(path) || file.info(path)$size<=0) stop('Missing RAW: ',path)
  if (requireNamespace('data.table',quietly=TRUE)) {
    d <- data.table::fread(path,header=TRUE,check.names=FALSE,data.table=FALSE,na.strings='NA')
  } else {
    d <- read.table(path,header=TRUE,check.names=FALSE,stringsAsFactors=FALSE,
                    comment.char='',quote='',na.strings='NA')
  }
  if (ncol(d)<7L || !identical(names(d)[1:2],c('FID','IID'))) stop('Invalid RAW columns: ',path)
  ids <- as.character(d[['IID']])
  if (anyDuplicated(ids) || any(as.character(d[['FID']])!=ids)) stop('Invalid RAW IDs: ',path)
  header <- colnames(d)[-(1:6)]
  if (anyDuplicated(header)) stop('Duplicate genotype columns: ',path)
  snps <- sub('_[^_]+$','',header)
  if (anyDuplicated(snps)) stop('Duplicate SNP IDs in RAW: ',path)
  G <- data.matrix(d[,-(1:6),drop=FALSE]); storage.mode(G) <- 'numeric'
  if (anyNA(G) || any(!is.finite(G))) stop('NA/nonfinite genotype: ',path)
  if (any(!(G %in% c(0,1,2,-9)))) stop('Non-dosage value: ',path)
  list(ids=ids,header=header,snps=snps,G=G)
}

manifest_file <- file.path(BASE,'final_genotype_manifest.tsv')
if (!file.exists(manifest_file) || !file.exists(file.path(BASE,'complete.flag'))) stop('Incomplete Step18b inputs')
manifest <- read.table(manifest_file,header=TRUE,sep='\t',check.names=FALSE,stringsAsFactors=FALSE)
if (nrow(manifest)!=1L || sprintf("%02d",as.integer(manifest$Outer))!=OUTER || as.character(manifest$Trait)!=TRAIT) stop('Bad manifest')
N_CHIP <- as.integer(manifest$N_chip); N_EXTRA <- as.integer(manifest$Best_Extra)
TRAIN_IDS <- read_ids(file.path(WORK,'02_CV',paste0('outer_',OUTER),'outer_train.samples'))
TEST_IDS <- read_ids(file.path(WORK,'02_CV',paste0('outer_',OUTER),'outer_test.samples'))
if (length(TRAIN_IDS)!=manifest$N_train || length(TEST_IDS)!=manifest$N_test) stop('CV sample count mismatch')
if (length(intersect(TRAIN_IDS,TEST_IDS))!=0L) stop('CV train/test overlap')

c_train <- read_raw(file.path(BASE,'chip_train.raw'))
c_test <- read_raw(file.path(BASE,'chip_test.raw'))
if (!identical(c_train$ids,TRAIN_IDS) || !identical(c_test$ids,TEST_IDS)) stop('CHIP RAW sample order differs from CV, refusing to infer order')
if (!identical(c_train$header,c_test$header) || ncol(c_train$G)!=N_CHIP) stop('CHIP RAW train/test columns differ')

if (N_EXTRA>0L) {
  e_train <- read_raw(file.path(BASE,'extra_train.raw'))
  e_test <- read_raw(file.path(BASE,'extra_test.raw'))
  if (!identical(e_train$ids,TRAIN_IDS) || !identical(e_test$ids,TEST_IDS)) stop('Extra RAW sample order differs from CV')
  if (!identical(e_train$header,e_test$header) || ncol(e_train$G)!=N_EXTRA) stop('Extra RAW columns differ or wrong count')
  if (length(intersect(c_train$snps,e_train$snps))) stop('CHIP/QTN overlap')
} else {
  e_train <- NULL; e_test <- NULL
}

PHENO <- file.path(WORK,'00_input',if (TRAIT=='BB') 'Trait-Adj_Budburst.txt' else 'Trait-Adj_H6.txt')
p <- read.table(PHENO,header=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
if (ncol(p)!=3L || any(as.character(p[[1]])!=as.character(p[[2]]))) stop('Phenotype format mismatch')
pid <- as.character(p[[2]]); y <- as.numeric(p[[3]])
if (anyDuplicated(pid) || any(!is.finite(y))) stop('Invalid phenotype IDs/values')
ph <- setNames(y,pid)
y_train <- unname(ph[TRAIN_IDS]); y_test <- unname(ph[TEST_IDS])
if (any(!is.finite(y_train)) || any(!is.finite(y_test))) stop('Missing phenotype in fold')
ph_sd <- sd(y); ph_mean <- mean(y)
if (!is.finite(ph_sd) || ph_sd<=0) stop('Invalid full phenotype SD')

fit_predict <- function(G_train,G_test) {
  # Exactly the rrBLUP::mixed.solve model used in prior new-server RRBLUP runs.
  fit <- rrBLUP::mixed.solve(y=y_train,Z=G_train,SE=FALSE,return.Hinv=FALSE)
  pred <- as.numeric(G_test %*% fit$u + as.numeric(fit$beta))
  if (length(pred)!=length(y_test) || any(!is.finite(pred))) stop('Invalid predictions')
  list(pred=pred, Vu=as.numeric(fit$Vu), Ve=as.numeric(fit$Ve))
}
metric <- function(pred) {
  z_test <- (y_test-ph_mean)/ph_sd
  z_pred <- (pred-ph_mean)/ph_sd
  c(Pearson=cor(y_test,pred,method='pearson'),
    Spearman=cor(y_test,pred,method='spearman'),
    RMSE=sqrt(mean((z_pred-z_test)^2)), MAE=mean(abs(z_pred-z_test)),
    Bias=mean(z_pred-z_test),RMSE_raw=sqrt(mean((pred-y_test)^2)),
    MAE_raw=mean(abs(pred-y_test)), Bias_raw=mean(pred-y_test))
}
cat(sprintf('[RUN] outer=%s %s chip=%s extra=%s, train=%d test=%d\n',OUTER,TRAIT,N_CHIP,N_EXTRA,length(TRAIN_IDS),length(TEST_IDS)))
baseline_fit <- fit_predict(c_train$G,c_test$G)
b <- metric(baseline_fit$pred)
# An independent parity check against the reported new-server CHIP baseline for outer01.
expected01 <- c(BB=0.742269974193988, H6=0.488048864182112)
if (OUTER=='01' && abs(b['Pearson']-expected01[[TRAIT]])>5e-5) {
  stop(sprintf('CHIP baseline Pearson mismatch: old-server=%.12f, prior-new-server=%.12f. DO NOT compare across pipelines until investigated.', b['Pearson'],expected01[[TRAIT]]))
}
if (N_EXTRA>0L) {
  hyb_fit <- fit_predict(cbind(c_train$G,e_train$G),cbind(c_test$G,e_test$G))
} else {
  hyb_fit <- baseline_fit  # 0 QTN is CHIP-only, no redundant model fit
}
h <- metric(hyb_fit$pred)

pred <- data.frame(Outer=OUTER,Trait=TRAIT,IID=TEST_IDS,Observed_raw=y_test,
                   RRBLUP_CHIP_raw=baseline_fit$pred,RRBLUP_all3positive_raw=hyb_fit$pred)
write.table(pred,file.path(OUT,'outer_test_predictions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
met <- data.frame(Outer=OUTER,Trait=TRAIT,Best_Chip=manifest$Best_Chip,
    N_chip_selected=N_CHIP,Best_Extra=N_EXTRA,N_total_selected=N_CHIP+N_EXTRA,
    N_train=length(TRAIN_IDS),N_test=length(TEST_IDS),Phenotype_SD_full=ph_sd,
    Pearson_CHIP=b['Pearson'],Pearson_COMBINED=h['Pearson'],
    Spearman_CHIP=b['Spearman'],Spearman_COMBINED=h['Spearman'],
    RMSE_CHIP=b['RMSE'],RMSE_COMBINED=h['RMSE'],
    MAE_CHIP=b['MAE'],MAE_COMBINED=h['MAE'],
    Bias_CHIP=b['Bias'],Bias_COMBINED=h['Bias'],
    Delta_Pearson=h['Pearson']-b['Pearson'],Delta_Spearman=h['Spearman']-b['Spearman'],
    Delta_RMSE=h['RMSE']-b['RMSE'],Delta_MAE=h['MAE']-b['MAE'],Delta_Bias=h['Bias']-b['Bias'],
    RMSE_raw_CHIP=b['RMSE_raw'],RMSE_raw_COMBINED=h['RMSE_raw'],
    stringsAsFactors=FALSE)
write.table(met,file.path(OUT,'outer_metrics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(paste('DONE',OUTER,TRAIT,N_EXTRA),file.path(OUT,'complete.flag'))
cat(sprintf('[PASS] %s outer_%s: extra=%d CHIP Pearson=%.8f combined Pearson=%.8f delta=%+.8f\n',
            TRAIT,OUTER,N_EXTRA,b['Pearson'],h['Pearson'],h['Pearson']-b['Pearson']))
