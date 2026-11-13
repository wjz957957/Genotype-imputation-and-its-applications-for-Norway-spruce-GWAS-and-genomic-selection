#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
WORK <- if(length(args))normalizePath(args[1],mustWork=TRUE) else getwd()
ROOT <- file.path(WORK,'10_outer_test/09_outer_prediction_crosschip_r2_0p8_all3positive')
DIR <- file.path(ROOT,'comparison')
dir.create(DIR,recursive=TRUE,showWarnings=FALSE)
metrics <- c('Pearson','Spearman','RMSE','MAE','Bias')
for(trait in c('BB','H6')) {
  rows <- lapply(sprintf('%02d',1:10),function(outer) {
    d <- file.path(ROOT,paste0('outer_',outer),trait)
    if (!file.exists(file.path(d,'complete.flag')) || !file.exists(file.path(d,'outer_test_predictions.tsv'))) stop('Missing completion/predictions: ',d)
    m <- read.table(file.path(d,'outer_metrics.tsv'),header=TRUE,sep='\t',check.names=FALSE)
    if(nrow(m)!=1L || sprintf("%02d",as.integer(m$Outer))!=outer || m$Trait!=trait) stop('Invalid output: ',d)
    if(m$Best_Extra==0 && any(abs(unlist(m[paste0('Delta_',metrics)]))>1e-12)) stop('Zero QTN produced nonzero deltas: ',d)
    m
  })
  dat <- do.call(rbind,rows)
  rownames(dat)<-NULL
  per_fold <- dat[,c('Trait','Outer','Best_Chip','N_chip_selected','Best_Extra','N_total_selected',
                     unlist(lapply(metrics,function(x)c(paste0(x,'_CHIP'),paste0(x,'_COMBINED'),paste0('Delta_',x)))))]
  write.table(per_fold,file.path(DIR,paste0(trait,'_CHIP_vs_all3positive_per_fold.tsv')),
              sep='\t',quote=FALSE,row.names=FALSE)
  summary_rows <- lapply(metrics,function(name) {
    x <- dat[[paste0(name,'_CHIP')]]; y <- dat[[paste0(name,'_COMBINED')]]; z <- y-x
    if(length(z)!=10L || any(!is.finite(z))) stop('Invalid 10 folds')
    if(all(z==0)) { t_p <- 1; wil_p <- 1; ci <- c(0,0) }
    else {
      test <- t.test(z,mu=0)
      ci <- test$conf.int; t_p <- test$p.value
      wil_p <- tryCatch(suppressWarnings(wilcox.test(z,mu=0,exact=FALSE)$p.value),error=function(e)NA_real_)
    }
    better <- if(name %in% c('RMSE','MAE')) sum(z<0) else if(name=='Bias') NA_integer_ else sum(z>0)
    worse <- if(name %in% c('RMSE','MAE')) sum(z>0) else if(name=='Bias') NA_integer_ else sum(z<0)
    data.frame(Trait=trait,Metric=name,N_folds=10,
      CHIP_mean=mean(x),CHIP_SD=sd(x),COMBINED_mean=mean(y),COMBINED_SD=sd(y),
      Mean_Delta=mean(z),SD_Delta=sd(z),CI95_low=ci[1],CI95_high=ci[2],
      Paired_t_p=t_p,Paired_Wilcoxon_p=wil_p,
      Better_folds=better,Worse_folds=worse,Equal_folds=sum(z==0))
  })
  sumdat <- do.call(rbind,summary_rows)
  write.table(sumdat,file.path(DIR,paste0(trait,'_CHIP_vs_all3positive_summary.tsv')),
              sep='\t',quote=FALSE,row.names=FALSE)
  cat('\n===== ',trait,' 10-fold summary =====\n',sep='')
  print(sumdat[,c('Metric','CHIP_mean','COMBINED_mean','Mean_Delta','Paired_t_p','Better_folds','Worse_folds','Equal_folds')],row.names=FALSE,digits=7)
}
cat('\nVERIFIED ALL 20 RR-BLUP OUTER RESULTS\nComparison files: ',DIR,'\n',sep='')
