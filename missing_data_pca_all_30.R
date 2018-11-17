
# This file performs the competing PCA approaches on the missing value data sets

#### Edit the below file_root to be the location of the saved folder:
file_root = "C:/........./paper_r_examples/"

source(file = paste0(file_root, "glrm_functions.R"))

# Load Libraries:
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("pcaMethods")
library(pcaMethods)

MCAR_offset_df = data.frame(baseline = rep(0,30), bpca = rep(0,30), ppca = rep(0,30), nipals = rep(0,30))
MCAR_MSE_df    = data.frame(baseline = rep(0,30), bpca = rep(0,30), ppca = rep(0,30), nipals = rep(0,30))
MAR_offset_df  = data.frame(baseline = rep(0,30), bpca = rep(0,30), ppca = rep(0,30), nipals = rep(0,30))
MAR_MSE_df     = data.frame(baseline = rep(0,30), bpca = rep(0,30), ppca = rep(0,30), nipals = rep(0,30))

set.seed(123)
for(dd in 1:30){
  #### MCAR ####
  # Load Data:
  MCAR       = read.csv(file=paste0(file_root,"missing_data/simulated_data/MCAR_data_",dd,".csv"))
  MCAR_truth = read.csv(file=paste0(file_root,"missing_data/simulated_data/MCAR_truth_",dd,".csv"))
  na_index = which(is.na(MCAR))
  bpca_fit = pca(MCAR, method = "bpca", center=TRUE, scale="uv", completeObs = TRUE,nPcs=4, seed=123)
  nipals_fit = pca(MCAR, method = "nipals", center=TRUE, scale="uv", completeObs = TRUE,nPcs=4, seed=123)
  ppca_fit = pca(MCAR, method = "ppca", center=TRUE, scale="uv", completeObs = TRUE,nPcs=4, seed=123)
  #summary(bpca_fit)
  #summary(nipals_fit)
  #summary(ppca_fit)
  ### Get the estimated complete observations
  bpca_est = completeObs(bpca_fit)
  nipals_est = completeObs(nipals_fit)
  ppca_est = completeObs(ppca_fit)
  MCAR_MSE_df$bpca[dd] = mean((bpca_est[na_index,1] - MCAR_truth[,1])^2)
  MCAR_offset_df$bpca[dd] = mean(bpca_est[,1])
  MCAR_MSE_df$nipals[dd] = mean((nipals_est[na_index,1] - MCAR_truth[,1])^2)
  MCAR_offset_df$nipals[dd] = mean(nipals_est[,1])
  MCAR_MSE_df$ppca[dd] = mean((ppca_est[na_index,1] - MCAR_truth[,1])^2)
  MCAR_offset_df$ppca[dd] =mean(ppca_est[,1])
  ## sample mean model:
  xbar = mean(MCAR[,1],na.rm = TRUE)
  MCAR_offset_df$baseline[dd] = xbar
  MCAR_MSE_df$baseline[dd] = mean((xbar - MCAR_truth[,1])^2)

  #### MAR ####
  # Load Data:
  MAR        = read.csv(file=paste0(file_root,"missing_data/simulated_data/MAR_data_",dd,".csv"))
  MAR_truth  = read.csv(file=paste0(file_root,"missing_data/simulated_data/MAR_truth_",dd,".csv"))
  na_index = which(is.na(MAR))
  bpca_fit = pca(MAR, method = "bpca", center=TRUE, scale="uv", completeObs = TRUE,nPcs=4, seed=123)
  nipals_fit = pca(MAR, method = "nipals", center=TRUE, scale="uv", completeObs = TRUE,nPcs=4, seed=123)
  ppca_fit = pca(MAR, method = "ppca", center=TRUE, scale="uv", completeObs = TRUE,nPcs=4, seed=123)
  #summary(bpca_fit)
  #summary(nipals_fit)
  #summary(ppca_fit)
  ### Get the estimated complete observations
  bpca_est = completeObs(bpca_fit)
  nipals_est = completeObs(nipals_fit)
  ppca_est = completeObs(ppca_fit)
  MAR_MSE_df$bpca[dd] = mean((bpca_est[na_index,1] - MAR_truth[,1])^2)
  MAR_offset_df$bpca[dd] = mean(bpca_est[,1])
  MAR_MSE_df$nipals[dd] = mean((nipals_est[na_index,1] - MAR_truth[,1])^2)
  MAR_offset_df$nipals[dd] = mean(nipals_est[,1])
  MAR_MSE_df$ppca[dd] = mean((ppca_est[na_index,1] - MAR_truth[,1])^2)
  MAR_offset_df$ppca[dd] = mean(ppca_est[,1])
  ## sample mean model:
  xbar = mean(MAR[,1],na.rm = TRUE)
  MAR_offset_df$baseline[dd] = xbar
  MAR_MSE_df$baseline[dd] = mean((xbar - MAR_truth[,1])^2)
}


#file_root2 = paste0(file_root,"missing_data/stored_results/"
#save(MCAR_MSE_df, file=paste0(file_root2,"MCAR_MSE_df.RData"))
#save(MCAR_offset_df, file=paste0(file_root2,"MCAR_offset_df.RData"))
#save(MAR_MSE_df, file=paste0(file_root2,"MAR_MSE_df.RData"))
#save(MAR_offset_df, file=paste0(file_root2,"MAR_offset_df.RData"))




