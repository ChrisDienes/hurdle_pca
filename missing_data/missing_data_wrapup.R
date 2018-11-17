
# This file wraps up the missing data analysis and reports objects used in final paper

#### Edit the below file_root to be the location of the saved folder:
file_root = "C:/........./paper_r_examples/"

file_root2 = pates0(file_root,"missing_data/stored_results/")

### load other PCA results:
load(file=paste0(file_root2,"MCAR_MSE_df.RData"))
load(file=paste0(file_root2,"MCAR_offset_df.RData"))
load(file=paste0(file_root2,"MAR_MSE_df.RData"))
load(file=paste0(file_root2,"MAR_offset_df.RData"))

### load hurdle results:
load(file=paste0(file_root2,"hurdle_MSE_MCAR.RData"))
load(file=paste0(file_root2,"hurdle_offset_MCAR.RData"))
load(file=paste0(file_root2,"hurdle_AUC_MCAR.RData"))
load(file=paste0(file_root2,"hurdle_MSE_MAR.RData"))
load(file=paste0(file_root2,"hurdle_offset_MAR.RData"))
load(file=paste0(file_root2,"hurdle_AUC_MAR.RData"))
load(file=paste0(file_root2,"hurdle_theta_MAR.RData"))
load(file=paste0(file_root2,"hurdle_dist_MAR.RData"))

### load penalty values:
load(file=paste0(file_root2,"penalty_summary_MCAR.RData"))
load(file=paste0(file_root2,"penalty_summary_MAR.RData"))

### combine results:
mcar_mse = cbind(MCAR_MSE_df, hurdle = hurdle_MSE_MCAR)
mcar_off = cbind(MCAR_offset_df, hurdle = hurdle_offset_MCAR)
mar_mse = cbind(MAR_MSE_df, hurdle = hurdle_MSE_MAR)
mar_off = cbind(MAR_offset_df, hurdle = hurdle_offset_MAR)

### explore results:

colMeans(mcar_mse)
apply(mcar_mse,2,sd)/sqrt(30)

colMeans((1 - mcar_off)^2)
apply((1 - mcar_off)^2,2,sd)/sqrt(30)

colMeans(mar_mse)
apply(mar_mse,2,sd)/sqrt(30)

colMeans((1 - mar_off)^2)
apply((1 - mar_off)^2,2,sd)/sqrt(30) 

mean(hurdle_AUC_MCAR)
sd(hurdle_AUC_MCAR)/sqrt(30)

mean(hurdle_AUC_MAR)
sd(hurdle_AUC_MAR)/sqrt(30)

hurdle_theta_MAR[1,]
colMeans(hurdle_theta_MAR)
apply(hurdle_theta_MAR,2,sd)/sqrt(30)

dist_ranks = t(apply(hurdle_dist_MAR,1,rank)-1)[,-2]
hold_tmp1 = 1:10
hold_tmp2 = 1:10
for(tx in 1:10){
  top_x = tx
  tmp = dist_ranks
  tmp[dist_ranks <= top_x] = 1
  tmp[dist_ranks > top_x] = 0
  hold_tmp1[tx] = 100*mean(rowSums(tmp[,2:3])>0)
  hold_tmp2[tx] = 100*mean(rowSums(tmp[,2:3])>1)
}
hold_tmp1[2]
hold_tmp1[3]

hurdle_dist_MAR[1,]
colMeans(hurdle_dist_MAR)
apply(hurdle_dist_MAR,2,sd)/sqrt(30)

improve_bpca_mcar = 100*(mcar_mse$bpca - mcar_mse$hurdle)/mcar_mse$bpca
improve_bpca_mar  = 100*(mar_mse$bpca - mar_mse$hurdle)/mar_mse$bpca 
improve_ppca_mar  = 100*(mar_mse$ppca - mar_mse$hurdle)/mar_mse$ppca 

I = 1000000
test_prob = rep(NA,30)
for(dd in 1:30){
  A = as.matrix(read.csv(file=paste0(file_root,"missing_data/simulated_data/MAR_data_",dd,".csv")))
  set1 = A[!is.na(A[,1]),1]
  set2 = read.csv(file=paste0(file_root,"missing_data/simulated_data/MAR_truth_",dd,".csv"))[,1]
  sample1 = sample(set1,size=I,replace=TRUE) 
  sample2 = sample(set2,size=I,replace=TRUE) 
  test_prob[dd] = mean(sample1 > sample2)
}

#pdf(paste0(file_root2,"improvement_plot.pdf"),height = 6, width = 6)
plot(test_prob*(1-test_prob), improve_bpca_mar,xlab="", ylab = "",ylim=c(-25,25),pch=16)
abline(h=0, lty="dashed",lwd=2)
#dev.off()


