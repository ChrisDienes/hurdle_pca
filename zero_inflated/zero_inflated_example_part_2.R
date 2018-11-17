
# This file summarizes the ZIP analysis and creates objects for final paper 

#### Edit the below file_root to be the location of the saved folder:
file_root = "C:/........./paper_r_examples/"

### Load requirements:
source(file = paste0(file_root,"glrm_functions.R"))
load(paste0(file_root, "zero_inflated/stored_results/pca.Rdata"))
load(file=paste0(file_root, "zero_inflated/stored_results/zifa_results.RData"))
zip_data = read.csv(paste0(file_root,"zero_inflated/zip_data.csv"))
# Prep data structures:
C   = matrix(1, ncol=ncol(zip_data),nrow=nrow(zip_data))
C[zip_data > 0] = 0

### ---------------------------------------------------------------------- ###
### create cum_var.pdf:
### ---------------------------------------------------------------------- ###

pca_cum_var = pca_saved$cum_var
hurdle_cum_var = rep(0,14)
for(kk in 1:14){
  load(paste0(file_root,"zero_inflated/stored_results/hurdle/hurdle_",kk,".Rdata"))
  hurdle_cum_var[kk] = 1 - hurdle_saved$model_loss/30786
}
#pdf(paste0(file_root,"zero_inflated/stored_results/plots/cum_var.pdf"),height = 6, width = 6)
plot(0:14, c(0,hurdle_cum_var), xlim = c(0,14), ylim= c(0,1), type="l", xlab = "", ylab = "",lwd=2)
lines(0:14, c(0,pca_cum_var), lty="dashed",lwd=2)
lines(0:14, zifa_results$zifa_loss, lty="dotdash",lwd=3)
#dev.off()

### ---------------------------------------------------------------------- ###
### comapare various reconstruction errors over range of k
### ---------------------------------------------------------------------- ###
A = as.matrix(zip_data)
B = prep_data(A = A)
nn = ncol(A)
A_embed = embed_data(A = B, embed_col = 1:nn, how = rep("Hurdle",nn), value = rep(0,nn))
embed_offset = colMeans(A_embed, na.rm = TRUE)
embed_scale = apply(A_embed, 2, sd, na.rm = TRUE)
load(paste0(file_root,"zero_inflated/stored_results/pca.Rdata"))
pca_sq_error = pca_compare_1 = rep(0,14)
for(kk in 1:14){
  if(kk == 1){Ahat   = scale( matrix(pca_saved$Xt[,1:kk],ncol=1)%*%matrix(pca_saved$Yt[1:kk,],nrow=1), center = -pca_saved$xbar, scale=FALSE)}
  if(kk > 1){Ahat    = scale( pca_saved$Xt[,1:kk]%*%pca_saved$Yt[1:kk,], center = -pca_saved$xbar, scale=FALSE)}
  pca_sq_error[kk]   = mean(scale(Ahat - A, center=FALSE, scale = pca_saved$std)^2)
  pca_compare_1[kk]  = mean(abs(C - ifelse(Ahat < 0.5, 1, 0))) 
}
loss = c(rep(c("Poisson", "Logistic"),nn))
offset = compute_offset(A = A_embed, loss = loss)
scale = 1:(2*nn)
for(jj in 1:nn){
  scale[(2*jj -1):(2*jj)] = 1/scale_zero_hurdle_eqn_36_v2(logistic_data=A_embed[,2*jj], poisson_data=A_embed[,2*jj-1], 
                                                          logistic_offset=offset[2*jj], poisson_offset=offset[2*jj -1], 
                                                          c="n_v")[2:1]
}
hurdle_sq_error = hurdle_compare_1 = hurdle_embed_sq_error = rep(0,14)
for(kk in 1:14){
  load(paste0(file_root,"zero_inflated/stored_results/hurdle/hurdle_",kk,".Rdata"))
  Yt = hurdle_saved$Yt
  Xt = hurdle_saved$Xt
  if(kk == 1){Zt = matrix(Xt,ncol=1)%*%matrix(Yt,nrow=1)}
  if(kk > 1){Zt = Xt%*%Yt}
  A_hat = A
  for(nn in 1:ncol(A)){
    A_hat[,nn] = hurdle_count_reconstruct(logistic_var=Zt[,2*nn], pois_var=Zt[,2*nn -1], logistic_offset=offset[2*nn], 
                           pois_offset=offset[2*nn-1], scale1=scale[2*nn], scale2=scale[2*nn-1], nu = 0)
  }
  A_hat_embed = scale(Zt, center = -offset, scale=FALSE)
  for(tt in 1:14){
    A_hat_embed[,2*tt]    = 2/(1+exp(-A_hat_embed[,2*tt])) - 1
    A_hat_embed[,2*tt -1] = exp(A_hat_embed[,2*tt -1])
  }
  hurdle_embed_sq_error[kk] = mean(scale(A_hat_embed - A_embed, center=FALSE, scale = embed_scale)^2, na.rm=TRUE)
  hurdle_sq_error[kk] = mean(scale(A_hat - A, center=FALSE, scale = pca_saved$std)^2)
  hurdle_compare_1[kk]  = mean(abs(C - ifelse(A_hat < 0.5, 1, 0))) 
}

ymax = max(c(max(hurdle_sq_error),max(pca_sq_error),max(zifa_results$zifa_sq_error)))
#pdf(paste0(file_root,"zero_inflated/stored_results/plots/reconst_original_sse.pdf"),height = 6, width = 6)
plot(1:14, hurdle_sq_error, xlim = c(1,14), ylim = c(0,ymax), type="l", xlab = "", ylab = "",lwd=2)
lines(1:14, pca_sq_error, lty="dashed",lwd=2)
lines(1:14, zifa_results$zifa_sq_error, lty="dotdash",lwd=3)
#dev.off()
ymax = max(c(max(hurdle_compare_1),max(pca_compare_1),max(zifa_results$zifa_compare_1)))
#pdf(paste0(file_root,"zero_inflated/stored_results/plots/compare_above_half.pdf"),height = 6, width = 6)
plot(1:14, hurdle_compare_1, xlim = c(1,14), ylim = c(0,ymax), type="l", xlab = "", ylab = "",lwd=2)
lines(1:14, pca_compare_1, lty="dashed",lwd=2)
lines(1:14, zifa_results$zifa_compare_1, lty="dotdash",lwd=3)
#dev.off()

