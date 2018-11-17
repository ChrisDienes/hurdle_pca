
# This file performs the hurdle analysis for all 30 missing value data sets 

#### Edit the below file_root to be the location of the saved folder:
file_root = "C:/........./paper_r_examples/"
source(file = paste0(file_root, "glrm_functions.R"))

library(ROCR)

### load other PCA results:
file_root2 = paste0(file_root,"missing_data/stored_results/")
load(file=paste0(file_root2,"MCAR_MSE_df.RData"))
load(file=paste0(file_root2,"MCAR_offset_df.RData"))
load(file=paste0(file_root2,"MAR_MSE_df.RData"))
load(file=paste0(file_root2,"MAR_offset_df.RData"))

### Prepare hurdle settings:
xmethod = "nm"
ymethod = "nm"
NT        = 250
stop_rule = 0.05
load(file=paste0(file_root2,"penalty_summary_MCAR.RData"))
load(file=paste0(file_root2,"penalty_summary_MAR.RData"))
k = 4
m = 5000
loss = c("Quadratic", "Logistic", rep("Quadratic",9))
hurdle_offset_MCAR = hurdle_offset_MAR = rep(NA, 30)
hurdle_MSE_MCAR = hurdle_MSE_MAR = rep(NA, 30)
hurdle_AUC_MCAR =hurdle_AUC_MAR = rep(NA,30)
hurdle_theta_MAR = as.data.frame(matrix(NA, ncol=11, nrow=30))
hurdle_dist_MAR  = as.data.frame(matrix(NA, ncol=11, nrow=30))  
names(hurdle_theta_MAR) = c("V1","V1_na", paste0("V", 2:10))
names(hurdle_dist_MAR) = c("V1","V1_na", paste0("V", 2:10))
### setup parallel backend:
library(doParallel)
library(foreach)
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
### Cycle trhough the data sets:
for(dd in 1:30){
  print(paste0("#### ", dd, " ####"))
  ### MCAR case:
  gamma_x = gamma_y = penalty_summary_MCAR$gamma_min[dd]
  print(paste0(dd," MCAR"))
  A = as.matrix(read.csv(file=paste0(file_root,"missing_data/simulated_data/MCAR_data_",dd,".csv")))
  MCAR_truth = read.csv(file=paste0(file_root,"missing_data/simulated_data/MCAR_truth_",dd,".csv"))
  my_na = ifelse(is.na(A[,1]),1,0)
  B = prep_data(A = A)
  A_embed = embed_data(A = B, embed_col = 1, how = "Hurdle", value = NA)
  offset = compute_offset(A = A_embed, loss = loss)
  scale = 1:length(offset)
  scale[3:length(offset)] = compute_scale(A = A_embed, loss = loss, offset = offset, nreturn = TRUE)[3:length(offset),1]
  scale[1:2] = 1/(scale_NA_hurdle_eqn_36_v2(logistic_data=A_embed[,2], quad_data=A_embed[,1],
                                          logistic_offset=offset[2], quad_offset=offset[1],
                                          c="n_v")[2:1])
  my_null_loss = null_loss(A=A_embed, loss=loss, offset=offset, scale=scale)
  n = ncol(A_embed)
  ### initial conditions:
  ### Random
  Yt = matrix(rnorm(k*n),nrow=k,ncol=n)
  Xt = matrix(1,nrow=m,ncol=k)
  ### start optimization:
  obj_fct = c(sum(generic_loss(loss=loss, offset=offset, scale=scale, A=A_embed, XY=Xt%*%Yt), na.rm=TRUE) +
                gamma_x*sum(Xt^2) + gamma_y*sum(Yt^2),rep(NA, NT))
  obj_diff = obj_fct[1]
  obj_diff
  t = 1
  while(t<=NT & obj_diff > stop_rule){
    #if(t == 1){cat(c("Iter","Loss Change", "Current Loss", "offset", "\n"))}
    #cat(c(t, obj_diff,obj_fct[t],offset[1],"\n"))
    # update Xt
    if(xmethod == "nm"){Xt <- foreach( i=1:m, .combine=rbind) %dopar% {
      predict_x_nm(loss=loss,new_a=A_embed[i,],old_x=Xt[i,],old_y=Yt,muj=offset,sigj=scale,gamma_x=gamma_x)
    }}
    # update offset
    offset[1] = mean((A_embed[,1] - Xt%*%Yt[,1])[my_na == 0])
    # update Yt
    if(ymethod =="nm"){Yt <- foreach( j=1:n, .combine=cbind) %dopar% {
      predict_y_nm(loss=loss[j],new_a=A_embed[,j],old_x=Xt,old_y=Yt[,j],muj=offset[j],sigj=scale[j],gamma_y=gamma_y)
    }}
    t = t+1
    obj_fct[t] = sum(generic_loss(loss=loss, offset=offset, scale=scale, A=A_embed, XY=Xt%*%Yt), na.rm = TRUE) +
      gamma_x*sum(Xt^2) + gamma_y*sum(Yt^2)
    obj_diff = obj_fct[t-1] - obj_fct[t]
  }
  ### percent reduction on baseline model:
  my_model_loss = sum(generic_loss(loss=loss, offset=offset, scale=scale, A=A_embed, XY=Xt%*%Yt), na.rm = TRUE)
  #print(paste0("Model loss: " ,100*(1 - my_model_loss/my_null_loss)))
  hurdle_offset_MCAR[dd] = offset[1]
  my_reconst = A_embed
  for(nn in 1:ncol(Yt)){my_reconst[,nn] = my_reconstruct(loss = loss[nn], X = Xt, Y = Yt[,nn], offset = offset[nn])}
  A_hat = my_reconst[,-2]
  hurdle_MSE_MCAR[dd] = mean((A_hat[my_na == 1,1] - MCAR_truth)^2)
  print(hurdle_MSE_MCAR[dd])
  prob_na = my_reconst[,2]
  my_pred = prediction(predictions=prob_na, labels=my_na)
  auc = performance(my_pred, measure = "auc")
  hurdle_AUC_MCAR[dd] = auc@y.values[[1]][1]
  ### MAR case:
  gamma_x = gamma_y = penalty_summary_MAR$gamma_min[dd]
  print(paste0(dd," MAR"))
  A = as.matrix(read.csv(file=paste0(file_root,"missing_data/simulated_data/MAR_data_",dd,".csv")))
  MAR_truth = read.csv(file=paste0(file_root,"missing_data/simulated_data/MAR_truth_",dd,".csv"))
  my_na = ifelse(is.na(A[,1]),1,0)
  B = prep_data(A = A)
  A_embed = embed_data(A = B, embed_col = 1, how = "Hurdle", value = NA)
  offset = compute_offset(A = A_embed, loss = loss)
  scale = 1:length(offset)
  scale[3:length(offset)] = compute_scale(A = A_embed, loss = loss, offset = offset, nreturn = TRUE)[3:length(offset),1]
  scale[1:2] = 1/(scale_NA_hurdle_eqn_36_v2(logistic_data=A_embed[,2], quad_data=A_embed[,1], 
                                            logistic_offset=offset[2], quad_offset=offset[1], 
                                            c="n_v")[2:1])
  my_null_loss = null_loss(A=A_embed, loss=loss, offset=offset, scale=scale)
  n = ncol(A_embed)
  ### initial conditions: 
  ### Random
  Yt = matrix(rnorm(k*n),nrow=k,ncol=n)
  Xt = matrix(1,nrow=m,ncol=k)
  ### start optimization:
  obj_fct = c(sum(generic_loss(loss=loss, offset=offset, scale=scale, A=A_embed, XY=Xt%*%Yt), na.rm=TRUE) + 
                gamma_x*sum(Xt^2) + gamma_y*sum(Yt^2),rep(NA, NT))
  obj_diff = obj_fct[1]
  obj_diff
  t = 1
  while(t<=NT & obj_diff > stop_rule){
    #if(t == 1){cat(c("Iter","Loss Change", "Current Loss", "offset", "\n"))}
    #cat(c(t, obj_diff,obj_fct[t],offset[1],"\n"))
    # update Xt
    if(xmethod == "nm"){Xt <- foreach( i=1:m, .combine=rbind) %dopar% {
      predict_x_nm(loss=loss,new_a=A_embed[i,],old_x=Xt[i,],old_y=Yt,muj=offset,sigj=scale,gamma_x=gamma_x)
    }}
    # update offset
    offset[1] = mean((A_embed[,1] - Xt%*%Yt[,1])[my_na == 0])
    # update Yt
    if(ymethod =="nm"){Yt <- foreach( j=1:n, .combine=cbind) %dopar% {
      predict_y_nm(loss=loss[j],new_a=A_embed[,j],old_x=Xt,old_y=Yt[,j],muj=offset[j],sigj=scale[j],gamma_y=gamma_y)
    }}
    t = t+1
    obj_fct[t] = sum(generic_loss(loss=loss, offset=offset, scale=scale, A=A_embed, XY=Xt%*%Yt), na.rm = TRUE) +
      gamma_x*sum(Xt^2) + gamma_y*sum(Yt^2)
    obj_diff = obj_fct[t-1] - obj_fct[t]
  }
  ### percent reduction on baseline model:
  my_model_loss = sum(generic_loss(loss=loss, offset=offset, scale=scale, A=A_embed, XY=Xt%*%Yt), na.rm = TRUE)
  #print(paste0("Model loss: " ,100*(1 - my_model_loss/my_null_loss)))
  hurdle_offset_MAR[dd] = offset[1]
  my_reconst = A_embed
  for(nn in 1:ncol(Yt)){my_reconst[,nn] = my_reconstruct(loss = loss[nn], X = Xt, Y = Yt[,nn], offset = offset[nn])}
  A_hat = my_reconst[,-2]
  hurdle_MSE_MAR[dd] = mean((A_hat[my_na == 1,1] - MAR_truth)^2)
  print(hurdle_MSE_MAR[dd])
  prob_na = my_reconst[,2]
  my_pred = prediction(predictions=prob_na, labels=my_na)
  auc = performance(my_pred, measure = "auc")
  hurdle_AUC_MAR[dd] = auc@y.values[[1]][1]
  gram_y = t(Yt)%*%Yt
  gram_y = gram_y*(diag(gram_y)^-.5)%*%t(diag(gram_y)^-.5)
  cosine_sim = 1 - acos((pmin(pmax(gram_y,-1.0),1.0)))/pi
  diag(cosine_sim) <- 1
  mydist = 1 - abs(2*(cosine_sim - 0.5))
  hurdle_theta_MAR[dd,] = cosine_sim[,2]
  hurdle_dist_MAR[dd,] = mydist[,2]
} 

# save(hurdle_MSE_MCAR, file=paste0(file_root2,"hurdle_MSE_MCAR.RData"))
# save(hurdle_offset_MCAR, file=paste0(file_root2,"hurdle_offset_MCAR.RData"))
# save(hurdle_AUC_MCAR, file=paste0(file_root2,"hurdle_AUC_MCAR.RData"))
# save(hurdle_MSE_MAR, file=paste0(file_root2,"hurdle_MSE_MAR.RData"))
# save(hurdle_offset_MAR, file=paste0(file_root2,"hurdle_offset_MAR.RData"))
# save(hurdle_AUC_MAR, file=paste0(file_root2,"hurdle_AUC_MAR.RData"))
# save(hurdle_theta_MAR, file=paste0(file_root2,"hurdle_theta_MAR.RData"))
# save(hurdle_dist_MAR, file=paste0(file_root2,"hurdle_dist_MAR.RData"))

### shutdown backend
stopCluster(cl)
