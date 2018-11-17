

# This file performs cross validation to select the penalty parameter values for
# the missing value data sets.

#### Edit the below file_root to be the location of the saved folder:
file_root = "C:/........./paper_r_examples/"

source(file = paste0(file_root, "glrm_functions.R"))

# ----------------------------------------- #
#      Choose regularization parameter
# ----------------------------------------- #

### setup parallel backend:
library(doParallel)
library(foreach)
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
### Tuning:
xmethod = "nm"
ymethod = "nm"
NT        = 50
stop_rule = 0.1
k = 4
### start optimization:
gamma_list = c(4:80)/4
gamma_min = rep(NA,30)
error_min = rep(NA,30)
set.seed(123)
for(dd in 1:30){
  print(paste0("##### ",dd," #####"))
  # MCAR case:
  # A = as.matrix(read.csv(file=paste0(file_root, "missing_data/simulated_data/MCAR_data_",dd,".csv")))
  # MAR case:
  A = as.matrix(read.csv(file=paste0(file_root, "missing_data/simulated_data/MAR_data_",dd,".csv")))
  na_prob = mean(is.na(A[,1]))
  A = A[!is.na(A[,1]), ]
  m = nrow(A)
  my_na = rbinom(m, size=1, prob=na_prob)
  my_obs = 1 - my_na 
  keep_for_later = A[my_na==1,1]
  A[my_na==1,1] = NA
  B = prep_data(A = A)
  A_embed = embed_data(A = B, embed_col = 1, how = "Hurdle", value = NA)
  loss = c("Quadratic", "Logistic", rep("Quadratic",9))
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
  reconstrution_error = rep(NA,length(gamma_list))
  for(gl in 1:length(gamma_list)){
    print(gamma_list[gl])
    gamma_x = gamma_y = gamma_list[gl]
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
      offset[1] = mean((A_embed[,1] - Xt%*%Yt[,1])[my_obs == 1])
      # update Yt
      if(ymethod =="nm"){Yt <- foreach( j=1:n, .combine=cbind) %dopar% {
        predict_y_nm(loss=loss[j],new_a=A_embed[,j],old_x=Xt,old_y=Yt[,j],muj=offset[j],sigj=scale[j],gamma_y=gamma_y)
      }}
      t = t+1
      obj_fct[t] = sum(generic_loss(loss=loss, offset=offset, scale=scale, A=A_embed, XY=Xt%*%Yt), na.rm = TRUE) +
        gamma_x*sum(Xt^2) + gamma_y*sum(Yt^2)
      obj_diff = obj_fct[t-1] - obj_fct[t]
    }
    print(t-1)
    my_reconst = A_embed
    for(nn in 1:ncol(Yt)){my_reconst[,nn] = my_reconstruct(loss = loss[nn], X = Xt, Y = Yt[,nn], offset = offset[nn])}
    A_hat = my_reconst[,-2]
    reconstrution_error[gl] = mean((A_hat[my_na == 1,1] - keep_for_later)^2)
  }
  my_min = which.min(reconstrution_error)[1]
  error_min[dd] = reconstrution_error[my_min]
  gamma_min[dd] = gamma_list[my_min]
  print(gamma_min[dd])
  print(error_min[dd])
}
penalty_summary_MCAR = list(gamma_min=gamma_min, error_min=error_min)
file_name = paste0(file_root,"missing_data/stored_results/penalty_summary_MCAR.RData")
#save(penalty_summary_MCAR,file = file_name)
penalty_summary_MAR = list(gamma_min=gamma_min, error_min=error_min)
file_name = paste0(file_root,"missing_data/stored_results/penalty_summary_MAR.RData")
#save(penalty_summary_MAR,file = file_name)
cbind(gamma_min, error_min)

### shutdown backend
stopCluster(cl)
