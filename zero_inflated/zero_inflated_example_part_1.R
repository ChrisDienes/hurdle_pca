
# This file performs hurdle analysis on the zip data set

#### Edit the below file_root to be the location of the saved folder:
file_root = "C:/........./paper_r_examples/"


# -------------------------------------------------- #
### Load custom glrm functions:
source(file = paste0(file_root, "glrm_functions.R"))

### Load example data:
zip_data = read.csv(paste0(file_root,"zero_inflated/zip_data.csv"))
A = as.matrix(zip_data)
m = nrow(A)
n = ncol(A)

# ---------------------------------------- #
#     Conduct regular PCA/SVD analysis:
# ---------------------------------------- #
xbar = colMeans(A)
std  = apply(A,2,sd)
s = svd(scale(A, center=xbar, scale=std))
pca_cum_var = cumsum(s$d^2/sum(s$d^2))
plot(1:length(pca_cum_var), pca_cum_var, type = "l", ylim=c(0,1),xlab="Number of PCs",ylab="Variance Explained")
Xt = s$u[,1:14] %*% diag(s$d[1:14]) 
Yt = t(s$v[,1:14])  %*% diag(std)
### Save data in stored results:
pca_saved = list(Xt = Xt, Yt = Yt, xbar = xbar, std = std, cum_var = pca_cum_var)
# save(pca_saved, file = paste0(file_root, "zero_inflated/stored_results/pca.RData")

# ---------------------------------------------------- #
#   Conduct full hurdle loss analysis of embedded data
# ---------------------------------------------------- #
B = prep_data(A = A)
nn = ncol(A)
A_embed = embed_data(A = B, embed_col = 1:nn, how = rep("Hurdle",nn), value = rep(0,nn))
loss = c(rep(c("Poisson", "Logistic"),nn))
offset = compute_offset(A = A_embed, loss = loss)
scale = 1:(2*nn)
for(jj in 1:nn){
  scale[(2*jj -1):(2*jj)] = 1/scale_zero_hurdle_eqn_36_v2(logistic_data=A_embed[,2*jj], poisson_data=A_embed[,2*jj-1], 
                            logistic_offset=offset[2*jj], poisson_offset=offset[2*jj -1], 
                            c="n_v")[2:1]
}
my_null_loss = null_loss(A=A_embed, loss=loss, offset=offset, scale=scale)
### Check: sum(apply(A, 2, function(x) sum(!is.na(x))) - 1)
# -------------------------------------------------- #
### Alternating Minimization
### setup parallel backend:
library(doParallel)
library(foreach)
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
k = 2 # reduced rank dimension
### initial conditions: 
m = nrow(A_embed)
n = ncol(A_embed)
### Random
Yt = matrix(rnorm(k*n),nrow=k,ncol=n)
Xt = matrix(1,nrow=m,ncol=k)
### Tuning:
xmethod = "nm"
ymethod = "nm"
NT        = 500
stop_rule = 0.1
### start optimization:
for(k in 2:14){
  print(paste0("#########  ",k,"    #########"))
  if(k > 2){
     Yt = rbind(Yt, rnorm(n))
     Xt = cbind(Xt, rep(1,nrow=m))
  }
  #Yt = matrix(rnorm(k*n),nrow=k,ncol=n)
  #Xt = matrix(1,nrow=m,ncol=k)
  obj_fct = c(sum(generic_loss(loss=loss, offset=offset, scale=scale, A=A_embed, XY=Xt%*%Yt), na.rm=TRUE),rep(NA, NT))
  obj_diff = obj_fct[1]
  obj_diff
  t = 1
  while(t<=NT & obj_diff > stop_rule){
    if(t == 1){cat(c("Iter","Loss Change", "Current Loss", "\n"))}
    cat(c(t, obj_diff,obj_fct[t],"\n"))
    if(xmethod == "nm"){Xt <- foreach( i=1:m, .combine=rbind) %dopar% {predict_x_nm(loss=loss,new_a=A_embed[i,],old_x=Xt[i,],old_y=Yt,muj=offset,sigj=scale)}}
    if(ymethod =="nm"){Yt <- foreach( j=1:n, .combine=cbind) %dopar% {predict_y_nm(loss=loss[j],new_a=A_embed[,j],old_x=Xt,old_y=Yt[,j],muj=offset[j],sigj=scale[j])}}
    t = t+1
    obj_fct[t] = sum(generic_loss(loss=loss, offset=offset, scale=scale, A=A_embed, XY=Xt%*%Yt), na.rm = TRUE)
    obj_diff = obj_fct[t-1] - obj_fct[t]
  }
  cat(c(k, t, obj_diff,obj_fct[t],"\n"))
  my_model_loss = sum(generic_loss(loss=loss, offset=offset, scale=scale, A=A_embed, XY=Xt%*%Yt), na.rm = TRUE)
  hurdle_saved = list(Xt=Xt, Yt=Yt, model_loss = my_model_loss)
  file_name = paste0(file_root,"zero_inflated/stored_results/hurdle_",k,".RData")
  #save(hurdle_saved, file = file_name)
  print(100*(1 - my_model_loss/my_null_loss))
}

stopCluster(cl)

##### ------------------------------ #####
##### k = 1 case:                    #####
##### ------------------------------ #####
B = prep_data(A = A)
nn = ncol(A)
A_embed = embed_data(A = B, embed_col = 1:nn, how = rep("Hurdle",nn), value = rep(0,nn))
loss = c(rep(c("Poisson", "Logistic"),nn))
offset = compute_offset(A = A_embed, loss = loss)
scale = 1:(2*nn)
for(jj in 1:nn){
  scale[(2*jj -1):(2*jj)] = 1/scale_zero_hurdle_eqn_36_v2(logistic_data=A_embed[,2*jj], poisson_data=A_embed[,2*jj-1], 
                                                          logistic_offset=offset[2*jj], poisson_offset=offset[2*jj -1], 
                                                          c="n_v")[2:1]
}
my_null_loss = null_loss(A=A_embed, loss=loss, offset=offset, scale=scale)
### Check: sum(apply(A, 2, function(x) sum(!is.na(x))) - 1)
# -------------------------------------------------- #
### Alternating Minimization
### setup parallel backend:
library(doParallel)
library(foreach)
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)
k = 1 # reduced rank dimension
### initial conditions: 
m = nrow(A_embed)
n = ncol(A_embed)
### Random
Yt = matrix(rnorm(m*n),nrow=k,ncol=n)
Xt = matrix(1,nrow=m,ncol=k)
### Tuning:
xmethod = "nm"
ymethod = "nm"
NT        = 500
stop_rule = 0.1
obj_fct = c(sum(generic_loss(loss=loss, offset=offset, scale=scale, A=A_embed, XY=Xt%*%Yt), na.rm=TRUE),rep(NA, NT))
obj_diff = obj_fct[1]
obj_diff
t = 1
while(t<=NT & obj_diff > stop_rule){
  if(t == 1){cat(c("Iter","Loss Change", "Current Loss", "\n"))}
  cat(c(t, obj_diff,obj_fct[t],"\n"))
  if(xmethod == "nm"){Xt <- foreach( i=1:m, .combine=rbind) %dopar% {predict_x_nm_k1(loss=loss,new_a=A_embed[i,],old_x=Xt[i,],old_y=matrix(Yt,nrow=1),muj=offset,sigj=scale)}}
  if(ymethod =="nm"){Yt <- foreach( j=1:n, .combine=cbind) %dopar% {predict_y_nm_k1(loss=loss[j],new_a=A_embed[,j],old_x=matrix(Xt,ncol=1),old_y=Yt[,j],muj=offset[j],sigj=scale[j])}}
  t = t+1
  obj_fct[t] = sum(generic_loss(loss=loss, offset=offset, scale=scale, A=A_embed, XY=Xt%*%Yt), na.rm = TRUE)
  obj_diff = obj_fct[t-1] - obj_fct[t]
}
my_model_loss = sum(generic_loss(loss=loss, offset=offset, scale=scale, A=A_embed, XY=Xt%*%Yt), na.rm = TRUE)
hurdle_saved = list(Xt=Xt, Yt=Yt, model_loss = my_model_loss)
#file_name = paste0(file_root,"zero_inflated/stored_results/hurdle_1.RData")
#save(hurdle_saved, file = file_name)
print(100*(1 - my_model_loss/my_null_loss))

### shutdown backend
stopCluster(cl)
