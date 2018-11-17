

# This file performs ZIFA style analysis

#### Edit the below file_root to be the location of the saved folder:
file_root = "C:/........./paper_r_examples/"


## Perform reduced hurdle loss on zip_data using ZIFA type model

### read-in data:
zip_data = read.csv(paste0(file_root,"zero_inflated/zip_data.csv"))
### load some functions
source(file = paste0(file_root,"glrm_functions.R"))

### data prep:
n = nrow(zip_data)
p = ncol(zip_data)
A = as.matrix(zip_data)
H = matrix(1, ncol=p,nrow=n)
H[zip_data > 0] = 0
zifa_loss = zifa_sq_error = zifa_compare_1 = rep(NA,p)
# compute offsets, scaling, and null_loss:
Z_0 = matrix(0, ncol=p,nrow=n)
quad_offset = apply(zip_data,2,function(x) mean(x))
quad_sd = apply(zip_data,2,function(x) sd(x))
initial_decay_coefs = zifa_offset = zifa_scale = rep(1,p)
for(pp in 1:p){
    initial_decay_coefs[pp] = optimize_lambda(lambda=1, X=rep(quad_offset[pp],n), H=H[,pp])
#   tmp = optimize_pair(mu = quad_offset[pp], lambda = 1, A = A[,pp], H=H[,pp], tol = 0.01, maxiter = 100)
#   initial_decay_coefs[pp] = tmp[2] 
#   zifa_offset[pp]         = tmp[1]
}
zifa_offset = quad_offset
for(pp in 1:p){zifa_scale[pp] = zifa_column_loss(A = A[,pp],H=H[,pp],X=Z_0[,pp],offset=zifa_offset[pp],scale=1,decay=initial_decay_coefs[pp])}
zifa_scale = zifa_scale/(n-1)
my_null_loss = rep(NA,p)
for(pp in 1:p){my_null_loss[pp] = zifa_column_loss(A = A[,pp],H=H[,pp],X=Z_0[,pp],offset=zifa_offset[pp],scale=zifa_scale[pp],decay=initial_decay_coefs[pp])}
my_null_loss = sum(my_null_loss)

# ### Alternating Minimization
# ### setup parallel backend:
# library(doParallel)
# library(foreach)
# num_cores <- detectCores() - 1
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)
# ### initial conditions: 
# #s = svd(scale(A, center=colMeans(A), scale=apply(A,2,sd)))
# s = svd(scale(A, center=colMeans(A),scale=FALSE))
# current_decay = initial_decay_coefs
# 
# ### k = 2 # reduced rank dimension
# for(k in 2:14){
#   Xt = s$u[,1:k] %*% diag(s$d[1:k]) 
#   Yt = t(s$v[,1:k])  #%*% diag(apply(A,2,sd))
#   ### get starting loss:
#   starting_loss = rep(NA,p)
#   Zt = Xt%*%Yt
#   for(pp in 1:p){starting_loss[pp] = zifa_column_loss(A = A[,pp],H=H[,pp],X=Zt[,pp],offset=zifa_offset[pp],scale=zifa_scale[pp],decay=current_decay[pp])}
#   current_loss = sum(starting_loss)/my_null_loss
#   obj = 10
#   my_tol = .01
#   while(obj > my_tol){
#     last_loss = current_loss
#     ### update Xt:
#     Xt <- foreach( i=1:n, .combine=rbind) %dopar% {zifa_predict_x_nm(new_a=A[i,],h=H[i,],old_x=Xt[i,],old_y=Yt,muj=zifa_offset,sigj=zifa_scale,decay=current_decay)}
#     ### update Yt:
#     Yt <- foreach( j=1:p, .combine=cbind) %dopar% {zifa_predict_y_nm(new_a=A[,j],h=H[,j],old_x=Xt,old_y=Yt[,j],muj=zifa_offset[j],sigj=zifa_scale[j],decay=current_decay[j])}
#     ### update decay:
#     Zt = scale(Xt%*%Yt,center=-zifa_offset,scale=FALSE)
#     for(pp in 1:p){current_decay[pp] = optimize_lambda(lambda=current_decay[pp], X=Zt[,pp], H=H[,pp])} 
#     ### get current loss
#     current_loss = rep(NA,p)
#     Zt = Xt%*%Yt
#     for(pp in 1:p){current_loss[pp] = zifa_column_loss(A = A[,pp],H=H[,pp],X=Zt[,pp],offset=zifa_offset[pp],scale=zifa_scale[pp],decay=current_decay[pp])}
#     current_loss = 100 - 100*sum(current_loss)/my_null_loss
#     obj = last_loss - current_loss
#     print(current_loss)
#   }
#   zifa_loss[k]      = current_loss
#   Zt = scale(Xt%*%Yt, center=-zifa_offset, scale=FALSE)
#   p0 = matrix(0, ncol=p,nrow=n)
#   for(pp in 1:p){p0[,pp] = exp(-current_decay[pp]*(Zt[,pp]^2))}
#   Ahat = Zt
#   Ahat[p0 > 0.50]   = 0
#   zifa_sq_error[k]  = mean(scale(Ahat - A, center=FALSE, scale = quad_sd)^2)  
#   zifa_compare_1[k] = mean(abs(H - ifelse(p0 > 0.5, 1, 0)))
# }


### initialized at SVD solution and tuning:
Xt = s$u[,1] *s$d[1] 
Yt = t(s$v[,1])
Zt = scale(Xt%*%Yt,center=-zifa_offset,scale=FALSE) 
for(pp in 1:p){current_decay[pp] = optimize_lambda(lambda=1, X=Zt[,pp], H=H[,pp])} 
current_loss = rep(NA,p)
for(pp in 1:p){current_loss[pp] = zifa_column_loss(A = A[,pp],H=H[,pp],X=Zt[,pp],offset=zifa_offset[pp],scale=zifa_scale[pp],decay=current_decay[pp])}
zifa_loss[1] = 100 - 100*sum(current_loss)/my_null_loss
p0 = matrix(0, ncol=p,nrow=n)
for(pp in 1:p){p0[,pp] = exp(-current_decay[pp]*(Zt[,pp]^2))}
Ahat = Zt
Ahat[p0 > 0.50]   = 0
zifa_sq_error[1]  = mean(scale(Ahat - A, center=FALSE, scale = quad_sd)^2)  
zifa_compare_1[1] = mean(abs(H - ifelse(p0 > 0.5, 1, 0)))
zifa_loss = zifa_sq_error = zifa_compare_1 = rep(NA,p)
Xt = s$u[,1] *s$d[1] 
Yt = t(s$v[,1])
Zt = scale(Xt%*%Yt,center=-zifa_offset,scale=FALSE) 
for(pp in 1:p){current_decay[pp] = optimize_lambda(lambda=1, X=Zt[,pp], H=H[,pp])} 
current_loss = rep(NA,p)
for(pp in 1:p){current_loss[pp] = zifa_column_loss(A = A[,pp],H=H[,pp],X=Zt[,pp],offset=zifa_offset[pp],scale=zifa_scale[pp],decay=current_decay[pp])}
zifa_loss[1] = 100 - 100*sum(current_loss)/my_null_loss
p0 = matrix(0, ncol=p,nrow=n)
for(pp in 1:p){p0[,pp] = exp(-current_decay[pp]*(Zt[,pp]^2))}
Ahat = Zt
Ahat[p0 > 0.50]   = 0
zifa_sq_error[1]  = mean(scale(Ahat - A, center=FALSE, scale = quad_sd)^2)  
zifa_compare_1[1] = mean(abs(H - ifelse(p0 > 0.5, 1, 0)))
for(k in 2:14){
  Xt = s$u[,1:k] %*% diag(s$d[1:k]) 
  Yt = t(s$v[,1:k])
  Zt = scale(Xt%*%Yt,center=-zifa_offset,scale=FALSE) 
  for(pp in 1:p){current_decay[pp] = optimize_lambda(lambda=1, X=Zt[,pp], H=H[,pp])} 
  current_loss = rep(NA,p)
  for(pp in 1:p){current_loss[pp] = zifa_column_loss(A = A[,pp],H=H[,pp],X=(Xt%*%Yt)[,pp],offset=zifa_offset[pp],scale=zifa_scale[pp],decay=current_decay[pp])}
  zifa_loss[k] = 100 - 100*sum(current_loss)/my_null_loss
  p0 = matrix(0, ncol=p,nrow=n)
  for(pp in 1:p){p0[,pp] = exp(-current_decay[pp]*(Zt[,pp]^2))}
  Ahat = Zt
  Ahat[p0 > 0.50]   = 0
  zifa_sq_error[k]  = mean(scale(Ahat - A, center=FALSE, scale = quad_sd)^2)  
  zifa_compare_1[k] = mean(abs(H - ifelse(p0 > 0.5, 1, 0)))
}
aaa = zifa_loss
bbb = zifa_sq_error
ccc = zifa_compare_1
plot(0:14,c(0,zifa_loss)/100,type="l",ylim=c(0,1))
plot(1:14,zifa_sq_error,type="l")
plot(1:14,zifa_compare_1,type="l",ylim=c(0,0.25))
zifa_results = list(zifa_loss = c(0,aaa)/100, zifa_sq_error=bbb, zifa_compare_1=ccc)
#save(zifa_results,file=paste0(file_root,"zero_inflated/stored_results/zifa_results.RData"))

