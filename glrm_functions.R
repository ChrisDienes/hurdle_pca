

# The below helper functions are used in the other R files:

# --------------------------------- #
#         Generic Wrapers           #
# --------------------------------- #

prep_data = function(A){
  colnames(A)  = paste0("A_", 1:ncol(A), ".0")
  return(A)
}

embed_data = function(A, embed_col, how, value){
  N = ncol(A)
  M = nrow(A)
  for(ec in 1:length(embed_col)){
    if(how[ec] == "Hurdle"){
      if(is.na(value[ec])){tmp = which(is.na(A[,embed_col[ec]]))}
      if(!is.na(value[ec])){tmp = which(A[,embed_col[ec]] == value[ec])}
      a1 = rep(-1, M)
      a1[tmp] =1
      A[tmp,paste0("A_",embed_col[ec],".0")] = NA
      A = cbind(A, a1)
      N = N + 1
      colnames(A)[N] <- paste0("A_",embed_col[ec],".1")
    }
  }
  A = A[,order(as.numeric(gsub("A_","",colnames(A))))]
  return(A)
}


generic_loss = function(loss, offset, scale, A, XY){
  flag = is.null(nrow(A))
  n = length(offset)
  m = ifelse(flag, ifelse(n==1, length(A), 1) ,nrow(A))
  if(flag){
    A  = matrix(A, ncol=n)
    XY = matrix(XY, ncol=n) 
  }
  loss_mat = matrix(NA, nrow=m, ncol=n)
  N = length(loss)
  for(nn in 1:N){
    if(loss[nn] == "Poisson"){loss_mat[,nn] = (1/scale[nn]) * pois_loss(u = XY[,nn] + offset[nn], a = A[,nn])}
    if(loss[nn] == "Quadratic"){loss_mat[,nn] = (1/scale[nn]) * l2_loss(u = XY[,nn] + offset[nn], a = A[,nn])}
    if(loss[nn] == "Logistic"){loss_mat[,nn] = (1/scale[nn]) * logistic_loss(u = XY[,nn] + offset[nn], a = A[,nn])}
    if(loss[nn] == "ZTQR"){
      tmp = which(A[,nn] > 0)
      tmp_a = ifelse(A[,nn] == 0, 1, -1)
      lambda = 1
      tmp_b = rep(0, m)
      tmp_b[tmp] = lambda*l2_loss(u = XY[tmp, nn] + offset[nn], a = A[tmp, nn])
      loss_mat[,nn] = (1/scale[nn]) *(logistic_loss(u = XY[,nn] + offset[nn], a = tmp_a) + tmp_b)
    }
    if(loss[nn] == "ZTPF"){
      tmp =which(!is.na(A[,nn]))
      loss_mat[tmp,nn] = (1/scale[nn]) * abs(hurdle_count_loss(u = XY[tmp,nn] + offset[nn], a = A[tmp,nn]))
    }
  }
  return(loss_mat)
}

compute_offset = function(A, loss, small_num =  10^-10, max_iter = 100){
  N = length(loss)
  hold = 1:N
  for(nn in 1:N){
    tmp = A[,nn]
    tmp = tmp[!is.na(tmp)]
    if(loss[nn] == "Quadratic"){hold[nn] = mean(tmp)}
    if(loss[nn] == "Poisson"){hold[nn] = log(mean(tmp))}
    if(loss[nn] == "Logistic"){
      n1 = sum(tmp == 1)
      n2 = sum(tmp == -1)
      hold[nn] = log(n1) - log(n2)
    }
    if(loss[nn] == "ZTQR"){
      a1 = ifelse(tmp == 0, 1, -1)
      a2 = tmp[tmp != 0]
      lambda = 1
      fn = function(mu, a1, a2, lambda){ sum(log(1 + exp(-a1*mu))) + lambda*sum((a2 - mu)^2)}
      gr = function(mu, a1, a2, lambda){ sum(-a1/(1+exp(a1*mu))) - 2*lambda*sum(a2-mu)}
      hold[nn] = optim(par=mean(a2), fn=fn, gr=gr, control = list(maxit=500), method="BFGS",a1=a1,a2=a2,lambda=lambda)$par
    }
    if(loss[nn] == "ZTPF"){
      a = mean(tmp)
      u0 = a
      test = 1
      last_iter = 0
      while(abs(test) > small_num & last_iter <= max_iter){
        pp = exp(-u0)
        top = a*(1 - pp) - u0
        bot = a*pp - 1
        un = u0 - top / bot
        test = top
        u0 = un
        last_iter = last_iter +1
      }
      hold[nn] = log(un)
    }
    if(loss[nn] == "ZTPR"){
      set0 = which(tmp == 0)
      n1 = length(set0)
      n = length(tmp)
      n2 =  n - n1
      a2 = sum(tmp[-set0])
      u0 = (n1/n)*log(n1/n2) + (n2/n)*log(mean(tmp))
      test = 1
      last_iter = 0
      while(abs(test) > small_num & last_iter <= max_iter){
        p1 = exp(u0)
        p2 = exp(-exp(u0))
        top = (n2*p1 - n1)/(1+p1) - a2 + (n2*p1)/(1-p2)
        bot = n*p1/((1+p1)^2) +(n2*p1*(1-2*p2)) / ( ((1-p2))^2 )
        un = u0 - top/bot
        test = top
        u0 = un
        last_iter = last_iter +1
      }
      hold[nn] = un
    }
  }
  return(hold)
}

scale_zero_hurdle_eqn_36_v2 = function(logistic_data, poisson_data, logistic_offset, poisson_offset, c){
  n_j = sum(!is.na(logistic_data))
  if(c == "n_v"){
    n_v = n_j - sum(!is.na(poisson_data))
    c = n_v/(n_j - n_v)
  }
  L_b = sum(logistic_loss(u = logistic_offset, a = logistic_data[!is.na(logistic_data)]))
  L_g = sum(pois_loss(u = poisson_offset, a = poisson_data[!is.na(poisson_data)]))
  lambda = solve(matrix(c(L_b, L_b, L_g, (-c)*L_g), nrow=2,ncol=2))%*%matrix(c(n_j - 1, 0),ncol=1)
  return(lambda)
}

scale_NA_hurdle_eqn_36_v2 = function(logistic_data, quad_data, logistic_offset, quad_offset, c){
  n_j = sum(!is.na(logistic_data))
  if(c == "n_v"){
    n_v = n_j - sum(!is.na(quad_data))
    c = n_v/(n_j - n_v)
  }
  L_b = sum(logistic_loss(u = logistic_offset, a = logistic_data[!is.na(logistic_data)]))
  L_g = sum(l2_loss(u = quad_offset, a = quad_data[!is.na(quad_data)]))
  lambda = solve(matrix(c(L_b, L_b, L_g, (-c)*L_g), nrow=2,ncol=2))%*%matrix(c(n_j - 1, 0),ncol=1)
  return(lambda)
}

compute_scale = function(A, loss, offset, nreturn = FALSE){
  N = length(loss)
  hold = 1:N
  hold_n = 1:N
  for(nn in 1:N){
    tmp = A[,nn]
    tmp = tmp[!is.na(tmp)]
    n = length(tmp) - 1
    hold_n[nn] = n +1
    if(loss[nn] == "Quadratic"){hold[nn] = sum(l2_loss(u = offset[nn], a = tmp))/n}
    if(loss[nn] == "Poisson"){hold[nn] = sum(pois_loss(u = offset[nn], a = tmp))/n}
    if(loss[nn] == "Logistic"){hold[nn] = sum(logistic_loss(u = offset[nn], a = tmp))/n}
    if(loss[nn] == "ZTQR"){
       a1 = ifelse(tmp == 0, 1, -1)
       a2 = tmp[tmp != 0]
       lambda = 1
       hold[nn] = (sum(logistic_loss(u = offset[nn], a = a1)) + lambda*sum(l2_loss(u = offset[nn], a = a2)))/n
    }
    if(loss[nn] == "ZTPF"){hold[nn] = sum(hurdle_count_loss(u = offset[nn], a = tmp[tmp>0]))/n}
  }
  if(nreturn){return(cbind(hold, hold_n))}
  if(!nreturn){return(hold)}
}

null_loss = function(A, loss, offset, scale){
  Zmat = matrix(0, ncol=ncol(A), nrow=nrow(A)) 
  return(sum(generic_loss(loss = loss, offset = offset, scale = scale, A = A, XY = Zmat), na.rm = TRUE))
}

generic_grad_x = function(loss, offset, scale, a, x, Y){
  n = length(offset)
  tmp = rep(0,length(x))
  nn_not_na = which(!is.na(a))
  for(nn in nn_not_na){
    if(loss[nn] == "Quadratic"){tmp = tmp + (2/scale[nn])*(Y[,nn]*(x%*%Y[,nn] + offset[nn] - a[nn]))}
    if(loss[nn] == "Poisson"){tmp = tmp + (1/scale[nn])*(Y[,nn]*(exp(x%*%Y[,nn] + offset[nn]) - a[nn]))}
    if(loss[nn] == "Logistic"){tmp = tmp + (-a[nn]/scale[nn])*Y[,nn]/(1 + exp(a[nn]*(x%*%Y[,nn] + offset[nn])))}
    if(loss[nn] == "ZTPF" & !is.na(a[nn])){
      ttmp = x%*%Y[,nn] + offset[nn]
      tmp = tmp + (1/scale[nn])*Y[,nn]*(exp(ttmp)/(1 - exp(-exp(ttmp))) - a[nn])
    }
  }
  return(tmp)
}

ver2_generic_grad_x = function(loss, offset, scale, a, x, Y){
  tmp = rep(0,length(x))
  n = length(offset)
  for(nn in 1:n){
    if(loss[nn] == "Quadratic"){tmp = tmp + (2/scale[nn])*(Y[,nn]*(x%*%Y[,nn] + offset[nn] - a[nn]))}
    if(loss[nn] == "Poisson"){tmp = tmp + (1/scale[nn])*(Y[,nn]*(exp(x%*%Y[,nn] + offset[nn]) - a[nn]))}
    if(loss[nn] == "Logistic"){tmp = tmp + (-a[nn]/scale[nn])*Y[,nn]/(1 + exp(a[nn]*(x%*%Y[,nn] + offset[nn])))}
    if(loss[nn] == "ZTQR"){
      ttmp = x%*%Y[,nn] + offset[nn]
      if(a[nn] == 0){
        tmp = tmp + (-1/scale[nn])*Y[,nn]/(1 + exp(ttmp))
      }else{
        lambda = 1
        tmp = tmp + (1/scale[nn])*Y[,nn]*(1/(1 + exp(-ttmp)) + 2*lambda*(x%*%Y[,nn] + offset[nn] - a[nn]))
      }
    }
    if(loss[nn] == "ZTPF" & !is.na(a[nn])){
      ttmp = x%*%Y[,nn] + offset[nn]
      tmp = tmp + (1/scale[nn])*Y[,nn]*(exp(ttmp)/(1 - exp(-exp(ttmp))) - a[nn])
    }
  }
  return(tmp)
}

generic_grad_Y = function(loss, offset, scale, a, X, y){
  a_na = is.na(a)
  if(sum(a_na) > 0){
    a = matrix(a[!a_na],ncol=1)
    X = X[!a_na,]
  }
  if(loss == "Quadratic"){tmp = (2/scale)*(t(X)%*%(X%*%y + offset - a ))}
  if(loss == "Poisson"){tmp = (1/scale)*(t(X)%*%(exp(X%*%y + offset) - a))}
  if(loss == "Logistic"){tmp = (-1/scale)*t(X)%*%(a/(1+exp((X%*%y + offset)*a)))}
  if(loss == "ZTPF"){
    ttmp = X%*%y + offset
    tmp = (1/scale)*t(X)%*%( exp(ttmp)/(1 - exp(-exp(ttmp)))  - a)
  }

  return(tmp)
}

ver2_generic_grad_y = function(loss, offset, scale, a, X, y){
  if(loss == "Quadratic"){tmp = (2/scale)*(t(X)%*%(X%*%y + offset - a ))}
  if(loss == "Poisson"){tmp = (1/scale)*(t(X)%*%(exp(X%*%y + offset) - a))}
  if(loss == "Logistic"){tmp = (-1/scale)*t(X)%*%(a/(1+exp((X%*%y + offset)*a)))}
  if(loss == "ZTQR"){
    a1 = ifelse(a == 0, 1, -1)
    ttmp = X%*%y + offset
    tmp1 = (-1/scale)*t(X)%*%(a1/(1+exp(ttmp*a1)))
    a2 = (a1 == 1)
    lambda = 1
    tmp2 = lambda*(2/scale)*t(X[a2,])%*%(ttmp[a2]  - a[a2])
    tmp = tmp1 + tmp2 
  }
  if(loss == "ZTPF"){
    ttmp = X%*%y + offset
    tmp = (1/scale)*t(X)%*%( exp(ttmp)/(1 - exp(-exp(ttmp)))  - a)
  }
  return(tmp)
}

predict_x_gd = function(loss,new_a,old_x=NULL,old_y,muj=offset,sigj=scale,gamma=1,small_num=10^-6,alpha=0.25,beta_down=0.7,beta_up=1.05,iter_n=10){
  k = nrow(old_y)
  gamma = gamma/ncol(old_y)
  if(is.null(old_x)){new_x = matrix(1, ncol=k)}
  if(!is.null(old_x)){new_x = old_x}
  old_loss_x = small_num + 1
  new_loss_x = 0
  TT = 0
  while(old_loss_x - new_loss_x > small_num & TT < iter_n){
    old_loss_x = sum(generic_loss(loss=loss, offset=muj, scale=sigj, A=new_a, XY=new_x%*%old_y), na.rm = TRUE)
    grad_x = generic_grad_x(loss=loss, offset=muj, scale=sigj, a=new_a, x=new_x, Y=old_y)
    tmp = new_x - gamma*grad_x
    new_loss_x = sum(generic_loss(loss=loss, offset=muj, scale=sigj, A=new_a, XY=tmp%*%old_y), na.rm = TRUE)
    while((is.finite(new_loss_x)==FALSE | new_loss_x > old_loss_x - alpha*gamma*sum(grad_x^2)) & old_loss_x !=0){
      # line search:
      gamma = beta_down*gamma
      tmp = new_x - gamma*grad_x
      new_loss_x = sum(generic_loss(loss=loss, offset=muj, scale=sigj, A=new_a, XY=tmp%*%old_y), na.rm = TRUE) 
    }
    gamma = beta_up*gamma
    new_x = tmp
    TT = TT + 1
  }
  return(new_x)
}

predict_y_gd = function(loss,new_a,old_x,old_y=NULL,muj,sigj,gamma=1,small_num=10^-6,alpha=0.25,beta_down=0.7,beta_up=1.05,iter_n=10){
  k = ncol(old_x)
  gamma = gamma/nrow(old_x)
  if(is.null(old_y)){new_y = matrix(1, nrow=k)}
  if(!is.null(old_y)){new_y = old_y}
  old_loss_y = small_num + 1
  new_loss_y = 0
  TT = 0
  while(old_loss_y - new_loss_y > small_num & TT < iter_n){
    old_loss_y = sum(generic_loss(loss=loss, offset=muj, scale=sigj, A=new_a, XY=old_x%*%new_y), na.rm = TRUE)
    grad_y = generic_grad_Y(loss=loss, offset=muj, scale=sigj, a=new_a, X=old_x, y=new_y)
    tmp = new_y - gamma*grad_y
    new_loss_y = sum(generic_loss(loss=loss, offset=muj, scale=sigj, A=new_a, XY=old_x%*%tmp), na.rm = TRUE)
    while((is.finite(new_loss_y)==FALSE | new_loss_y > old_loss_y - alpha*gamma*sum(grad_y^2)) & old_loss_y != 0){
      # line search:
      gamma = beta_down*gamma
      tmp = new_y - gamma*grad_y
      new_loss_y = sum(generic_loss(loss=loss, offset=muj, scale=sigj, A=new_a, XY=old_x%*%tmp), na.rm = TRUE)
    }
    TT = TT + 1
    gamma = beta_up*gamma
    new_y = tmp
  }
  return(new_y)
}

predict_x_nm = function(loss,new_a,old_x=NULL,old_y,muj,sigj,gamma_x=0){
 # For null x:
 if(is.null(old_x)){new_x = matrix(1, ncol=k)}
 if(!is.null(old_x)){new_x = old_x}
 # For NAs:
 n = length(muj)
 nn_not_na = which(!is.na(new_a))
 NN = length(nn_not_na)
 if(NN < n){
     muj   = muj[nn_not_na]
     sigj  = sigj[nn_not_na]
     loss  = loss[nn_not_na]
     old_y = old_y[,nn_not_na]
     new_a = new_a[nn_not_na]
 }
 if(gamma_x > 0){
  fn = function(loss,new_a,old_x,old_y,muj,sigj,gamma_x){sum(generic_loss(loss=loss, offset=muj, scale=sigj, A=new_a, XY=old_x%*%old_y)) + gamma_x*sum(old_x^2)}
  gr = function(loss,new_a,old_x,old_y,muj,sigj,gamma_x){ver2_generic_grad_x(loss=loss, offset=muj, scale=sigj, a=new_a, x=old_x, Y=old_y) + 2*gamma_x*old_x}
  tmp = optim(par=old_x, fn=fn, gr=gr, control = list(maxit=500), method="BFGS",loss=loss,new_a=new_a,old_y=old_y,muj=muj,sigj=sigj,gamma_x=gamma_x)$par 
 }
 if(gamma_x == 0){
   fn = function(loss,new_a,old_x,old_y,muj,sigj){sum(generic_loss(loss=loss, offset=muj, scale=sigj, A=new_a, XY=old_x%*%old_y))}
   gr = function(loss,new_a,old_x,old_y,muj,sigj){ver2_generic_grad_x(loss=loss, offset=muj, scale=sigj, a=new_a, x=old_x, Y=old_y)}
   tmp = optim(par=old_x, fn=fn, gr=gr, control = list(maxit=500), method="BFGS",loss=loss,new_a=new_a,old_y=old_y,muj=muj,sigj=sigj)$par
 }
 return(tmp)
}

predict_x_nm_k1 = function(loss,new_a,old_x=NULL,old_y,muj,sigj){
  # For null x:
  if(is.null(old_x)){new_x = matrix(1, ncol=k)}
  if(!is.null(old_x)){new_x = old_x}
  # For NAs:
  n = length(muj)
  nn_not_na = which(!is.na(new_a))
  NN = length(nn_not_na)
  if(NN < n){
    muj   = muj[nn_not_na]
    sigj  = sigj[nn_not_na]
    loss  = loss[nn_not_na]
    old_y = old_y[,nn_not_na]
    new_a = new_a[nn_not_na]
  } 
  fn = function(loss,new_a,old_x,old_y,muj,sigj){sum(generic_loss(loss=loss, offset=muj, scale=sigj, A=new_a, XY=old_x%*%old_y))}
  gr = function(loss,new_a,old_x,old_y,muj,sigj){ver2_generic_grad_x(loss=loss, offset=muj, scale=sigj, a=new_a, x=old_x, Y=matrix(old_y,nrow=1))}
  tmp = optim(par=old_x, fn=fn, gr=gr, control = list(maxit=500), method="BFGS",loss=loss,new_a=new_a,old_y=matrix(old_y,nrow=1),muj=muj,sigj=sigj)$par 
  return(tmp)
}

    
predict_y_nm = function(loss,new_a,old_x,old_y,muj,sigj,gamma_y=0){
  # For NAs:
  a_na = is.na(new_a)
  if(sum(a_na) > 0){
    new_a = matrix(new_a[!a_na],ncol=1)
    old_x = old_x[!a_na,]
  }
  if(gamma_y > 0){
    fn = function(loss,new_a,old_x,old_y,muj,sigj,gamma_y){sum(generic_loss(loss=loss, offset=muj, scale=sigj, A=new_a, XY=old_x%*%old_y)) + gamma_y*sum(old_y^2)}
    gr = function(loss,new_a,old_x,old_y,muj,sigj,gamma_y){ver2_generic_grad_y(loss=loss, offset=muj, scale=sigj, a=new_a, X=old_x, y=old_y) + 2*gamma_y*old_y}
    tmp = optim(par=old_y, fn=fn, gr=gr, control = list(maxit=500), method="BFGS", loss=loss,new_a=new_a,old_x=old_x,muj=muj,sigj=sigj,gamma_y=gamma_y)$par
  }
  if(gamma_y == 0){
    fn = function(loss,new_a,old_x,old_y,muj,sigj){sum(generic_loss(loss=loss, offset=muj, scale=sigj, A=new_a, XY=old_x%*%old_y))}
    gr = function(loss,new_a,old_x,old_y,muj,sigj){ver2_generic_grad_y(loss=loss, offset=muj, scale=sigj, a=new_a, X=old_x, y=old_y)}
    tmp = optim(par=old_y, fn=fn, gr=gr, control = list(maxit=500), method="BFGS", loss=loss,new_a=new_a,old_x=old_x,muj=muj,sigj=sigj)$par
  }
  return(tmp)
}

predict_y_nm_k1 = function(loss,new_a,old_x,old_y,muj,sigj){
  # For NAs:
  a_na = is.na(new_a)
  if(sum(a_na) > 0){
    new_a = matrix(new_a[!a_na],ncol=1)
    old_x = old_x[!a_na,]
  }
  fn = function(loss,new_a,old_x,old_y,muj,sigj){sum(generic_loss(loss=loss, offset=muj, scale=sigj, A=new_a, XY=matrix(old_x,ncol=1)%*%matrix(old_y,nrow=1)))}
  gr = function(loss,new_a,old_x,old_y,muj,sigj){ver2_generic_grad_y(loss=loss, offset=muj, scale=sigj, a=new_a, X=matrix(old_x,ncol=1), y=matrix(old_y,nrow=1))}
  tmp = optim(par=old_y, fn=fn, gr=gr, control = list(maxit=500), method="BFGS", loss=loss,new_a=new_a,old_x=matrix(old_x,ncol=1),muj=muj,sigj=sigj)$par
  return(tmp)
}

XY_pca_transform = function(X,Y){
  s = svd(X)
  newX = s$u %*% diag(s$d^0.5)
  newY = diag(s$d^0.5) %*% t(s$v) %*% Y
  return(list(newX,newY))
} 

hurdle_count_reconstruct = function(logistic_var, pois_var, logistic_offset, pois_offset, scale1, scale2, nu = 0){
  loss_logit_pos1 = (1/scale1) * logistic_loss(u = logistic_var + logistic_offset, a = 1)
  loss_logit_neg1 = (1/scale1) * logistic_loss(u = logistic_var + logistic_offset, a = -1)
  pois_min = exp(pois_var + pois_offset)
  # since pois_loss(u = x, a = e^x) = 0:
  loss_pois_min   = 0 
  tmp = (loss_logit_neg1 + loss_pois_min)/loss_logit_pos1
  my_hat = rep(nu, length(logistic_var))
  my_hat[tmp < 1] = pois_min[tmp < 1] 
  return(my_hat)
}

my_reconstruct = function(loss, X, Y, offset){
  if(loss == "Quadratic"){return(scale(X%*%Y, center=-offset, scale=FALSE))}
  if(loss == "Poisson"){return(exp(scale(X%*%Y, center=-offset, scale=FALSE)))}
  if(loss == "Logistic"){return(1/(1+exp(-scale(X%*%Y, center=-offset, scale=FALSE))))}
  if(loss == "ZTQR"){
    ttmp = scale(X%*%Y, center=-offset, scale=FALSE)
    return(ttmp)
  }
  if(loss == "ZTPF"){
    ttmp = exp(scale(X%*%Y, center=-offset, scale=FALSE))
    return(sapply(ttmp, pred_hc_a))
  }
}

# --------------------------------- #
#    Specific Loss Functions        #
# --------------------------------- #
l2_loss = function(u, a){return( (a - u)^2 ) } 
logistic_loss = function(u,a){
  tmp = -a*u
  set1 = which(tmp < 20)
  tmp[set1] = log(1 + exp(tmp[set1]))
  return(tmp)
}
pois_loss = function(u,a){
  tmp = a*log(a)
  tmp[a==0] = 0
  return(exp(u) - u*a - a + tmp)
}
partial_pois_loss = function(u,a){return(exp(u) - u*a)}
### Hurdle count look up values when a is small (1<= a <=14)
Ahh = 17
mle_u = 1:Ahh
for(aa in 1:Ahh){
  u0 = aa
  test = 1
  while(abs(test) > 0.0000001){
    top = aa - u0 - aa*exp(-u0)
    bot = aa - u0 - 1
    un = u0 - (top/bot)
    test = un - u0
    u0 = un
  }
  mle_u[aa] = u0
}
lookup_loss_hc = log((mle_u^(1:Ahh)) / (exp(mle_u) - 1))
lookup_loss_hc[1] = 0
### Need loss calculation function:
hurdle_count_loss = function(u, a, look_up_a = lookup_loss_hc){
  norm_const = l_metric = 1:length(a)
  seta = (a > 17)
  norm_const[!seta] = look_up_a[a[!seta]]
  norm_const[seta] = a[seta]*log(a[seta]) - a[seta]
  setu = (u > 17)
  l_metric[setu] = a[setu]*u[setu] - exp(u[setu])
  l_metric[!setu] = a[!setu]*u[!setu] - log(exp(exp(u[!setu])) - 1)
  return(round(norm_const - l_metric,digits = 7))  
} 
### Need function to predict a based on fitted u #####
pred_hc_a = function(x){
  if(x > 0 & x <= 17){ return(x/(1 - exp(-x)))}
  if(x > 17){return(x)}
  if(x == 0){return(1)}
}


### ZIFA optimization functions
optimize_lambda = function(lambda, X, H){
  fn = function(lambda, X, H){
    p0 = exp(-lambda*(X^2))
    return(sum(lambda*(X^2)[H == 1]) - sum(log(1-p0[H==0])))
  }
  gr = function(lambda, X, H){
    part1 = sum((X^2)[H == 1])
    p0 = exp(-lambda*(X^2))
    part2= -1*sum(((X^2)*(p0/(1-p0)))[H==0])
    return(part1 + part2)
  }
  return(optim(par=lambda, f=fn, gr=NULL,control = list(maxit=500),lower=0,upper=50, method="Brent",X=X,H=H)$par)
}
optimize_mu = function(mu, lambda, A, H){
  fn = function(mu, lambda, A, H, X){
    part1 = lambda*(mu^2)*sum(H == 1) - sum(H == 0)*log(1 - exp(-lambda*(mu^2)))
    part2 = sum((A[H==0] - mu)^2)
    return(part1 + part2)
  }
  return(optim(par=mu, fn=fn, gr=NULL,control = list(maxit=500),lower=0,upper=100,method="Brent",lambda=lambda,A=A,H=H)$par)
}
optimize_pair = function(mu, lambda, A, H, tol = 0.1, maxiter = 100){
  cur_tol = 1+tol
  my_mu = mu
  my_lambda = lambda
  tt = 1
  while(cur_tol > tol | tt > maxiter){
    my_lambda_t = optimize_lambda(lambda=my_lambda, X=rep(my_mu,length(H)), H=H)
    my_mu_t = optimize_mu(mu=my_mu, lambda = my_lambda, A=A, H=H)
    cur_tol = max( c(abs(my_lambda - my_lambda_t) , abs(my_mu - my_mu_t)) )
    my_mu = my_mu_t
    my_lambda = my_lambda_t
    tt = tt + 1
  }
  return(c(my_mu, my_lambda))
}
zifa_column_loss = function(A,H,X,offset,scale,decay){
  Z = X + offset
  p0 = exp(-decay*(Z^2))
  binom_part = sum(decay*(Z^2)[H == 1]) - sum(log(1-p0[H==0]))
  quad_part  = sum((A[H==0] - Z[H==0])^2)
  return( (1/scale)*(binom_part + quad_part) )
}

zifa_grad_x = function(offset, scale, decay, a, x, Y){
  ### assumes no NAs:
  tmp = rep(0,length(x))
  p = length(offset)
  for(pp in 1:p){
    zij = x%*%Y[,pp] + offset[pp]
    hold = 2*decay[pp]*zij*Y[,pp]
    if(a[pp] == 0){tmp_new = hold/scale[pp]}
    if(a[pp] != 0){
      p0 = exp(-decay[pp]*(zij^2))
      tmp_new =  (1/scale[pp])*( -(p0/(1-p0))*hold  +   2*(Y[,pp]*(zij - a[pp])))
    }
    tmp = tmp + tmp_new
  }
  return(tmp)
}
zifa_grad_y = function(offset, scale, decay, a, X, y){
  a_index_0 = which(a == 0)
  a_index_1 = which(a != 0)
  Z = X%*%y + offset
  tmp_0 = (2/scale)*decay*(t(X[a_index_0,])%*%(Z[a_index_0,]))
  tmp_1 = (2/scale)*decay*(t(X[a_index_1,])%*%(Z[a_index_1,]))
  p0 = exp(-decay*Z^2)
  p_frac = p0/(1-p0)
  tmp_quad  = (2/scale)*(t(X[a_index_1,])%*%(p_frac[a_index_1]*(Z[a_index_1] - a[a_index_1])))
  return(tmp_0 + tmp_1 + tmp_quad)
}
zifa_predict_x_nm = function(new_a,h,old_x,old_y,muj,sigj,decay){
  fn = function(new_a,h,old_x,old_y,muj,sigj,decay){
    tmp_loss = 0
    zi = old_x%*%old_y
    for(pp in 1:length(old_x)){tmp_loss = tmp_loss + zifa_column_loss(A=new_a[pp],H=h[pp],X = zi[pp],offset=muj[pp],scale=sigj[pp],decay=decay[pp])}
    return(tmp_loss)
  }
  gr = function(new_a,h,old_x,old_y,muj,sigj,decay){
    zifa_grad_x(offset=muj, scale=sigj, decay=decay, a=new_a, x=old_x, Y=old_y)
  }
  tmp = optim(par=old_x, fn=fn, gr=gr, control = list(maxit=500), method="BFGS",new_a=new_a,h=h,old_y=old_y,muj=muj,sigj=sigj,decay=decay)$par
  return(tmp)
}
zifa_predict_y_nm = function(new_a,h,old_x,old_y,muj,sigj,decay){
  fn = function(new_a,h,old_x,old_y,muj,sigj,decay){
    Z = old_x%*%old_y
    return(zifa_column_loss(A=new_a,H=h,X = Z,offset=muj,scale=sigj,decay=decay))
  }
  gr = function(new_a,h,old_x,old_y,muj,sigj,decay){zifa_grad_y(offset=muj, scale=sigj, decay=decay, a=new_a, X=old_x, y=old_y)}
  tmp = optim(par=old_y, fn=fn, gr=gr, control = list(maxit=500), method="BFGS",new_a=new_a,h=h,old_x=old_x,muj=muj,sigj=sigj,decay=decay)$par
  return(tmp)
}

