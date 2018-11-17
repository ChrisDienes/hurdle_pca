
# This file generates the 30 data sets for the missing data analysis


#### Edit the below file_root to be the location of the saved folder:
file_root = "C:/........./paper_r_examples/"

# Create data sets 2-30:
# Data set 1 was created previously, with MAR probability parameter 
# rounded to c = 0.40 (root solve would be 0.3855668)

library(rootSolve)
source(file = paste0(file_root, "glrm_functions.R"))

# ------------------------------------------------- #
# Simulate data:
# ------------------------------------------------- #
# Following similar format to ZIFA paper:
set.seed(23)
m = 5000
n = 10
k = 4
my_c  = c(0.4, rep(NA,29))
my_a2 = c(2.031158, rep(NA,29))
my_a3 = c(2.993387, rep(NA,29))  
for(dd in 2:30){
  A_mean  = matrix(rnorm(n*k, 0, 1), nrow = n, ncol=k)
  W = matrix(0, ncol = n, nrow = n)
  sigma_W = 1
  diag(W) = sigma_W*runif(n, 0.9, 1.1)
  my_mu = 1:n
  my_data = matrix(0, nrow = m, ncol = n)
  for(ii in 1:m){
    zi = matrix(rnorm(k,0,1), ncol=1)
    wi = matrix(rnorm(n,0,diag(W)), ncol=1)
    my_data[ii,] =  A_mean%*%zi + my_mu + wi
  }
  theory_cov = A_mean%*%t(A_mean) + W
  theory_cor = (A_mean%*%t(A_mean) + W)*(matrix(diag(theory_cov)^-.5,ncol=1) %*% matrix(diag(theory_cov)^-.5,nrow=1))
  #cor(my_data)
  
  # Missing completely at random (MCAR) data
  MCAR = list()
  A = my_data
  c = log(1/0.1537191 - 1)
  prob = 1/(1 + exp(1.7))
  MCAR$mean_prob = prob
  MCAR$my_na = rbinom(n=m, size=1, prob=prob)
  MCAR$obs_mean_prob = mean(MCAR$my_na)
  MCAR$keep_for_later = A[MCAR$my_na == 1,1]
  A[MCAR$my_na == 1,1] = NA 
  MCAR$A = A

  # Missing at random (MAR) data
  MAR = list()
  A = my_data
  f = function(c){1/(1 + exp(1.7)) - mean(1/(1 + exp(c + A[,2] + A[,3])))}
  c = uniroot(f, c(-10, 10))$root
  my_c[dd] = c
  prob = 1/(1 + exp(c + A[,2] + A[,3]))
  print(mean(prob))
  my_a2[dd] = mean(A[,2])
  my_a3[dd] = mean(A[,3])
  MAR$mean_prob = mean(prob)
  MAR$my_na = rbinom(n=m, size=1, prob=prob)
  #summary(prob[MAR$my_na == 0])
  #summary(prob[MAR$my_na == 1])
  MAR$obs_mean_prob = mean(MAR$my_na)
  MAR$keep_for_later = A[MAR$my_na == 1,1]
  A[MAR$my_na == 1,1] = NA 
  MAR$A = A

 # write.csv(MCAR$A, file=paste0(file_root,"missing_data/simulated_data/MCAR_data_",dd,".csv"), row.names=FALSE)
 # write.csv(MCAR$keep_for_later, file=paste0(file_root,"missing_data/simulated_data/MCAR_truth_",dd,".csv"), row.names=FALSE)
 # write.csv(MAR$A, file=paste0(file_root,"missing_data/simulated_data/MAR_data_",dd,".csv"), row.names=FALSE)
 # write.csv(MAR$keep_for_later, file=paste0(file_root,"missing_data/simulated_data/MAR_truth_",dd,".csv"), row.names=FALSE)
}

# my_c
# my_ave_lin_pred = my_c + my_a2 + my_a3
# my_ave_lin_pred

