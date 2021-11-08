# Replication codes for 
# "Using generalized estimating equations to estimate nonlinear models with spatial data"
# by Weining Wang, Jeffrey M. Wooldridge, Mengshan Xu, Cuicui Lu

# Authors of the codes: Mengshan Xu, Weining Wang
# We thank Dr. Xuening Zhu for her participation and insightful discussions.

rm(list=ls())
rm(list=ls(all=TRUE))

# The following library is required
library(mvtnorm)
library(CovTools)
library(Matrix)
library(corrplot)
library(MASS)

#sample size
N=n=400
#Group number
G=100
# Number of Parameters
m=2
# Number of Monte-Carlo simulations
NumIt=10
# The true value of Beta
beta_0=c(1,1)
# Spatial correlation parameter
rho_0=0.1
#The type of the true spatial correlation matrix: 0 identity matrix, 1, Grouped
class=1
#N_sim_w=1000

member = rep(1:G, each = N/G)

# Data matrix to collect the result.
estimates_PMLE<- matrix(0, nrow= NumIt, ncol= m)
estimates_P_E_IV<- matrix(0, nrow= NumIt, ncol= m)
estimates_P_IV<- matrix(0, nrow= NumIt, ncol= m)


setwd("C:\\Users\\Uwe-x\\Dropbox\\NonlinerGEE\\code\\Uploaded version") 
sim = 10
set.seed(sim)# seed for the distance (closeness) matrix. It will be generated in the following source file
# source files needed 
source("simulatorW2.R")
source("generating_GEE_spatial.R")
source("functions_GEE_spatial.R") 

# closeness matrix (1-d)
Sig=Sig_g(1,rho_0) 

#vw=wv_uniform(1)
#vw=wv_log_normal(1)
vw=wv_normal(class)
#vw
v=vw$V
#vvvv=vw$W
# mean(v)
#var(v)
W=vw$W

starter_sse = proc.time()
time_index=starter_sse
# data generating process
for (j in 1: NumIt){
  #j=1
  sim = 100+j
  set.seed(sim)
  # Sig=Sig_g(0,rho_0)
  # vw=wv_uniform(1)
  X = matrix(rnorm(m*n),n,m, byrow = FALSE)
  y0 = X%*%beta_0+v
 
  Y = (y0>0)
  # mean(Y)
  # Y= y0
# initial values
init=c(0.8,0.5)

# QMLE estimator

#p_mle <- optim(init, PMLE)
p_mle <- optim(init, PMLE_b)
estimates_PMLE[j,]=p_mle$par

# p_iv  <- optim(p_mle$par, P_IV)
# estimates_P_IV[j,]=p_iv$par

#Working variance-covariance matrix

#W_in=w_inv_logN_E(p_mle$par)
W_in=w_inv_binary_E(p_mle$par)

#simulated weighting matrix.(just for test)
#W_in=w_inv_binary_sim(p_mle$par,N_sim_w)

p_e_iv  <- optim(p_mle$par, P_E_IV_bi)
estimates_P_E_IV[j,]=p_e_iv$par
}
time_sse = (proc.time() -starter_sse)[3]
ones=matrix(1,NumIt,1, byrow = FALSE)

Pmle_sd=sqrt(colMeans((estimates_PMLE-ones%*%beta_0)^2))
GEE_sd=sqrt(colMeans((estimates_P_E_IV-ones%*%beta_0)^2))

z_Pmle=estimates_PMLE/Pmle_sd
z_GEE=estimates_P_E_IV/GEE_sd

q_Pmle1=quantile(z_Pmle[,1],probs=c(0.25, 0.75))
q_Pmle2=quantile(z_Pmle[,2],probs=c(0.25, 0.75))
q_GEE1=quantile(z_GEE[,1],probs=c(0.25, 0.75))
q_GEE2=quantile(z_GEE[,2],probs=c(0.25, 0.75))

p_Pmle1=1-pnorm(q_Pmle1)
p_Pmle2=1-pnorm(q_Pmle2)
p_GEE1=1-pnorm(q_GEE1)
p_GEE2=1-pnorm(q_GEE2)


rho_0
noquote(format(colMeans(estimates_PMLE), digits=1, nsmall=4))
noquote(format(colMeans(estimates_P_E_IV), digits=1, nsmall=4))
#sd
noquote(format(Pmle_sd,digits=2, nsmall=4))
noquote(format(GEE_sd,digits=2, nsmall=4))
#MSE
format((colMeans((estimates_PMLE-ones%*%beta_0)^2)*n),digits=2, nsmall=4)
format((colMeans((estimates_P_E_IV-ones%*%beta_0)^2)*n),digits=2, nsmall=4)


time_sse