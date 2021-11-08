

# functions


PMLE <- function(x_ini) {  
  beta=x_ini
  temp=t(X)%*%(Y-exp(X%*%(beta)))
  #norm(matrix(3,2,1, byrow = FALSE),"2")
  norm(temp,"2")
}
# the inverse function
w_inv_logN_0 <- function(w,beta) {  
  #beta=beta_0
  #w=diag(1/as.vector(exp(X%*%(beta))^2))
  #w=W
  W1= matrix(0,N,N)
  W2= matrix(0,N,N)
  W0=w
  for    (group in 1:G)
  {
    nG1= ((group-1)*nss+1)
    nG2= ((group-1)*nss+nss)
    
    #C=exp(1)*(exp(1)-1)
    C=exp(1)-1
    # lets do it first then check others.
    W1[nG1,nG1]=C*exp(X[nG1,]%*%(beta))^2+exp(X[nG1,]%*%(beta))
    W1[nG2,nG2]=C*exp(X[nG2,]%*%(beta))^2+exp(X[nG2,]%*%(beta))
    W1[nG1,nG2]=W0[nG1,nG2]*exp(X[nG2,]%*%(beta))*exp(X[nG1,]%*%(beta))
    W1[nG2,nG1]=W0[nG1,nG2]*exp(X[nG2,]%*%(beta))*exp(X[nG1,]%*%(beta))
    W1[nG1,nG2]=0
    W1[nG2,nG1]=0
    W2[nG1:nG2,nG1:nG2]= solve( W1[nG1:nG2,nG1:nG2])
    # W[nG1:nG2,nG1:nG2]=  xxx[nG1:nG2,nG1:nG2]*vargamma.true[nG1:nG2,nG1:nG2]
  }
  #KK=W2%*%W1
  W2
}

# this is the case we only estimate variance
w_inv_logN_E0 <- function(w,beta) {  
  #beta=beta_0
  #w=diag(1/as.vector(exp(X%*%(beta))^2))
  #w=W
  W1= matrix(0,N,N)
  W2= matrix(0,N,N)
  W0=w
  
  # the calculation of C
  yc=exp(X%*%beta)
  u_s=Y-yc
  # the simplst case, it is good
  y=u_s*u_s-exp(X%*%beta)
  x=exp(2*X%*%beta)
  
  # adjust for hetro-case
  y=y/exp(X%*%beta)
  x=x/exp(X%*%beta)
  
  aaa=lm(y~x+0)
  C=aaa$coefficients
  #C=exp(3)-1
  for    (group in 1:G)
  {
    nG1= ((group-1)*nss+1)
    nG2= ((group-1)*nss+nss)
    #C=exp(1)*(exp(1)-1)
    #C=exp(1)-1
    # lets do it first then check others.
    W1[nG1,nG1]=C*exp(X[nG1,]%*%(beta))^2+exp(X[nG1,]%*%(beta))
    W1[nG2,nG2]=C*exp(X[nG2,]%*%(beta))^2+exp(X[nG2,]%*%(beta))
    W1[nG1,nG2]=C*W0[nG1,nG2]*exp(X[nG2,]%*%(beta))*exp(X[nG1,]%*%(beta))
    W1[nG2,nG1]=C*W0[nG1,nG2]*exp(X[nG2,]%*%(beta))*exp(X[nG1,]%*%(beta))
    # W1[nG1,nG2]=0
    # W1[nG2,nG1]=0
    W2[nG1:nG2,nG1:nG2]= solve( W1[nG1:nG2,nG1:nG2])
    # W[nG1:nG2,nG1:nG2]=  xxx[nG1:nG2,nG1:nG2]*vargamma.true[nG1:nG2,nG1:nG2]
  }
  #KK=W2%*%W1
  W2
}

# estimated inverse function.
w_inv_unif_E0 <- function(w,beta) {  
  #beta=beta_0
  #w=diag(1/as.vector(exp(X%*%(beta))^2))
  #w=W
  W1= matrix(0,N,N)
  W2= matrix(0,N,N)
  W0=w
  
  # the calculation of C
  yc=exp(X%*%beta)
  u_s=Y-yc
  # the simplst case, it is good
  y=u_s*u_s-exp(X%*%beta)
  x=exp(2*X%*%beta)
  
  # adjust for hetro-case
  y=y/exp(X%*%beta)
  x=x/exp(X%*%beta)
  
  aaa=lm(y~x+0)
  C=aaa$coefficients
 
  # end of the calculation of rho
  
  for    (group in 1:G)
  {
    nG1= ((group-1)*nss+1)
    nG2= ((group-1)*nss+nss)
    # Now I need to estimate this C.
    # lets do it first then check others.
    W1[nG1,nG1]=C*exp(X[nG1,]%*%(beta))^2+exp(X[nG1,]%*%(beta))
    W1[nG2,nG2]=C*exp(X[nG2,]%*%(beta))^2+exp(X[nG2,]%*%(beta))
    
    # the following works with numIt=200, note the following is still not the estimate, only C
    # is estimated.
     W1[nG1,nG2]=W0[nG1,nG2]*C*exp(X[nG2,]%*%(beta))*exp(X[nG1,]%*%(beta))
     W1[nG2,nG1]=W0[nG1,nG2]*C*exp(X[nG2,]%*%(beta))*exp(X[nG1,]%*%(beta))
    
    # W1[nG1,nG2]=0
    # W1[nG2,nG1]=0
    
    W2[nG1:nG2,nG1:nG2]= solve( W1[nG1:nG2,nG1:nG2])# why should they multiply it
    # W[nG1:nG2,nG1:nG2]=  xxx[nG1:nG2,nG1:nG2]*vargamma.true[nG1:nG2,nG1:nG2]
  }
  #KK=W2%*%W1
  W2
}


w_inv_logN_E <- function(beta) {  
  #beta=beta_0

  W1= matrix(0,N,N)
  W2= matrix(0,N,N)
  #W0=w
  
  # the calculation of C
  yc=exp(X%*%beta)
  u_s=Y-yc
  # the simplst case, it is good
  y=u_s*u_s-exp(X%*%beta)
  x=exp(2*X%*%beta)
  
  # adjust for hetro-case
  y=y/exp(X%*%beta)
  x=x/exp(X%*%beta)
  
  aaa=lm(y~x+0)
  C=aaa$coefficients
  
  u_sc=u_s/exp(X%*%beta)
  #u_sc=v-1

  e_rho <- optimize(rho_E, c(0,1),tol = 0.0001,ux=u_sc,C_rho=C) # rho should be constraint into 0 and 1
  rho_e=e_rho$minimum
  #rho_e=rho_0
  #C=exp(3)-1
  for    (group in 1:G)
  {
    #group=1
    #beta=beta_0
    #C=1/3
    #rho_e=rho_0
    nG1= ((group-1)*nss+1)
    nG2= ((group-1)*nss+2)
    nG3= ((group-1)*nss+3)
    nG4= ((group-1)*nss+4)
    # C=exp(1)*(exp(1)-1)
    # C=exp(1)-1
    # lets do it first then check others.
    W1[nG1,nG1]=C*exp(X[nG1,]%*%(beta))^2+exp(X[nG1,]%*%(beta))
    W1[nG2,nG2]=C*exp(X[nG2,]%*%(beta))^2+exp(X[nG2,]%*%(beta))
    W1[nG3,nG3]=C*exp(X[nG3,]%*%(beta))^2+exp(X[nG3,]%*%(beta))
    W1[nG4,nG4]=C*exp(X[nG4,]%*%(beta))^2+exp(X[nG4,]%*%(beta))
    
  
    W1[nG1,nG2]=W1[nG2,nG1]=C*rho_e*(1-loc_diff_abs[nG1,nG2])*exp(X[nG2,]%*%(beta))*exp(X[nG1,]%*%(beta))
    W1[nG1,nG3]=W1[nG3,nG1]=C*rho_e*(1-loc_diff_abs[nG1,nG3])*exp(X[nG3,]%*%(beta))*exp(X[nG1,]%*%(beta))
    W1[nG1,nG4]=W1[nG4,nG1]=C*rho_e*(1-loc_diff_abs[nG1,nG4])*exp(X[nG4,]%*%(beta))*exp(X[nG1,]%*%(beta))

    W1[nG2,nG3]=W1[nG3,nG2]=C*rho_e*(1-loc_diff_abs[nG2,nG3])*exp(X[nG2,]%*%(beta))*exp(X[nG3,]%*%(beta))
    W1[nG2,nG4]=W1[nG4,nG2]=C*rho_e*(1-loc_diff_abs[nG2,nG4])*exp(X[nG2,]%*%(beta))*exp(X[nG4,]%*%(beta))

    W1[nG3,nG4]=W1[nG4,nG3]=C*rho_e*(1-loc_diff_abs[nG3,nG4])*exp(X[nG3,]%*%(beta))*exp(X[nG4,]%*%(beta))
    #
    # W1[nG1,nG2]=0
    # W1[nG2,nG1]=0
    W2[nG1:nG4,nG1:nG4]= solve( W1[nG1:nG4,nG1:nG4])# why should they multiply it
    # W[nG1:nG2,nG1:nG2]=  xxx[nG1:nG2,nG1:nG2]*vargamma.true[nG1:nG2,nG1:nG2]
    }
  
  #KK=W2%*%W1
  W2
}

# uniform with estiamted version.
w_inv_unif_E <- function(beta) {  
  #beta=p_iv$par
  W1= matrix(0,N,N)
  W2= matrix(0,N,N)

  #W0=w
  
  # the calculation of C
  yc=exp(X%*%beta)
  u_s=Y-yc
  # the simplst case, it is good
  y=u_s*u_s-exp(X%*%beta)
  x=exp(2*X%*%beta)
  
  # adjust for hetro-case
  y=y/exp(X%*%beta)
  x=x/exp(X%*%beta)
  
  aaa=lm(y~x+0)
  C=aaa$coefficients
  #C=1/3
  # the calculation of rho.
  u_sc=u_s/exp(X%*%beta)
  #u_sc=v-1
  # u_jc=u_sc[2*(1:G)-1]
  # u_ic=u_sc[2*(1:G)]
  e_rho <- optimize(rho_E, c(0,1),tol = 0.0001,ux=u_sc,C_rho=C)
  rho_e=e_rho$minimum
  # end of the calculation of rho
  
  for    (group in 1:G)
  {
    nG1= ((group-1)*nss+1)
    nG2= ((group-1)*nss+nss)
    # Now I need to estimate this C.
    
    # lets do it first then check others.
    W1[nG1,nG1]=C*exp(X[nG1,]%*%(beta))^2+exp(X[nG1,]%*%(beta))
    W1[nG2,nG2]=C*exp(X[nG2,]%*%(beta))^2+exp(X[nG2,]%*%(beta))
    
    # the following works with numIt=200, note the following is still not the estimate, only C
    # is estimated.
    W1[nG1,nG2]= rho_e*C*exp(X[nG2,]%*%(beta))*exp(X[nG1,]%*%(beta))
    W1[nG2,nG1]= rho_e*C*exp(X[nG2,]%*%(beta))*exp(X[nG1,]%*%(beta))
    # W1[nG1,nG2]=W0[nG1,nG2]*C*exp(X[nG2,]%*%(beta))*exp(X[nG1,]%*%(beta))
    # W1[nG2,nG1]=W0[nG1,nG2]*C*exp(X[nG2,]%*%(beta))*exp(X[nG1,]%*%(beta))
    
    # W1[nG1,nG2]=0
    # W1[nG2,nG1]=0
    # 
    W2[nG1:nG2,nG1:nG2]= solve( W1[nG1:nG2,nG1:nG2])# why should they multiply it
    # W[nG1:nG2,nG1:nG2]=  xxx[nG1:nG2,nG1:nG2]*vargamma.true[nG1:nG2,nG1:nG2]
  }
  #KK=W2%*%W1
  W2
}

# tau_e <- function(beta) { 
#   # beta=beta_0
#   u_s=(Y-exp(X%*%(beta)))*(Y-exp(X%*%(beta)))
#   y=u_s-exp(X%*%(beta))
#   x=exp(2*X%*%(beta))
#   lm(y~x+0)
#   exp(1)*(exp(1)-1)
# }

# This is only uniform case with true values 1/3
w_inv_unif <- function(w,beta) {  
  #beta=beta_0
  #w=diag(1/as.vector(exp(X%*%(beta))^2))
  #w=W
  W1= matrix(0,N,N)
  W2= matrix(0,N,N)
  W0=w
  for    (group in 1:G)
  {
    nG1= ((group-1)*nss+1)
    nG2= ((group-1)*nss+nss)
    
    #C=exp(1)*(exp(1)-1)
    C=1/3
    # lets do it first then check others.
    W1[nG1,nG1]=C*exp(X[nG1,]%*%(beta))^2+exp(X[nG1,]%*%(beta))
    W1[nG2,nG2]=C*exp(X[nG2,]%*%(beta))^2+exp(X[nG2,]%*%(beta))
    # W1[nG1,nG2]=exp(1)*(exp(W0[nG1,nG2])-1)*exp(X[nG2,]%*%(beta))*exp(X[nG1,]%*%(beta))
    # W1[nG2,nG1]=exp(1)*(exp(W0[nG2,nG1])-1)*exp(X[nG2,]%*%(beta))*exp(X[nG1,]%*%(beta))
    W1[nG1,nG2]=0
    W1[nG2,nG1]=0
    W2[nG1:nG2,nG1:nG2]= solve( W1[nG1:nG2,nG1:nG2])# why should they multiply it
    # W[nG1:nG2,nG1:nG2]=  xxx[nG1:nG2,nG1:nG2]*vargamma.true[nG1:nG2,nG1:nG2]
  }
  #KK=W2%*%W1
  W2
}

# this is the function only for simulating v.
rho_v <- function(init) {  
  
  rho=init
  #rho=0.9
  d=matrix(0, nrow= G, ncol= 1)
  for (dd in 1:G){
    d[dd]=loc_diff_abs[(2*dd-1),2*dd]
  }
  C=1/3
  
  temp=(u_jc)*(u_ic)-C*rho*(1-d) 
  
  norm(temp,"2")
}

rho_E <- function(init,ux,C_rho) {  
  #rho=p_iv$par
  rho=init
  u_sc=ux
  C=C_rho
  #C=1/3
  u_jc=u_sc[2*(1:G)-1]
  u_ic=u_sc[2*(1:G)]
  
  
  d=matrix(0, nrow= G, ncol= 1)
  for (dd in 1:G){
    d[dd]=loc_diff_abs[(2*dd-1),2*dd]
  }
  # C=1/3
  # FOC target function, see sketch it seem 
  # temp=(u_jc*u_ic-C*exp(-rho*d))*exp(rho*d)*d
  # temp=(u_jc*u_ic-C*exp(-rho*d))*d
  #temp=(u_jc)*(u_ic)-C*rho 
  temp=(u_jc)*(u_ic)-C*rho*(1-d) 
  
  norm(temp,"2")
}

P_E_IV <- function(x_ini) {  #very very slow.
  beta=x_ini
  #W_in=w_inv(W,beta)
  #W=diag(1/as.vector(exp(X%*%(beta))^2))
  #temp=t(X)%*%W_in%*%(exp(X%*%(beta))*(Y-exp(X%*%(beta)))) # this is wrong, see sketch
  temp=t(diag(as.vector(exp(X%*%(beta))))%*%X)%*%W_in%*%(Y-exp(X%*%(beta)))
  # kk=t(X)
  # kkk=exp(X%*%(beta))
  #norm(matrix(3,2,1, byrow = FALSE),"2")
  norm(temp,"2")
}


P_IV <- function(x_ini) {  
  beta=x_ini
  temp=t(X)%*%((Y-exp(X%*%(beta)))*exp(-1*X%*%(beta)))
  #norm(matrix(3,2,1, byrow = FALSE),"2")
  norm(temp,"2")
}

P_IV0 <- function(x_ini) {  
  beta=x_ini
  temp=t(X)%*%((Y-exp(X%*%(beta))*exp(0.5))*exp(-1*X%*%(beta)))
  #norm(matrix(3,2,1, byrow = FALSE),"2")
  norm(temp,"2")
}
####################################################################################
# next is binary PLME

PMLE_b <- function(x_ini) {  
  #beta=beta_0
  beta=x_ini
  #temp=t(X)%*%(Y-exp(X%*%(beta)))
  temp=Y*log( pnorm(X%*%beta) )+(1-Y)*log(1-pnorm(X%*%beta))
  # norm(matrix(3,2,1, byrow = FALSE),"2")
  -sum(temp)
  #norm(-temp,"2")
}


# the following can be used in both binary and count data.
rho_E1 <- function(init,ux,C_rho,G_rho,k_rho) {  
  #rho=p_iv$par
  rho=init
  u_sc=ux
  C=C_rho
  #C=1/3
  # get a better calculation
  G1=G_rho
  k=k_rho
  u_jc=u_sc[k*(1:G1)]
  u_ic=u_sc[k*(1:G1)+1]
  
  d=matrix(0, nrow= G1, ncol= 1)
  for (dd in 1:G1){
    d[dd]=loc_diff_abs[(k*dd),k*dd+1]
  }
  
  temp=(u_jc)*(u_ic)-C*rho*(1-d) 
  #temp=(u_jc)*(u_ic)-C*rho 
  
  norm(temp,"2")
}
rho_E2<- function(init,ux,C_rho,G_rho,k_rho) {  
  #rho=p_iv$par
  rho=init
  u_sc=ux
  C=C_rho
  # get a better calculation
  G1=G_rho
  k=k_rho
  u_jc=u_sc[k*(1:G1)]
  u_ic=u_sc[k*(1:G1)+2]
  
  d=matrix(0, nrow= G1, ncol= 1)
  for (dd in 1:G1){
    d[dd]=loc_diff_abs[(k*dd),(k*dd+2)]
  }
  
  temp=(u_jc)*(u_ic)-C*rho*(1-d) 
  #temp=(u_jc)*(u_ic)-C*rho 
  
  norm(temp,"2")
}
rho_E3<- function(init,ux,C_rho,G_rho,k_rho) {  
  #rho=p_iv$par
  rho=init
  u_sc=ux
  C=C_rho
  #C=1/3
  # get a better calculation
  G1=G_rho-1
  k=k_rho
  u_jc=u_sc[k*(1:G1)]
  u_ic=u_sc[k*(1:G1)+3]
  
  d=matrix(0, nrow= G1, ncol= 1)
  for (dd in 1:G1){
    d[dd]=loc_diff_abs[(k*dd),(k*dd+3)]
  }
  
  temp=(u_jc)*(u_ic)-C*rho*(1-d) 
  #temp=(u_jc)*(u_ic)-C*rho 
  
  norm(temp,"2")
}

P_E_IV_bi <- function(x_ini) {  #very very slow.
  #beta=beta_0
  beta=x_ini
  #W_in=w_inv(W,beta)
  # test , the test result is the exactly the same
  #W_in=diag(1/as.vector(pnorm(X%*%beta_0)*(1-pnorm(X%*%beta_0)))) 
  #temp=t(X)%*%W_in%*%(exp(X%*%(beta))*(Y-exp(X%*%(beta)))) # this is wrong, see sketch
  temp=t(diag(as.vector(dnorm(X%*%(beta))))%*%X)%*%W_in%*%(Y-pnorm(X%*%(beta)))
  # kk=t(X)
  # kkk=exp(X%*%(beta))
  #norm(matrix(3,2,1, byrow = FALSE),"2")
  norm(temp,"2")
}

w_inv_binary_E <- function(beta) {
  
  #beta=beta_0
  #w=diag(1/as.vector(exp(X%*%(beta))^2))
  #w=W
  W1= matrix(0,N,N)
  W2= matrix(0,N,N)
  #W0=W
  
  #Residual
  yc=pnorm(X%*%beta)
  u_s=Y-yc
  # mean(u_s)
  
  # adjust for hetro-case. Now everything with variance 1
  u_sc=u_s/sqrt( pnorm(X%*%beta)*(1-pnorm(X%*%beta)) )
  #x=x/exp(X%*%beta)
  
  # how fine is the regression
  # G0_rho=N/2-1
  # k0_rho=2
  
  G0_rho=N-2
  k0_rho=1
  
  # not good
  # G0_rho=G-1
  # k0_rho=nss
  
  e_rho <- optimize(rho_E1, c(0,1),tol = 0.0001,ux=u_sc,C_rho=1,G_rho=G0_rho,k_rho=k0_rho)
  rho_e=e_rho$minimum
  rho_e2=rho_e3=rho_e
  
  e_rho2 <- optimize(rho_E2, c(0,1),tol = 0.0001,ux=u_sc,C_rho=1,G_rho=G0_rho,k_rho=k0_rho)
  e_rho3<- optimize(rho_E3, c(0,1),tol = 0.0001,ux=u_sc,C_rho=1,G_rho=G0_rho,k_rho=k0_rho)
  rho_e2=e_rho2$minimum
  rho_e3=e_rho3$minimum
  

  #var(u_sc)
  for    (group in 1:G)
  { 
    #The following is the version that
    #group=1
    #beta=beta_0
    nG1= ((group-1)*nss+1)
    nG2= ((group-1)*nss+2)
    nG3= ((group-1)*nss+3)
    nG4= ((group-1)*nss+4)
    # le
    W1[nG1,nG1]=pnorm(X[nG1,]%*%beta)*(1-pnorm(X[nG1,]%*%beta))
    W1[nG2,nG2]=pnorm(X[nG2,]%*%beta)*(1-pnorm(X[nG2,]%*%beta))
    W1[nG3,nG3]=pnorm(X[nG3,]%*%beta)*(1-pnorm(X[nG3,]%*%beta))
    W1[nG4,nG4]=pnorm(X[nG4,]%*%beta)*(1-pnorm(X[nG4,]%*%beta))
    
    # covariance 
  
    
    W1[nG1,nG2]=W1[nG2,nG1]=sqrt(W1[nG1,nG1]*W1[nG2,nG2])*rho_e*(1-loc_diff_abs[nG1,nG2])
    W1[nG1,nG3]=W1[nG3,nG1]=sqrt(W1[nG1,nG1]*W1[nG3,nG3])*rho_e2*(1-loc_diff_abs[nG1,nG3])
    W1[nG1,nG4]=W1[nG4,nG1]=sqrt(W1[nG1,nG1]*W1[nG4,nG4])*rho_e3*(1-loc_diff_abs[nG1,nG4])
    # W1[nG1,nG2]=W1[nG2,nG1]=0
    W1[nG2,nG3]=W1[nG3,nG2]=sqrt(W1[nG2,nG2]*W1[nG3,nG3])*rho_e*(1-loc_diff_abs[nG2,nG3])
    W1[nG2,nG4]=W1[nG4,nG2]=sqrt(W1[nG2,nG2]*W1[nG4,nG4])*rho_e2*(1-loc_diff_abs[nG2,nG4])
    
    W1[nG3,nG4]=W1[nG4,nG3]=sqrt(W1[nG3,nG3]*W1[nG4,nG4])*rho_e*(1-loc_diff_abs[nG3,nG4])
    
    W2[nG1:nG4,nG1:nG4]= solve( W1[nG1:nG4,nG1:nG4])# why should they multiply it
    # W[nG1:nG2,nG1:nG2]=  xxx[nG1:nG2,nG1:nG2]*vargamma.true[nG1:nG2,nG1:nG2]
  }
  #KK=W2%*%W1
  W2
  # calculate the correlation of y
  
}

w_inv_binary_sim<- function(beta,m1) {
  
  #m1=10000
  #beta=beta_0
  
  W1= matrix(0,N,N)
  W2= matrix(0,N,N)
  NumIt1=m1
  y1=matrix(0, nrow= n, ncol=   NumIt1)
  Y1=matrix(0, nrow= n, ncol=   NumIt1)
  u1=matrix(0, nrow= n, ncol=   NumIt1)
  v1=matrix(0, nrow= n, ncol=   NumIt1)
  for (jj in 1: NumIt1){
    #j=1
    sim = 100+jj
    #set.seed(sim)
    vw2=wv_normal(class)
    v1[,jj] =vw2$V
    y1[,jj] = X%*%beta+v1[,jj]
    Y1[,jj] = (y1[,jj]>0)
    #u1[,j] =( Y1[,j]-pnorm(X%*%beta) )/sqrt( pnorm(X%*%beta)*(1-pnorm(X%*%beta)) )
    u1[,jj] =( Y1[,jj]-pnorm(X%*%beta) )
    mm=pnorm(X%*%beta)
    mmm=X%*%beta
    #mean(u1)
  }
  for    (group in 1:G)
  { 
    #The following is the version that
    #group=1
    #beta=beta_0
    nG1= ((group-1)*nss+1)
    nG2= ((group-1)*nss+2)
    nG3= ((group-1)*nss+3)
    nG4= ((group-1)*nss+4)
  
    W1[nG1,nG1]=pnorm(X[nG1,]%*%beta)*(1-pnorm(X[nG1,]%*%beta))
    W1[nG2,nG2]=pnorm(X[nG2,]%*%beta)*(1-pnorm(X[nG2,]%*%beta))
    W1[nG3,nG3]=pnorm(X[nG3,]%*%beta)*(1-pnorm(X[nG3,]%*%beta))
    W1[nG4,nG4]=pnorm(X[nG4,]%*%beta)*(1-pnorm(X[nG4,]%*%beta))
    
    # covariance
    W1[nG1,nG2]=W1[nG2,nG1]=cov(u1[nG1,],u1[nG2,])
    W1[nG1,nG3]=W1[nG3,nG1]=cov(u1[nG1,],u1[nG3,])
    W1[nG1,nG4]=W1[nG4,nG1]=cov(u1[nG1,],u1[nG4,])
    # W1[nG1,nG2]=W1[nG2,nG1]=0
    W1[nG2,nG3]=W1[nG3,nG2]=cov(u1[nG2,],u1[nG3,])
    W1[nG2,nG4]=W1[nG4,nG2]=cov(u1[nG2,],u1[nG4,])

    W1[nG3,nG4]=W1[nG4,nG3]=cov(u1[nG3,],u1[nG4,])
  
    
    W2[nG1:nG4,nG1:nG4]= solve( W1[nG1:nG4,nG1:nG4])# why should they multiply it
    
  }

W2
   
}
