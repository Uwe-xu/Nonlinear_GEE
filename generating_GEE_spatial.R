# simulator functions

nss = N/G
locs = get.loc(N, G, max_range = 10, sd = 0.01)
loc_diff_abs = get.loc_dist(locs) # distance matrix.

Sig_g<- function(klass,r) {  
  Sig=matrix(0, nrow= N, ncol= N)
  
  #r=0.9
  if (klass==0) {
    #block=matrix(r, nrow= nss, ncol= nss)
    block=matrix(r, nrow= nss, ncol= nss)
    for (i in 1:nss)
    {
      block[i,i]=1
    }
    for (group in 1:G)
    {
      nG1= ((group-1)*nss+1)
      nG2= ((group-1)*nss+nss)
      Sig[nG1:nG2,nG1:nG2]=  block 
    } 
    
  } else 
  {
    #r=10
    Sig = get.loc.cov(loc_diff_abs, rho = r) # original is 20
    
  }
  Sig
}

wv_normal <- function(klass) { 
  W = matrix(0,N,N)
  if (klass==0) {
    W=diag(n) # we need vvv~N(-1,2)
  } else if(klass==1)
  {
    for (group in 1:G)
    {
      nG1= ((group-1)*nss+1)
      nG2= ((group-1)*nss+nss)
      W[nG1:nG2,nG1:nG2]=  Sig[nG1:nG2,nG1:nG2] # why should they multiply it
    } 
  }
  else{W=pmax(Sig,0) # some problem here, Sig is generated in this way it could be smaller than zero
  }
  #option 1 very slow£¬
  vvv=matrix(mvrnorm(1,rep(0,n),W),n,1, byrow = FALSE)
  
  #option 2 not stable, seem's problem
  
  #normalized w
  # wei=1/colSums(W)#weight
  # #typeof(wei)
  # W=diag(sqrt(wei))%*%W
  # #colSums(W)
  # vv=matrix(rnorm(1*n,mean=0,sd=1),n,1, byrow = FALSE) # this one is far more efficient
  # vvv=W%*%vv

  v=vvv
  list(V = v,W=W)
}
# simulating v with log_normal from N(-1,2)

wv_log_normal <- function(klass) {  
  W = matrix(0,N,N)
  if (klass==0) {
    W=diag(n) # we need vvv~N(-1,2)
  } else 
  {
    
    for (group in 1:G)
    {
      nG1= ((group-1)*nss+1)
      nG2= ((group-1)*nss+nss)
      W[nG1:nG2,nG1:nG2]=  Sig[nG1:nG2,nG1:nG2] # why should they multiply it
    } 
  }
  # vv=matrix(rnorm(1*n,mean=0,sd=1),n,1, byrow = FALSE)
  # vvv=W%*%vv # this is wrong
  vvv=matrix(mvrnorm(1,rep(0,n),W),n,1, byrow = FALSE)
  vvv=vvv-1/2
  
  # vv=matrix(rnorm(1*n,mean=0,sd=sqrt(0.2)),n,1, byrow = FALSE) # this one is far more efficient
  # vvv=W%*%vv
  # vvv=sqrt(0.2)*matrix(mvrnorm(1,rep(0,n),W),n,1, byrow = FALSE)
  # vv=vvv-0.1
  
  v=exp(vvv)
  list(V = v,W=W)
  
}

# simulating v with uniform from 1+U(-1,1)
wv_uniform <- function(klass) {  
  W = matrix(0,N,N)
  cw=matrix(0,G,1)
  if (klass==0) {
    W=diag(n) # we need vvv~N(-1,2)
    vvv=matrix(runif(n, min = -1, max = 1),n,1, byrow = FALSE)
  } else 
  {
    temp=matrix(runif(n, min = -1, max = 1),n,1, byrow = FALSE)
    for (group in 1:G)
    {
      nG1= ((group-1)*nss+1)
      nG2= ((group-1)*nss+nss)
      W[nG1:nG2,nG1:nG2]=  Sig[nG1:nG2,nG1:nG2] # why should they multiply it
      # the following is only for the 2 member cases.
      cw[group]=W[nG2,nG1]
    } 
    U=rbinom(G,1,cw)
    #temp1=temp[2*(1:G)-1]
    # the problem of this structure is we got series with the same element at many places.
    temp[2*(1:G)]=temp[2*(1:G)-1]*U+temp[2*(1:G)]*(1-U)
    #cov(temp[2*(1:G)-1],temp[2*(1:G)])
    vvv=temp
  }
  v=vvv+1
  list(V = v,W=W)
  
}