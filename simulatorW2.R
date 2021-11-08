#install.packages("CovTools")
#install.packages("rgl")



# 1. generate locations

get.loc<-function(N, G, max_range, sd)
{
  n_g = floor(N/G)
  loc_center = seq(0, max_range, length.out = G)
  locs = sapply(loc_center, function(x){
    rnorm(n_g, mean = x, sd = sd)
  })
  return(as.vector(locs))
}

get.loc_dist<-function(locs)
{
  loc_diff_abs = abs(outer(locs, locs, "-"))
  return(loc_diff_abs)
}

get.loc.cov<-function(loc_diff_abs, rho)
{ 
  #rho=0.9
  #Sig = 1-rho
  #Sig = exp(-rho*loc_diff_abs) # this is the equation defines the sigma, before blocking.
  #Sig = 1-rho*loc_diff_abs # this is somehow similar to the one above
  
  # we let rho multiplied with it from outside.
  Sig = rho*(1-loc_diff_abs)
  Sig=Sig+diag(rep(1-rho,n))
  
  
  return(Sig)
}

eigen.Sig<-function(Sig)
{
  eig_S = eigen(Sig)
  
}

get.gamma<-function(eig_S)
{
  N = length(eig_S$values)
  eps = rnorm(N)
  eps_cov = as.vector(eig_S$vectors%*%(sqrt(eig_S$values)*eps))
  probs = pnorm(eps_cov)
  eps_gam = qgamma(probs, shape = 1, rate = 1)
  return(eps_gam)
}


simu.poi<-function(X, beta, eps_gam)
{
  y0 = eps_gam*exp(X%*%beta)
  y = sapply(y0, rpois, n = 1)
  return(y)
}


simu.expCov.trueW<-function(eig_S, X, beta0, sd_X = 0.5)
{
  Nrep = 500
  N = length(eig_S$values)
  
  ymat = matrix(0, nrow = N, ncol = Nrep)
  set.seed(1234)
  for (r in 1:Nrep)
  {
    cat(r, "\r")
    eps_gam = get.gamma(eig_S)
    ymat[,r] = simu.poi(X, beta0, eps_gam)
  }
  ymat[is.na(ymat)] = max(ymat, na.rm = T) +5
  y_dmean = ymat - rowMeans(ymat)
  W = y_dmean%*%t(y_dmean)/Nrep
  return(W)
}





simu.expCov.trueW.w2<-function(eig_S, X, beta0, sd_X = 0.5)
{
  Nrep = 2000
  N = length(eig_S$values)
  
  # ymat = matrix(0, nrow = N, ncol = Nrep)
  
  ns = 500
  eps_gam_m = matrix(1, N, ns)
  for (sample in 1:ns)
  {
    
    eps_gam_m[,sample] = get.gamma(eig_S)
    
  } 
  
  vargamma.true = Sig
  vargamma = cov(t(eps_gam_m))
  
  cc       = abs(as.vector(vargamma))
  thresh   = 0.5
  
   vargamma[abs(vargamma)< thresh] <- 0
  
  #   hist(cc)
  uu        = exp(X%*%beta0) 
  thresh2   = mean(uu)
        uu[abs(uu)>thresh/2] <- thresh2/2
       uu[abs(uu)> 4*min(uu)] <- 4*min(uu)
   xxx= uu%*% t(uu)
   
   nss = N/G
  W = matrix(0,N,N)
   for    (group in 1:G)
   {
     nG1= ((group-1)*nss+1)
     nG2= ((group-1)*nss+nss)
     W[nG1:nG2,nG1:nG2]=  xxx[nG1:nG2,nG1:nG2]*vargamma[nG1:nG2,nG1:nG2]
     # W[nG1:nG2,nG1:nG2]=  xxx[nG1:nG2,nG1:nG2]*vargamma.true[nG1:nG2,nG1:nG2]
     
   }   
   
  #varr = diag(W)
  #invv = solve(sqrt(diag(varr)))
  #corre = invv%*%W%*%invv
  #corrplot(corre[1:100,1:100])
  #for (r in 1:Nrep)
  # {
  #  cat(r, "\r")
  # eps_gam = get.gamma(eig_S)
  # ymat[,r] = simu.poi(X, beta0, eps_gam)
  # }
  # ymat[is.na(ymat)] = max(ymat, na.rm = T) +5
  # y_dmean = ymat - rowMeans(ymat)
  # W = y_dmean%*%t(y_dmean)/Nrep
  return(W)
}




#####################################

get.M<-function(loc_diff_abs, rho_M = 10)
{
  M0 = exp(-rho_M*loc_diff_abs)
  M1 = M0/rowSums(M0)
  return(M1)
}

simu.e<-function(M1, rho_e = 0.5, sig_eps = 1)
{
  N = nrow(M1)
  eps = rnorm(N, sd = sig_eps)
  I = diag(N)
  e = M1%*%eps
  return(e)
}

get.v<-function(e, M1, rho_e, sig_eps)
{
  N = nrow(M1)
  I = diag(N)
  cov_e = sig_eps^2*M1%*%t(M1)
  var_e = diag(cov_e)
  v = e^2/var_e
  return(v)
}

simu.y<-function(X, beta, v)
{
  y0 = v*exp(X%*%beta)
  y = sapply(y0, rpois, n = 1)
  return(y)
}



simu.trueW<-function(M1, rho_e, X, beta0, sig_eps)
{
  Nrep = 500
  N = nrow(M1)
  
  ymat = matrix(0, nrow = N, ncol = Nrep)
  set.seed(1234)
  for (r in 1:Nrep)
  {
    cat(r, "\r")
    e = simu.e(M1, rho_e, sig_eps)
    v = get.v(e, M1, rho_e, sig_eps)
    ymat[,r] = simu.y(X, beta0, v)
  }
  y_dmean = ymat - rowMeans(ymat)
  W = y_dmean%*%t(y_dmean)/Nrep
  return(W)
}