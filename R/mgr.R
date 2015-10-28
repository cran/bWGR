mgr = function(Y,gen,it=1000,bi=250,th=3,
               df=5,R2=0.5,pi=0,verb=TRUE){
  
  X = gen
  xx = colSums(X^2)
  MSx = sum(apply(X,2,var,na.rm=T))
  S_prior = R2*apply(Y,2,var,na.rm=T)*(df+2)/MSx/(1-pi)
  k = ncol(Y)
  n = nrow(Y)
  p = ncol(X)
  s = diag(S_prior,k)
  N = crossprod(!is.na(Y))
  B = G = D = matrix(0,p,k)
  E = apply(Y,2,function(x)x-mean(x,na.rm = T))
  E[is.na(E)]=0
  md = rep(pi,k)
  mu = colMeans(Y,na.rm = T)
  v = var(Y,na.rm = T)
  VA = 0.4*diag(diag(var(Y,na.rm = T)))
  VE = 0.1*diag(diag(var(Y,na.rm = T)))
  
  # Store
  mc = seq(bi,it,th)
  lmc = length(mc)
  mcMu = rep(0,k)
  mcB = mcG = mcD = B
  mcVA = mcVE = matrix(0,k,k)
  
  # Indicators
  y = z = list()
  for(i in 1:k){
    z[[i]] = which(!is.na(Y[,i]))
    y[[i]] = Y[z[[i]],i]
  }
  
  # MCMC
  if(verb) pb = txtProgressBar(style=3)
  
  for(j in 1:it){
    
    mEA = colSums(solve(VA)%*%VE)
    #mE = colSums(VA);mA = colSums(VE)
    if(pi>0){PI = rbeta(1,10*pi+md+1,10*(1-pi)-md+1)
    }else{PI=0}
    
    # Update regression coefficients
    for(i in 1:k){
      
      lambda = max(mEA[i],0)
      #lambda = max(mE[i]/mA[i],0)      
      L = rep(lambda,p)
      up = KMUP(X[z[[i]],],B[,i],xx,E[z[[i]],i],L,p,VE[i,i],PI)
      b = up[[1]]
      d = up[[2]]
      e = up[[3]]
      if(pi>0) d[is.nan(d)] = 1
      B[,i] = b
      D[,i] = d
      G[,i] = up[[1]]*d     
      mu[i] = rnorm(1,mu[i]+mean(e),VE[i,i]/N[i,i])
      E[z[[i]],i] = y[[i]]-X[z[[i]],]%*%B[,i]-mu[i]
      md[i] = mean(D[,i])
    }
    
    # Update variance components
    VA = solve(rWishart(1,n+df,solve(crossprod(B)+s))[,,1])
    VE = solve(rWishart(1,N+2,solve(crossprod(E)+s))[,,1])
    
    if(j%in%mc){
      mcMu = mcMu + mu
      mcB = mcB + B
      mcG = mcG + G
      mcD = mcD + D
      mcVA = mcVA + VA
      mcVE = mcVE + VE
    }
    
    if(verb) setTxtProgressBar(pb, j/it)
  }
  
  if(verb) close(pb)
  
  # Posterior means
  mcMu = mcMu/lmc
  mcB = mcB/lmc
  mcG = mcG/lmc
  mcD = mcD/lmc
  mcVA = mcVA/lmc
  mcVE = mcVE/lmc
  A = t(mcMu+t(X%*%mcB))
  
  final = list('Mu'=mcMu,'B'=mcB,'G'=mcG,'D'=mcD,
               'VA'=mcVA,'VE'=mcVE,'Fit'=A)
  return(final)
}
