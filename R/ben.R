
BEN = function(y,gen,it=750,bi=250,th=1,bag=0.80,alpha=0.5,wpe=50,MH=FALSE,verb=TRUE){
  
  X=gen;
  rm(gen);
  
  # Function to update beta
  upB = function(e,mu,X,b,l,a,xx,p,E2,X2,bag,pi,wpe,ve,O,mh=TRUE){
    xx = xx*bag
    pi = (0.5+a*99)/100
    if(mh){
      a_new = rbeta(1,200*pi,200*(1-pi))
      b1 = b2 = rep(NA,p)
      e1 = e2 = e
      L1_1 = (1-a)*l/(2*xx)
      L2_1 = (a)*(xx/(xx*l))
      L1_2 = (1-a_new)*l/(2*xx)
      L2_2 = (1-a_new)*(xx/(xx*l))
      # Update regression coefficients
      for(j in O){
        # Old alpha
        xy1 = (crossprod(e1,X[,j])+b[j]*xx[j])/xx[j]
        s1 = sign(xy1)
        beta = abs(xy1)-L1_1[j]
        ben = s1*beta*L2_1[j]
        b1[j] = ifelse(beta>0,rnorm(1,ben,ve/(xx[j])),0)
        e1 = e1 - X[,j]*(b1[j]-b[j])
        # New alpha
        xy2 = (crossprod(e2,X[,j])+b[j]*xx[j])/xx[j]
        s2 = sign(xy2)
        beta = abs(xy2)-L1_2[j]
        ben = s2*beta*L2_2[j]
        b2[j] = ifelse(beta>0,rnorm(1,ben,ve/(xx[j])),0)
        e2 = e2 - X[,j]*(b2[j]-b[j])
      }
      # Loss function
      SSPE_1 = sum(as.vector(tcrossprod(b1,X2)-E2)^2)
      SSPE_2 = sum(as.vector(tcrossprod(b2,X2)-E2)^2)
      LOSS1 = wpe*SSPE_1+crossprod(e1)+l*(0.5*crossprod(b1)*(1-a)+sum(abs(b1))*a)
      LOSS2 = wpe*SSPE_2+crossprod(e2)+l*(0.5*crossprod(b2)*(1-a_new)+sum(abs(b2))*a_new)
      LR = LOSS2/LOSS1
      if(is.na(LR)|is.nan(LR)) LR=0
      if(LR>1){
        P=list('b'=b2,'a'=a_new,'e'=e2,'oob'=SSPE_2)
      }else{
        if(MH){
          # Metropolis-Hastings
          if(LR>runif(1)){
            P=list('b'=b2,'a'=a_new,'e'=e2,'oob'=SSPE_2)
          }else{
            P=list('b'=b1,'a'=a,'e'=e1,'oob'=SSPE_1)
          }
        }else{
          # Acceptance-Rejection
          P=list('b'=b1,'a'=a,'e'=e1,'oob'=SSPE_1)
        }
      }
    }else{
      b1 = rep(NA,p)
      e1 = e
      L1_1 = (1-a)*l/(2*xx)
      L2_1 = (a)*(xx/(xx*l))
      # Update regression coefficients
      for(j in O){
        # Old alpha
        xy1 = (crossprod(e1,X[,j])+b[j]*xx[j])/xx[j]
        s1 = sign(xy1)
        beta = abs(xy1)-L1_1[j]
        ben = s1*beta*L2_1[j]
        b1[j] = ifelse(beta>0,rnorm(1,ben,ve/(xx[j])),0)
        e1 = e1 - X[,j]*(b1[j]-b[j])}
      # Loss function
      SSPE = sum(as.vector(tcrossprod(b1,X2)-E2)^2)
      P=list('b'=b1,'a'=a,'e'=e1,'oob'=SSPE)
    }
    return(P)
  }
  
  
  # Missing
  if(anyNA(y)){
    mis = which(is.na(y))
    Y = y
    XX = X[mis,]
    y = y[-mis]
    X = X[-mis,]
    MISS = TRUE
  }else{
    MISS = FALSE
  }
  # Data
  xx = apply(X,2,function(x)crossprod(x))
  b0 = crossprod(X,y)[,1]/xx
  O = order(b0^2,decreasing = TRUE)
  n = nrow(X)
  p = ncol(X)
  bn = round(n*bag)
  MCMC = seq(bi,it,th)
  MC = length(MCMC)
  # Parameters
  mu = mean(y)
  e = y-mu
  b = rep(0,p)
  a = alpha
  l = 1
  ve = 0.1
  # Store posterior
  B = rep(0,p)
  MU = A = L = SSPE = 0
  if(verb) pb = txtProgressBar(style = 3)
  # Loop
  for(i in 1:it){
    
    # Bagging
    s = sort(sample(1:n,n-bn,replace=FALSE))
    # UPDATE
    UP = upB(e[-s],mu,X[-s,],b,l,a,xx*bag,p,e[s],
             X[s,],bag,pi=a,wpe,ve,O=O,mh=i%%10==0)
    b = UP[[1]]
    a = UP[[2]]
    e = UP[[3]]
    
    mu = mu + mean(e)
    df_prior = 2+rpois(1,3)
    Se = runif(1,0,1)
    Sb = runif(1,0,1)*df_prior
    ve = (crossprod(e)+Se)/(n+2)
    vb = (crossprod(b)+Sb)/(p+df_prior)
    l = ve/vb
    e = as.vector(y-(mu+tcrossprod(b,X)))
    
    # STORE
    if(i%in%MCMC){
      B = B+b
      MU = MU+mu
      A = A+a
      L = L+l
      if(bag<1) SSPE = SSPE+UP$oob/(n-bn)
    }
    if(verb) setTxtProgressBar(pb, i/it)
  }
  if(verb) close(pb)
  # Posterior
  Bhat = B/MC
  MUhat = MU/MC
  Ahat = A/MC
  Lhat = L/MC
  MSPEout = mean(SSPE)/MC
  # Prediction
  if(MISS){
    Yhat = Y
    Yhat[-mis] = as.vector(mu+tcrossprod(Bhat,X))
    Yhat[mis] = as.vector(mu+tcrossprod(Bhat,XX))
  }else{
    Yhat = as.vector(mu+tcrossprod(Bhat,X))
  }
  # OUTPUT
  LIST = list('hat'=Yhat,'coef'=Bhat,'b0'=MUhat,
              'alp'=Ahat,'lmb'=Lhat,'MSPEoob'=MSPEout)
  
}

GELA = function(y,x,It=500,CL=3,PC=15,PC2=5,Par=15,Bg=0.75,rdg=1e-6,weighted=TRUE){
  if(!any(is.na(y))) stop("There are no missing values in Y to be predicted")
  if(any(is.na(x))) stop("There are missing values among predictors")
  mi = which(is.na(y))
  nmi = length(mi)
  vx = apply(x,2,var)
  if(any(vx<0.1)) x=x[,-which(vx<0.1)]
  K = tcrossprod(x)
  K = K/mean(diag(K))
  eig = eigen(K,TRUE)
  dr = which(cumsum(eig$values)/ncol(K)>0.98)[1]
  E = eig$vectors[-mi,1:dr]
  ee = eig$vectors[mi,1:dr]
  eig$values = eig$values[1:dr]
  Y = y[-mi]
  X = x[-mi,]
  xx = x[mi,]
  R2x = abs(cor(Y,X)^2)
  R2x = as.vector(R2x/sum(R2x))
  R2e = abs(cor(Y,E)^2)
  R2e = as.vector(R2e/sum(R2e))
  R2v = as.vector(eig$values/sum(eig$values))
  Bag = round(Bg*length(Y))
  yhat = matrix(NA,nmi,It)
  Xb = matrix(NA,Bag,Par)
  Eb = matrix(NA,Bag,PC)
  Pb = matrix(NA,Bag+nmi,PC2)
  n = length(Y)
  nW = (1+PC+CL*Par)
  Wp = matrix(NA,Bag,nW)
  Wv = matrix(NA,nmi,nW)
  WR2 = rep(NA,It)
  pb = txtProgressBar(style = 3)
  for(i in 1:It){
    o = sample(1:n,Bag)
    t = sample(1:ncol(X),Par,prob=R2x)
    g = sample(1:dr,PC,prob=R2e)
    p = sample(1:dr,PC2,prob=R2v)
    Xb[1:Bag,1:Par] = X[o,t]
    Eb[1:Bag,1:PC] = E[o,g]
    Pb[1:Bag,1:PC2] = E[o,p]
    Pb[(Bag+1):(Bag+nmi),1:PC2] = ee[,p]
    K = factor(suppressWarnings(kmeans(Pb,CL))$cluster)
    K1 = K[1:Bag]
    K2 = K[(Bag+1):(Bag+nmi)]
    Wp[1:Bag,1:nW] = model.matrix(~Eb+K1:Xb)
    Wv[1:nmi,1:nW] = model.matrix(~ee[,g]+K2:xx[,t])
    WW = crossprod(Wp)
    diag(WW) = diag(WW)+c(0,rep(rdg,nW-1))
    Wy = crossprod(Wp,Y[o])
    fit = suppressWarnings(solve(WW,Wy))
    WR2[i] = mean((Wp%*%fit-Y[o])^2)
    yhat[,i] = Wv %*% fit
    setTxtProgressBar(pb, i/It)}
  close(pb)
  WR2 = 1-WR2/max(WR2)
  if(weighted){
    Hat = apply(yhat,1,weighted.mean,w=WR2)
  }else{
    Hat = apply(yhat,1,mean)
  }
  yFit = y
  yFit[mi] = Hat
  return(yFit)
}
