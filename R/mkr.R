# Multivariate Kernel Regression through RKHS

mkr = function(Y,K=NULL,eK=NULL,it=500,bu=200,th=3,
               df=-2,S=0,EigT=0.05,verb=TRUE){

 # it=100;bu=50;th=1;K=NULL;CE=TRUE;EigT=NULL;verb=TRUE
  if(is.null(K)&is.null(eK)) stop('Either K or eK have to be provided')
  if(is.null(eK)) eK = eigen(K)
  
# Inputs
U = eK$vectors
D = eK$values
if(!is.null(EigT)) {U=U[,1:sum(D>EigT)];D=D[D>EigT]}
k = ncol(Y)
n = nrow(Y)
p = ncol(U)
s = diag(S*df,k)
N = crossprod(!is.na(Y))
B = matrix(0,p,k)
E = apply(Y,2,function(x)x-mean(x,na.rm = T))
E[is.na(E)]=0
xx = rep(1,p)
mu = colMeans(Y,na.rm = T)
v = var(Y,na.rm = T)
VA = 0.1*diag(diag(var(Y,na.rm = T)))
VE = 0.1*diag(diag(var(Y,na.rm = T)))

# Store
mc = seq(bu,it,th)
lmc = length(mc)
mcMu = rep(0,k)
mcB = B
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
  
  mE = colSums(VA)
  mA = colSums(VE)
  
  # Update regression coefficients
  for(i in 1:k){
    
    ### TRY SOLVE AS IN BGLR
    lambda = max(mE[i]/mA[i],0)
    L = rep(lambda,p)
    up = KMUP(X=U[z[[i]],],b=B[,i],xx=xx,
              E=E[z[[i]],i],L=L,p=p,Ve=VE[i,i],pi=0)
    B[,i] = up[[1]]    
    err = up[[3]]
    mu[i] = rnorm(1,mu[i]+mean(err),VE[i,i]/N[i,i])
    E[z[[i]],i] = y[[i]]-U[z[[i]],]%*%B[,i]-mu[i]    
  }
  
  # Update variance components
  
    SSa = crossprod(B/D,B)
    VA = solve(rWishart(1,n+df,solve(SSa+s))[,,1])
    #VA = diag(diag(crossprod(B))/rchisq(k,diag(n-2)))
  
    VE = solve(rWishart(1,N+df,solve(crossprod(E)+s))[,,1])
    #VE = diag(diag(crossprod(E))/rchisq(k,diag(N+s)))

  
  if(j%in%mc){
    mcMu = mcMu + mu
    mcB = mcB + B
    mcVA = mcVA + VA
    mcVE = mcVE + VE
  }
  
  if(verb) setTxtProgressBar(pb, j/it)
}

if(verb) close(pb)

# Posterior means
mcMu = mcMu/lmc
mcB = mcB/lmc
mcVA = mcVA/lmc
mcVE = mcVE/lmc
A = U%*%mcB

final = list('BV'=A,'VA'=mcVA,'VE'=mcVE)
return(final)
}