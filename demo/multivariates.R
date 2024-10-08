
Z = SimZ()
GC = SimGC()
S = SimY(Z=Z, GC=GC)
Y = S$Y
tbv = S$tbv

system.time(fit_mega <- MEGA(Y,Z)$gebv)[3]
system.time(fit_gsem <- GSEM(Y,Z)$hat)[3]
system.time(fit_fuvbeta <- Z %*% FUVBETA(Y,Z))[3]
system.time(fit_mrrf <- MRR3F(Y,Z)$hat)[3]
system.time(fit_xfa <- MRR3(Y,Z,XFA=T)$hat)[3]
system.time(fit_sem <- SEM(Y,Z)$hat)[3]
system.time(fit_xsemf <- XSEMF(Y,Z)$hat)[3]
system.time(fit_zsemf <- ZSEMF(Y,Z)$hat)[3]
system.time(fit_ysemf <- YSEMF(Y,Z)$hat)[3]

sort(sapply(grep('fit',ls(),value=T),function(x) mean(diag(cor(get(x),tbv)))))
