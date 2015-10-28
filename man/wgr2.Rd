\name{WGR2}
\alias{wgr2}
\title{
Whole-genome Regression 2
}
\description{
Multivariate model to find breeding values through whole-genome regression.
}
\usage{
wgr2(Y,gen,it=1000,bi=250,th=3,df=5,eigK=NULL,EigT=0.05,R2=0.5,pi=0,verb=FALSE)
}
\arguments{
  \item{Y}{
Numeric matrix of observations (\eqn{n,k}) describing the trait to be analyzed. \code{NA} is allowed.
}
  \item{gen}{
Numeric matrix containing the genotypic data. A matrix with \eqn{n}
rows of observations and (\eqn{m}) columns of molecular markers.
}
  \item{it}{
Integer. Number of iterations or samples to be generated.
}
  \item{bi}{
Integer. Burn-in, the number of iterations or samples to be discarted.
}
  \item{th}{
Integer. Thinning parameter, used to save memory by storing only one every 'th' samples.
}
  \item{df}{
Hyperprior degrees of freedom of variance components.
}
  \item{eigK}{
Output of function 'eigen'. Spectral decomposition of the kernel used to compute the polygenic term.
}
  \item{EigT}{
Null or numeric. If provided, the model uses just Eigenpairs with Eigenvalues above the specified theshold.
}
  \item{R2}{
Expected R2, used to calculate the prior shape as proposed by de los Campos et al. (2013).
}
  \item{pi}{
Value between 0 and 1. If greater than zero it activates variable selection, where markers have expected probability pi of having null effect. The model conjugates variable selection from a Beta-Binomial distribution.
}
  \item{verb}{
Logical. If verbose is TRUE, function displays MCMC progress bar.
}

}
\details{
This function solve multiple traits at once. It does not make use of kronecker products, but sample variance all components in multivariate fashion. The model for the whole-genome regression is as follows:

\eqn{Y = Mu + XG+ Z(UH) + E},

where \eqn{Y} is a matrix of response variables, \eqn{Mu} is the vector of intercepts, \eqn{X} is the genotypic matrix, \eqn{G} is a matrix that is the product of two terms (\eqn{g = BG}), \eqn{B} is the matrix of marker effects, \eqn{D} is a matrix of indicator variables that define whether or not the marker should be included into the model for each trait, \eqn{Z} is the design matrix, \eqn{U} is the matrix of Eigenvector of K, \eqn{H} is the matrix of regression coefficients of the polygenic term for different traits and \eqn{E} is the residual matrix.

Users can obtain two WGR methods out of this function: BRR (pi=0) and BayesC (pi>0). These methods can be fitted with or without a polygenic term solve via RKHS. Variance components are sampled from a inverse Wishart distribution (Sorensen and Gianola 2002). The variable selection works through the unconditional prior algorithm proposed by Kuo and Mallick (1998). Gibbs sampler that updates regression coefficients is adapted from GSRU algorithm (Legarra and Misztal 2008). Regression coefficients of the polygenic term are solved with the algorithm proposed by de los Campos et al. (2010).

}
\value{
The function wgr2 returns a list with the matrix of the expected marker effects for each trait (\eqn{B}), matrix of probability of each marker be included into the model of each trait (\eqn{D}), matrix of regression coefficients (\eqn{G}), matrix of polygenic coefficients (\eqn{U}), markers effect covariance matrix (\eqn{VA}), polygenic effect covariance matrix (\eqn{VK}), the vector of intercepts (\eqn{Mu}), residual covariance matrix (\eqn{VE}) and the fitted value (\eqn{Fit}).

}
\references{

Kuo, L., & Mallick, B. (1998). Variable selection for regression models. Sankhya: The Indian Journal of Statistics, Series B, 65-81.

de los Campos, G., Hickey, J. M., Pong-Wong, R., Daetwyler, H. D., and Calus, M. P. (2013). Whole-genome regression and prediction methods applied to plant and animal breeding. Genetics, 193(2), 327-345.

de los Campos, G., Gianola, D., Rosa, G. J., Weigel, K. A., & Crossa, J. (2010). Semi-parametric genomic-enabled prediction of genetic values using reproducing kernel Hilbert spaces methods. Genetics Research, 92(04), 295-308.

Legarra, A., & Misztal, I. (2008). Technical note: Computing strategies in genome-wide selection. Journal of dairy science, 91(1), 360-366.

Sorensen D., and Gianola D. (2002) Likelihood, Bayesian, and MCMC methods in quantitative genetics. Springer.

}
\author{
Alencar Xavier
}
\examples{

# G matrix
data(tpod)
G = tcrossprod(gen)
G = G/mean(diag(G))
eG = eigen(G,symmetric = TRUE)

# Phenotypes
Y1 = rnorm(196,y,.1)
Y2 = rnorm(196,y,.2)
Y3 = rnorm(196,y,.3)
Y = cbind(Y1,Y2,Y3)

# BRR
test1 = wgr2(Y,gen,it=200,bi=50)
diag(test1$VA/(test1$VA+test1$VE)) #h2
cov2cor(test1$VA) #Gen.cor.

# BayesC
test2 = wgr2(Y,gen,it=200,bi=50,pi=0.8)
diag(test2$VA/(test2$VA+test2$VE)) #h2
cov2cor(test2$VA) #Gen.cor.

# BayesC + RKHS
test3 = wgr2(Y,gen,eigK=eG,it=200,bi=50,pi=0.8)
diag((test3$VA+test3$VK)/(test3$VA+test3$VK+test3$VE)) #h2
cov2cor(test3$VA+test3$VK) #Gen.cor.

}