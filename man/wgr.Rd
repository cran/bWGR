\name{WGR}
\alias{wgr}
\alias{KMUP}
\alias{.Random.seed}
\title{
Whole-genome Regression
}
\description{
Univariate model to find breeding values through regression with optional resampling techniques.
}
\usage{
wgr(y,gen,it=1500,bi=500,th=1,bag=1,rp=TRUE,iv=FALSE,
    pi=0,df=5,R2=0.5,eigK=NULL,rankK=0.25,verb=FALSE)
}
\arguments{
  \item{y}{
Numeric vector of observations (\eqn{n}) describing the trait to be analyzed. \code{NA} is allowed.
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
  \item{bag}{
Proportion of data used to bag MCMC.
}
  \item{rp}{
Logical. Use replacement (bootstrap) samples when bag is different than one.
}
  \item{iv}{
Logical. Assign markers independent variance. If true, turns the default model BLUP into BayesA. For this model, the shape parameter is conjugated by a gamma with hyperpriors calculated based on the R2 rule.
}
  \item{pi}{
Value between 0 and 1. If greater than zero it activates variable selection, where markers have expected probability pi of having null effect. The model conjugates variable selection from a Beta-Binomial distribution.
}
  \item{df}{
Hyperprior degrees of freedom of variance components.
}
  \item{R2}{
Expected R2, used to calculate the prior shape as proposed by de los Campos et al. (2013).
}
  \item{eigK}{
Output of function 'eigen'. Spectral decomposition of the kernel used to compute the polygenic term.
}
  \item{rankK}{
Numeric between 0 and 1. Indicates the proportion of Eigenpairs used to fit the polygenic term.
}
  \item{verb}{
Logical. If verbose is TRUE, function displays MCMC progress bar.
}

}
\details{
The model for the whole-genome regression is as follows:

\eqn{y = mu + Xg + u + e},

where \eqn{y} is the response variable, \eqn{mu} is the intercept, \eqn{X} is the genotypic matrix, \eqn{g} is the product of two terms (\eqn{g = bg}), \eqn{b} is the effect of an allele substitution, \eqn{d} is an indicator variable that define whether or not the marker should be included into the model, \eqn{u} is the polygenic term and \eqn{e} is the residual term.

Users can obtain four WGR methods out of this function: BRR (pi=0,iv=F), BayesA (pi=0,iv=T), BayesB (pi=0.8,iv=T) and BayesC (pi=0.8,iv=F). The full theoretical basis of each model is described by de los Campos et al. (2013).

Gibbs sampler that updates regression coefficients is adapted from GSRU algorithm (Legarra and Misztal 2008). The variable selection works through the unconditional prior algorithm proposed by Kuo and Mallick (1998). The polygenic term is solved by Bayesian algorithm of reproducing kernel Hilbert Spaces proposed by de los Campos et al. (2010).
}
\value{
The function wgr returns a list with expected value from the marker effect (\eqn{b}), probability of marker being in the model (\eqn{d}), regression coefficient (\eqn{g}), variance of each marker (\eqn{Vb}), the intercept (\eqn{mu}), the polygene (\eqn{u}) and polygenic variance (\eqn{Vk}), residual variance (\eqn{Ve}) and the fitted value (\eqn{hat}).

variance components distribution a posteriori (Posterior.VC) and mode estimated (VC.estimate), a list with the posterior distribution of regression coefficients (Posterior.Coef) and the posterior mean (Coef.estimate), and the fitted values using the mean (Fit.mean) of posterior coefficients. 
}
\references{

de los Campos, G., Hickey, J. M., Pong-Wong, R., Daetwyler, H. D., and Calus, M. P. (2013). Whole-genome regression and prediction methods applied to plant and animal breeding. Genetics, 193(2), 327-345.

de los Campos, G., Gianola, D., Rosa, G. J., Weigel, K. A., & Crossa, J. (2010). Semi-parametric genomic-enabled prediction of genetic values using reproducing kernel Hilbert spaces methods. Genetics Research, 92(04), 295-308.

Kuo, L., & Mallick, B. (1998). Variable selection for regression models. Sankhya: The Indian Journal of Statistics, Series B, 65-81.

Legarra, A., & Misztal, I. (2008). Technical note: Computing strategies in genome-wide selection. Journal of dairy science, 91(1), 360-366.

}
\author{
Alencar Xavier
}
\examples{

data(tpod)
gen = gen[,seq(1,376,5)]

# BLUP
BRR = wgr(y,gen,iv=FALSE,pi=0,bag=0.5,rp=TRUE,it=400,bi=50)
cor(y,BRR$hat)

# BayesA
BA = wgr(y,gen,iv=TRUE,pi=0,bag=0.5,rp=TRUE,it=400,bi=50)
cor(y,BA$hat)

# BayesB
BB = wgr(y,gen,iv=TRUE,pi=.5,bag=0.5,rp=TRUE,it=400,bi=50)
cor(y,BB$hat)

# BayesC
BC = wgr(y,gen,iv=FALSE,pi=.5,bag=0.5,rp=TRUE,it=400,bi=50)
cor(y,BC$hat)
}