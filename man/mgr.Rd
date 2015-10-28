\name{MGR}
\alias{mgr}
\title{
Multivariate Genome-wide Regression
}
\description{
Multivariate model to find breeding values through regression with optional resampling techniques.
}
\usage{
mgr(Y,gen,it=1000,bi=250,th=3,df=5,R2=0.5,pi=0,verb=TRUE)
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
  \item{pi}{
Value between 0 and 1. If greater than zero it activates variable selection, where markers have expected probability pi of having null effect. The model conjugates variable selection from a Beta-Binomial distribution.
}
  \item{df}{
Hyperprior degrees of freedom of variance components.
}
  \item{R2}{
Expected R2, used to calculate the prior shape as proposed by de los Campos et al. (2013).
}
  \item{verb}{
Logical. If verbose is TRUE, function displays MCMC progress bar.
}

}
\details{
The model for the whole-genome regression is as follows:

\eqn{Y = Mu + XG + E},

where \eqn{Y} is a matrix of response variables, \eqn{Mu} is the vector of intercepts, \eqn{X} is the genotypic matrix, \eqn{G} is a matrix that is the product of two terms (\eqn{g = BG}), \eqn{B} is the matrix of marker effects, \eqn{D} is a matrix of indicator variables that define whether or not the marker should be included into the model for each trait, \eqn{E} is the residual matrix.

Users can obtain two WGR methods out of this function: BRR (pi=0) and BayesC (pi>0). Variance components are sampled from a inverse Wishart distribution (Sorensen and Gianola 2002). The variable selection works through the unconditional prior algorithm proposed by Kuo and Mallick (1998). Gibbs sampler that updates regression coefficients is adapted from GSRU algorithm (Legarra and Misztal 2008).
}
\value{
The function mgr returns a list with the matrix of the expected marker effects for each trait (\eqn{B}), matrix of probability of each marker be included into the model of each trait (\eqn{D}), matrix of regression coefficients (\eqn{G}), random effect covariance matrix (\eqn{VA}), the vector of intercepts (\eqn{Mu}), residual covariance matrix (\eqn{VE}) and the fitted value (\eqn{Fit}).

}
\references{

Kuo, L., & Mallick, B. (1998). Variable selection for regression models. Sankhya: The Indian Journal of Statistics, Series B, 65-81.

de los Campos, G., Gianola, D., Rosa, G. J., Weigel, K. A., & Crossa, J. (2010). Semi-parametric genomic-enabled prediction of genetic values using reproducing kernel Hilbert spaces methods. Genetics Research, 92(04), 295-308.

}
\author{
Alencar Xavier
}
\examples{

# Data
data(tpod)
Y = cbind(rnorm(196,y,.1),rnorm(196,y,.2))

# BRR
test1 = mgr(Y,gen,it=250,bi=50)

# BayesC
test2 = mgr(Y,gen,it=250,bi=50,pi=0.8)

}