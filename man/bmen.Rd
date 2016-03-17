\name{BMEN}
\alias{BagMEN}
\alias{GSEN}
\alias{GELA}
\alias{wGELA}
\title{
Bagging Metropolis Elastic Net (BMEN)
}
\description{
Univariate model regularized via elastic net to find breeding values through regression with resampling techniques.
}
\usage{
BagMEN(y,X,it=500,bi=100,bag=0.5,alpha=0.5,wpe=1)
}
\arguments{
  \item{y}{
Numeric vector of observations (\eqn{n}) describing the trait to be analyzed. \code{NA} is allowed.
}
  \item{X}{
Numeric matrix containing the genotypic data. A matrix with \eqn{n}
rows of observations and (\eqn{m}) columns of molecular markers.
}
  \item{it}{
Integer. Number of iterations or samples to be generated.
}
  \item{bi}{
Integer. Burn-in, the number of iterations or samples to be discarted.
}
  \item{bag}{
Proportion of data used to bag MCMC.
}
  \item{alpha}{
Value between 0 and 1. An elastic net with alpha = 1 is a ridge regression, and alpha = 0 is a LASSO. Proposed values of alphas for the Manhattan algorithm are sampled from a beta distribution with expectation provided by this argument.
}
  \item{wpe}{
Weight assigned for out-of-bag sum of squared prediction error. Important for the update of alpha via Metropolis algorithm.
}

}
\details{
The model for the whole-genome regression is as follows:

\eqn{y = mu + Xb + e},

where \eqn{y} is the response variable, \eqn{mu} is the intercept, \eqn{X} is the genotypic matrix, \eqn{b} is the effect of an allele substitution and \eqn{e} is the residual term.

This function minimizes the loss:

\eqn{l = (e'e)+(w*f'f)+L(A*(b'b)+(1-A)*sum(abs(b))}

Thus minimizing both goodness-of-fit error (e) and out-of-bag prediction error (f), where f can be present (w=1) or absent (w=0). Regression coefficients are updated via coordinate descent. L is computed using variances (VarE/VarB) and A is found through a Metropolis algorithm.

}
\value{
The function BagMEN returns a list with regression coefficient (\eqn{coef}), regularization parameters alpha (\eqn{alp}) and lambda (\eqn{lmb}), the intercept (\eqn{b0}), the out-of-bag mean square prediction error (\eqn{MSPEoob}) and the fitted value (\eqn{hat}).

}

\author{
Alencar Xavier
}
\examples{

data(tpod)
gen = gen[,seq(1,376,5)]

# Fit BagMEN 
fit = BagMEN(y,gen)
cor(y,fit$hat)

}
