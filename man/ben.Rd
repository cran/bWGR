\name{BEN}
\alias{BEN}
\alias{GSEN}
\alias{GELA}
\title{
Bagging Elastic Net (BEN)
}
\description{
Univariate model regularized via elastic net to find breeding values through regression with resampling techniques.
}
\usage{
BEN(y,gen,it=750,bi=250,th=1,bag=0.80,
    alpha=0.5,wpe=50,MH=FALSE,verb=TRUE)
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
  \item{verb}{
Logical. If verbose is TRUE, function displays MCMC progress bar.
}
  \item{alpha}{
Numeric between 0 and 1. Starting value of alpha parameter for the bagging elastic net model, where 0 is ridge (L2) and 1 is lasso (L1).
}
  \item{wpe}{
Weight of prediction error sum of squared for the jump function that defines alpha in the bagging elastic net model.
}
  \item{MH}{
Logical. If TRUE, the search for alpha is performed via Metropolis-Hastings. If FALSE, acceptance-rejection.
}
}
\details{
The model for the whole-genome regression is as follows:

\eqn{y = mu + Xb + e},

where \eqn{y} is the response variable, \eqn{mu} is the intercept, \eqn{X} is the genotypic matrix, \eqn{b} is the effect of an allele substitution and \eqn{e} is the residual term.

This function minimizes the loss:

\eqn{l = (e'e)+(w*f'f)+L(A*(b'b)+(1-A)*sum(abs(b))}

Thus minimizing both goodness-of-fit error (e) and out-of-bag prediction error (f), where f can be present (w>=1) or absent (w=0). Regression coefficients are updated via coordinate descent. L is computed using variances (VarE/VarB) and A is found through Acceptance-Rejection or Metropolis algorithm.

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

# Fit BEN
fit = BEN(y,gen,200,50)
cor(y,fit$hat)
plot(abs(fit$coef))

}
