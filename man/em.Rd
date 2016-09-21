\name{WGR (EM)}
\alias{emBA}
\alias{emBB}
\alias{emBC}
\alias{emRR}
\title{
Whole-genome Regression
}
\description{
Univariate model to find breeding values through regression fitted via expectation-maximization implemented in C++.
}
\usage{
emRR(y, gen, df = 5, R2 = 0.5, it = 75L)
emBA(y, gen, df = 5, R2 = 0.5, it = 75L)
emBB(y, gen, df = 5, R2 = 0.5, it = 75L, Pi = 0.95)
emBC(y, gen, df = 5, R2 = 0.5, it = 75L, Pi = 0.95)
}
\arguments{
  \item{y}{
Numeric vector of observations (\eqn{n}) describing the trait to be analyzed. \code{NA} is not allowed.
}
  \item{gen}{
Numeric matrix containing the genotypic data. A matrix with \eqn{n}
rows of observations and (\eqn{m}) columns of molecular markers.
}
  \item{df}{
Hyperprior degrees of freedom of variance components.
}
  \item{R2}{
Expected R2, used to calculate the prior shape as proposed by de los Campos et al. (2013).
}
  \item{it}{
Integer. Number of iterations to fit the model. Most models converge within 75 interations. Large datasets (p>10k and n>20k) may need more.
}
  \item{Pi}{
Value between 0 and 1. Expected probability pi of having null effect (or 1-Pi if Pi>0.5).
}
}
\details{
The model for the whole-genome regression is as follows:

\deqn{y = mu + Xb + e}

where \eqn{y} is the response variable, \eqn{mu} is the intercept, \eqn{X} is the genotypic matrix, \eqn{b} is the effect of an allele substitution (or regression coefficient) and \eqn{e} is the residual term.

Users can obtain four WGR methods out of this function: BRR (pi=0,iv=F), BayesA (pi=0,iv=T), BayesB (pi=0.01,iv=T) and BayesC (pi=0.01,iv=F). The full theoretical basis of each model is described by de los Campos et al. (2013).

Regression coefficients are updated using GSRU algorithm (Legarra and Misztal 2008). The variable selection works through the unconditional prior algorithm proposed by Kuo and Mallick (1998).
}
\value{
The EM functions returns a list with the intercept (\eqn{mu}), the regression coefficient (\eqn{b}), the fitted value (\eqn{hat}), and the estimated intraclass-correlation (\eqn{h2}).
}
\references{

de los Campos, G., Hickey, J. M., Pong-Wong, R., Daetwyler, H. D., and Calus, M. P. (2013). Whole-genome regression and prediction methods applied to plant and animal breeding. Genetics, 193(2), 327-345.

Kuo, L., & Mallick, B. (1998). Variable selection for regression models. Sankhya: The Indian Journal of Statistics, Series B, 65-81.

Legarra, A., & Misztal, I. (2008). Technical note: Computing strategies in genome-wide selection. Journal of dairy science, 91(1), 360-366.

}
\author{
Alencar Xavier
}
\examples{

data(tpod)

# BLUP
BRR = emRR(y,gen)
cor(y,BRR$hat)

# BayesA
BayesA = emBA(y,gen)
cor(y,BayesA$hat)

# BayesB
BayesB = emBB(y,gen)
cor(y,BayesB$hat)

# BayesC
BayesC = emBC(y,gen)
cor(y,BayesC$hat)
}