\name{spicov-package}
\alias{spicov}
\docType{package}
\title{
High-dimensional Spiked Covariance Matrix Estimation
}
\description{
Computes a high-dimensional spiked covariance matrix estimate by correcting the eigenvalues of a sample covariance matrix. The spiked rank can be specified, or selected automatically (default). See Ch\'etelat and Wells (2014) (arxiv.org/abs/1408.6440).
}
\details{
\tabular{ll}{
Package: \tab spicov\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014-06-11\cr
License: \tab GPL-2\cr
}
There is currently only one function, \code{\link{spiCov}} which computes the estimator. 
}
\author{
Maintainer: Didier Chételat <dc623@cornell.edu>
}
\references{
D. Ch\'etelat and M.T. Wells (2014). Noise estimation in the spiked covariance model. Preprint arXiv:1408.6440 [stat.TH].
}
\examples{
library(spicov)
library(MASS)

p=20;
c=0.5;
n=floor(p/c);
Sig=diag(c(10,9,8,7,6,rep(1,p-5)));
X=mvrnorm(n,rep(0,p),Sig);
S=t(X)\%*\%X/n;

hS=spiCov(S,n)
}
\keyword{ package }
