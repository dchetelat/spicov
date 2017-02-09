\name{spiCov}
\alias{spiCov}
\title{
    Compute the Covariance Estimator
}
\description{
    Computes the high-dimensional spiked covariance estimator, given a sample covariance
    matrix based on the data. The spiked rank can be automatically selected (default) or
    specified.
}
\usage{
    spiCov(S,n,r=NULL,verbose=FALSE)
}
\arguments{
    \item{S}{
	a sample covariance matrix. Must be symmetric positive-definite.
    }
    \item{n}{
        the sample size. Must be greater than the dimension of S.
    }
    \item{r}{
	the spiked rank. If unspecified or NULL, the rank is selected based on an unbiased risk
	estimation procedure.
    }
    \item{verbose}{
	logical. Prints additional information.
    }
}
\details{
    The spectral decomposition is done through Armadillo.
}
\value{
    The covariance estimator, a symmetric positive-definite matrix of the same dimensions as S.
}
\author{
    Didier Chételat (dc623@cornell.edu)
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
\keyword{file}
