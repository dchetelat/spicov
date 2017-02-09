pkgname <- "spicov"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('spicov')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("spiCov")
### * spiCov

flush(stderr()); flush(stdout())

### Name: spiCov
### Title: Compute the Covariance Estimator
### Aliases: spiCov
### Keywords: file

### ** Examples

library(spicov)
library(MASS)

p=20;
c=0.5;
n=floor(p/c);
Sig=diag(c(10,9,8,7,6,rep(1,p-5)));
X=mvrnorm(n,rep(0,p),Sig);
S=t(X)%*%X/n;

hS=spiCov(S,n)



cleanEx()
nameEx("spicov-package")
### * spicov-package

flush(stderr()); flush(stdout())

### Name: spicov-package
### Title: High-dimensional Spiked Covariance Matrix Estimation
### Aliases: spicov
### Keywords: package

### ** Examples

library(spicov)
library(MASS)

p=20;
c=0.5;
n=floor(p/c);
Sig=diag(c(10,9,8,7,6,rep(1,p-5)));
X=mvrnorm(n,rep(0,p),Sig);
S=t(X)%*%X/n;

hS=spiCov(S,n)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
