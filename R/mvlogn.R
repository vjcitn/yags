

rmvlnorm = function(n=100, lmu = c(2,3,4), lcov=diag(c(2,3,4),3,3)) {
#
# algorithm for correlated lognormals, from Wang and Carey 2004 JASA
# Sept 2004  lmu is the intended mean value for the lognormal response
# lcov is the intended covariance matrix for the lognormal response
# returns an n x p matrix, where p is the length of lmu
#
# NB -- to recover the means, one does not take
#
 require(MASS)
 cplus = function(x) outer(x,x,"+")
 dlc = dim(lcov)
 if (!all(dlc == length(lmu))) stop("row dim of lcov and length lmu must coincide")
 S = diag(diag(log(1+diag(lcov)/lmu^2)))
 eta = log(lmu) - S/2
 sig = log(lcov + exp(cplus(eta) + .5*cplus(S))) - (cplus(eta) + .5*cplus(S))
# list(eta, sig)
 exp(mvrnorm(n, eta, sig))
}
