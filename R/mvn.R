
gauPatCorLL <- function(y, x, tim, beta, var, id, clusCorFun, corPar, invlink=function(x)x,
    hetfac = function(m)rep(1,length(m)))
 {
#
# use hetfac = function(m) m for linear het, function(m)m^2 for quadratic
# actually -2LL
# clusCorFun has form function( par, coord )
 require(mvtnorm)
 ys <- split(y,id)
 inds <- split(1:length(y),id)
 ts <- split(tim, id)
 ll <- 0
 for (i in 1:length(inds))
  {
  mn = invlink(x[inds[[i]],,drop=FALSE] %*% beta)
  hf = hetfac(mn)
  HM = diag(as.numeric(hf),length(hf), length(hf))
  sigi <- var*sqrt(HM)%*%clusCorFun(  corPar, ts[[i]] )%*%sqrt(HM)
  dimnames(sigi) = NULL
  ll <- ll+dmvnorm( ys[[i]], invlink(x[inds[[i]],,drop=FALSE] %*% beta), sigi, log=TRUE )
  }
 -2*ll
 }

m2LG = function(gmod,response,x,id,tim,invlink=function(x)x,hetfac=function(m)rep(1,length(m))) {
     #
     # some pattern definitions
     exchCor = function (par, coord) 
     {
         n <- length(coord)
         o <- matrix(par, nr = n, nc = n)
         diag(o) <- 1
         o
     }
     ar1Cor = function (par, coord) 
     par^as.matrix(dist(coord))
     #
     # end
     #
  if (inherits(gmod, "gee")) {
    beta = gmod$coef
    var = gmod$scale
    rho = gmod$working.corr[1,2]
    wmod = gmod$model$corstr
    if (wmod == "Exchangeable") ccfun = exchCor
    else if (wmod == "AR-M") ccfun = ar1Cor
    else if (wmod == "Independent") ccfun = function( par, coord) diag(length(coord))
    else {
      warning("m2LG only defined for exch, AR-M[M=1], independent working structures.  Returning NA for m2LG")
      return(as(NA, "numeric"))
      }
    }
  else if (inherits(gmod, "yagsResult")) {
    beta = coef(gmod)
    var = gmod@phi
    rho = gmod@alpha
    wmod = gmod@corstruct.tag
    if (wmod == "exchangeable") ccfun = exchCor
    else if (wmod == "ar1" | wmod == "UQ.fom") ccfun = ar1Cor
    else if (wmod == "independence") ccfun = function( par, coord) diag(length(coord))
    else {
      warning("m2LG only defined for exch, ar1, independence working structures.  Returning NA for m2LG")
      return(as(NA, "numeric"))
      }
   }
  else stop("requires either gee or yags result")
  return(gauPatCorLL(response, x, tim, beta, var, id, ccfun, rho, invlink, hetfac))
}

