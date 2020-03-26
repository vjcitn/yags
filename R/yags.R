yags.make.libu <- function(lib.loc=.lib.loc,cxx="g++ -c") {
acat <- function(x) {cat(x,file="user_wcor.cc",append=T); 
  cat("\n",file="user_wcor.cc",append=T); invisible(0)}
cat('#include "MC++.h"',file="user_wcor.cc")
cat('\n',file="user_wcor.cc",append=T)
acat('#include "MC++class.h"')
acat('matrix alpfun_user(matrix PRin, matrix ID, matrix TIMin, ')
acat('       double phi, int p, matrix alpin,')
acat('       double atol, int amaxit)')
acat('{return mat11(0.5);}')
acat('')
acat('matrix wcorinv_user(matrix alp, int ni, matrix tim) {')
acat('// alp, tim included for fixed signature')
acat(' return ident(ni);')
acat(' }')
inclpathsw <- paste("-I",lib.loc,"/yags/scripts",sep="")
system(paste(cxx,inclpathsw,"user_wcor.cc"))
system(paste("R SHLIB user_wcor.o"))
system("mv user_wcor.so libu.so")
invisible(0)
}

.First.lib <- function(lib, pkg) {
   library.dynam("yags", pkg, lib) 
   }

setOldClass("family")
setClass("yagsResult",
 representation(coefficients="numeric", coefnames="character",
	naive.parmvar="matrix", robust.parmvar="matrix",
	corstruct.tag="character", alpha="numeric",
	phi="numeric", linear.predictors="numeric",
	fitted.values="numeric", residuals="numeric",
	iter="numeric", family="family",
	rank="numeric", errorcode="numeric",
	sealp="numeric", qls="numeric", pan.aic="numeric",
	sealp.OK="logical", varnames="character",
	n="numeric", nclus="numeric", maxni="numeric",
        wcovmat="matrix", wcormat="matrix", m2LG="numeric", del1="numeric", del2="numeric",
        critcoefs = "list",
	Call="call") )

setGeneric("coef", function(object, ...)standardGeneric("coef"))
setMethod("coef", "yagsResult", function(object, ...) {
 object@coefficients
})

setGeneric("zTable", function(object, ...)standardGeneric("zTable"))
setMethod("zTable", "yagsResult", function(object, ...) {
  yags.glmReport(object)
})

setGeneric("alpha", function(object, ...) standardGeneric("alpha"))
setMethod("alpha", "yagsResult", function(object, ...) {
 object@alpha
})

setGeneric("sealp", function(object, ...) standardGeneric("sealp"))
setMethod("sealp", "yagsResult", function(object, ...) {
 object@sealp
})

setGeneric("wcor", function(object, ...) standardGeneric("wcor"))
setMethod("wcor", "yagsResult", function(object, ...) {
 alp = object@alpha
 maxni = object@maxni
 if (maxni == 1) return(1)
 if (maxni > 10) {
     warning("max cluster size > 10; trimming cor. report to 10x10")
     maxni = 10
     }
 cortag = object@corstruct.tag
 if (cortag == "independence") return(diag(maxni))
 if (cortag == "exchangeable") {
    ans = matrix(alp, maxni, maxni)
    diag(ans) = 1
    return(ans)
    }
 if (cortag %in% c("ar1", "UQ.fom", "UJ.fom") ) {
    ans = as.matrix(dist(1:maxni))
    return(alp^ans)
    }
 if (cortag == "unstructured") {
    if (nrow(object@wcormat)>10)
       return(object@wcormat[1:10,1:10])
    else return(object@wcormat)
    }
 else stop(paste("no wcor report for model with corstruct", cortag))
})

setClass("yagsAdeq",
 representation(qls="numeric", pan.aic="numeric",
        rj1="numeric", rj2="numeric"))

print.yagsResult <- function(x,...)
 {
 cat("YAGS (yet another GEE solver) $Revision: 5.6 $\n\n")
 cat("Call:\n")
 print(slot(x,"Call"))
 cat("Regression estimates:\n\n")
 print(yags.glmReport(x))
 cat("\n")
 print(yags.wcorReport(x))
 cat("\n")
 cat(show(yags.adeqReport(x)))
 cat("\n")
 cat("yags/R: $Id: yags.R,v 5.6 2008/03/14 18:39:13 stvjc Exp $\n")
 invisible(0)
 }

yags.glmReport <- function(x,...)
 {
 coef <- slot(x,"coefficients")
 nse <- sqrt(diag(slot(x,"naive.parmvar")))
 nz <- coef/nse
 rse <- sqrt(diag(slot(x,"robust.parmvar")))
 rz <- coef/rse
 omat <- cbind(coef,nse,nz,rse,rz)
 dimnames(omat) <- list(slot(x,"varnames"),
	c("est.","naive s.e.", "naive z",
         "sand. s.e.", "sand. z"))
 omat
 }

yags.wcorReport <- function(x,...)
 {
 alpha <- slot(x,"alpha")
 maxni <- slot(x,"maxni")
 cortag <- slot(x,"corstruct.tag")
 sealpOK <- slot(x,"sealp.OK")
 sealp <- slot(x, "sealp")
 cat(paste("Working correlation model:", cortag))
 if (cortag=="unstructured") 
     {
     cat("\nworking covariance:\n") 
     print(round(slot(x,"wcovmat"),4))
     cat("\nworking correlation:\n") 
     print(round(slot(x,"wcormat"),4))
     }
 else cat(paste("\nalpha est:", round(alpha,4)))
 cat(ifelse(sealpOK,paste("; s.e. alpha",round(sealp,4))," "))
 cat("\n")
 #invisible(list(alpha=alpha,sealp=sealp))
 invisible(NULL)
 }

yags.adeqReport <- function(x,...)
 {
 qls <- slot(x,"qls")
 pan.aic <- slot(x,"pan.aic")
 H0 <- slot(x,"naive.parmvar")
 H1 <- slot(x,"robust.parmvar")
 H <- solve(H0,tol=1e-12)%*%H1
 c1 <- mean(eigen(H)$values)
 c2 <- mean(eigen(H%*%H)$values)
 RJ.crit <- c(c1=c1,c2=c2)
 new("yagsAdeq", qls=qls, pan.aic=pan.aic, 
      rj1=RJ.crit[1], rj2=RJ.crit[2])
 }

setGeneric("adeq", function(x) standardGeneric("adeq"))
setMethod("adeq", "yagsResult", function(x){
 yags.adeqReport(x)})

setMethod("show", "yagsAdeq", function(object) {
 cat(paste("Pan QIC(R):", round(object@pan.aic,3),"\nQLS:",
  round(object@qls,3),
  "\nRotnitzky-Jewell:", paste(round(c(object@rj1,object@rj2),3),collapse=", ")))})
 
setMethod("show", "yagsResult", 
     function(object) print.yagsResult(object))
  
yags.control <- function(maxiter=15, tol=.0001, verbose=FALSE,
                  Ua.maxit=20, Ua.tol=0.001, 
                  Ua.gridlo=0.1, Ua.gridhi=0.9,
                  Ua.gridnpts = 10, Ua.secantdel = 0.01,
                  fixscale=FALSE)
    list(maxiter=maxiter, tol=tol, verbose=verbose,
                  Ua.maxit=Ua.maxit, Ua.tol=Ua.tol, 
                  Ua.gridlo=Ua.gridlo, Ua.gridhi=Ua.gridhi,
                  Ua.gridnpts = Ua.gridnpts, 
                  Ua.secantdel = Ua.secantdel,
                  fixscale=fixscale)

yags <- function(formula, id, 
        cor.met=NULL, family=gaussian(),
        corstruct="independence", control=yags.control(), 
        weights=NULL, betainit=NULL, alphainit=NULL, data=list(), subset=NULL,
          allcrit=FALSE, lhetfam=NULL, qhetfam=NULL , icritalp=NULL, critar1tag="ar1")
{
#
# $Header: /udd/stvjc/VCROOT/yags/R/yags.R,v 5.6 2008/03/14 18:39:13 stvjc Exp $
#
    headstring <- "$Header"
#
# need id for error checking
#
	m <- match.call(expand=FALSE)
        m$family <- m$corstruct <- m$control <- m$betainit <-
            m$alphainit <- m$allcrit <- m$lhetfam <- m$qhetfam <- m$icritalp <- m$critar1tag <- NULL
        m[[1]] <- as.name("model.frame")
        TMP <- eval(m, sys.parent())
	id <- TMP[["(id)"]]
	cor.met <- TMP[["(cor.met)"]]
	nclus <- length(sid <- split(id,id))
        lsid <- sapply(sid,length)
        maxni <- max(lsid)
        if (corstruct == "unstructured") {
            if(is.null(cor.met)) stop("for unstructured working corr., cor.met must be supplied")
            if (min(cor.met)!=0) stop("for unstructured working corr., cor.met must have zero origin, but no element of cor.met is 0")
            if ((maxcm <- max(cor.met))>maxni) warning(paste("for unstructured working corr., cor.met must range from 0 to number of sampling times\nbut max(cor.met) is", maxcm," and maxni =",maxni))
        }
        if(corstruct != "independence" &&
               is.null(alphainit))
                    stop("a q-vector must be supplied for alphainit")
	if(corstruct == "independence") alphainit <- 0
#
# defer cor.met check until reading from data is complete	
#
# need to revise id/cor.met/weights in case of subsetting
#
	m <- match.call(expand=FALSE)
        m$family <- m$corstruct <- m$control <- m$betainit <-
            m$alphainit <- m$allcrit <- m$lhetfam <- m$qhetfam <- m$icritalp <- m$critar1tag <- NULL
        m[[1]] <- as.name("model.frame")
        TMP <- eval(m, sys.parent())
	id <- TMP[["(id)"]]
	cor.met <- TMP[["(cor.met)"]]
#
# defer cor.met check until reading from data/subset is complete	
#
        if(corstruct != "independence" && corstruct != "exchangeable" &&
               is.null(cor.met))
                    stop(paste("must supply cor.met for all",
                       "but independence and exch. working corstructs"))
	weights <- TMP[["(weights)"]]
#
# -- 2009.09.27 -- why?
#
#        if(corstruct == "independence" || corstruct == "exchangeable")
#		cor.met <- rep(0, length(id))
#
# now get customary quantities
#	
	Terms <- attr(TMP, "terms")
  	y <- as.matrix(model.extract(TMP, response))
	if (is.null(weights)) weights <- rep(1,nrow(y)) else
           if (length(weights)!=nrow(y)) stop("lengths of y and weights incompatible")
	if (ncol(y)==2)
		{
		n <- apply(y,1,sum)
		y <- y[,1,drop=FALSE]/n
		weights <- 1/(weights*n)
		}
	else weights <- 1/(weights)
#
	x <- model.matrix(Terms, TMP, contrasts)
	varnames <- dimnames(x)[[2]]
	OFFSET <- model.extract(TMP,offset)
	SUBSET <<- model.extract(TMP,subset)
	if (length(OFFSET) < length(y))
		OFFSET <- rep(0.0, length(y))

#
# coerce y to vector if appropriate
#
	if(is.matrix(y) && ncol(y) == 1) y <- as.double(y)
#
#	assign("subset",subset,pos=1)
	if(is.null(betainit)) betainit <- glm(formula, family = family, 
                          subset=SUBSET, data=data)$coef
	n <- length(y)
	p <- ncol(x)
	q <- length(alphainit)
if (corstruct == "unstructured")
         {
         q <- (max(cor.met)+1)^2  # assumes lattice-valued timings
# supplied with zero origin!
         }
	phi <- 1.
	b0 <- betainit
	bout <- b0
	alpout <- rep(0,q)
	bcov.naive <- matrix(0, p, p)
	bcov.rob <- matrix(0, p, p)
	famtag <- try(tolower(family$family[1]))
        if (inherits(famtag, "try-error")) stop("use glm family as evaluated function (e.g., binomial())")
        if (famtag == "quasi") {
           ltag = family$link
           convar = FALSE
           if(isTRUE(all.equal(family$variance, quasi(variance="constant")$variance)))
		convar = TRUE
           if (!(ltag %in% c("log", "identity")))
                   stop("quasi family must have log or identity link specified")
           vtag = deparse(body(family$variance))
        #   if (ltag == "log" & ((vtag != "mu^2") & !convar)) stop("quasi family with log link in yags must have variance constant or mu^2")
           if (ltag == "log" & isTRUE(all.equal(family$variance, quasi(variance="constant")$variance))) famcode = 8
           else if (ltag == "log") famcode = 5
           else if (ltag == "identity" & vtag == "mu") famcode = 6
           else if (ltag == "identity" & vtag == "mu^2") famcode = 7
	#   else stop("quasi family in yags must have link identity and variance mu or mu^2, or link log and variance mu^2 or constant")
           else if (ltag == "log" & vtag == "mu") famcode = 9
           else stop("link and variance combination not yet handled")
           }
	else if(is.na(famcode <- match(famtag, c("gaussian", "binomial", "poisson",
		"gamma"))))
		stop("in yags, only gaussian, binomial, poisson, gamma (or quasi with restrictions, see doc) fams supported")
	CORSTR.OPTS <- c("independence", "exchangeable", 
		"UJ.fom", "ar1", "UQ.fom","unstructured", 
		"UJ.equi", "UJ.equimart", "user")
	if(is.na(corcode <- match(corstruct, CORSTR.OPTS)))
		stop(paste("only",CORSTR.OPTS,"corstrs supported"))
            sealp.OK <- FALSE
# at this point we know hardcode==T, so UJ.fom selection gets sealp
       if (corstruct=="UJ.fom" | corstruct=="UJ.equi")
            sealp.OK <- TRUE
        ctl <- control
        if (is.null(cor.met)) cor.met = rep(0,n)
	out <- .C("yags_engine",
		as.integer(n),
		as.integer(p),
		as.integer(q),
		as.double(x),
		as.double(y),
		as.double(id),
		as.double(cor.met),
		as.double(b0),
		bout = as.double(bout),
		outiter = as.integer(ctl$maxiter),
		tol = as.double(ctl$tol),
		as.integer(famcode),
		as.integer(corcode),
		as.double(alphainit),
		aout = as.double(alpout),
		phiout = as.double(phi),
		bcov.naive = as.double(bcov.naive),
		bcov.rob = as.double(bcov.rob),
		ua = as.double(0),
		duda = as.double(rep(0, q^2)),
		sum.uut = as.double(rep(0, q^2)),
		ua.tol = as.double(ctl$Ua.tol),
		uamaxit = as.integer(ctl$Ua.maxit),
		as.integer(ctl$verbose),
		as.double(ctl$Ua.gridlo),
		as.double(ctl$Ua.gridhi),
		as.integer(ctl$Ua.gridnpts),
		as.double(ctl$Ua.secantdel), as.double(weights),
                qls=as.double(0), pan.aic=as.double(0), gau.ll=as.double(0),
                as.double(OFFSET), as.integer(ctl$fixscale), PACKAGE="yags")
        wcovmat <- wcormat <- new("matrix")
        if (corstruct=="unstructured")
               {
               wcovmat <- matrix(out$aout, maxni)
               wcormat <- cov2cor(wcovmat)
               }
	tmpsig <- list(coef = out$bout, alpest = out$aout, phiest = out$phiout,
		bcov.naive = out$bcov.naive, bcov.rob = out$bcov.rob, ua = out[[
		"ua"]], duda = out$duda, sum.uut = out$sum.uut, var.alp = out$
		sum.uut/out$duda^2, qls=out$qls, pan.aic=out$pan.aic,
                wcormat = wcormat, wcovmat=wcovmat)
	#print(tmpsig)
   if (corcode != 3) tmpsig$var.alp <- NA
	tmp.eta <- as.double(x %*% tmpsig$coef)
	tmp.mu <- family$linkinv(tmp.eta)
	mpp <- function(x, p)
	matrix(x, p, p)
	Call <- match.call()
        if (out$outiter >= ctl$maxiter) warning("maximum number of yags iterations consumed, consider incrementing maxit in control parameter")
	final.out <- new("yagsResult", 
                coefficients = tmpsig$coef, naive.parmvar = mpp(
		tmpsig$bcov.naive, p), robust.parmvar = mpp(tmpsig$bcov.rob,
		p), alpha = tmpsig$alpest, phi = tmpsig$phiest, 
		linear.predictors = tmp.eta, fitted.values = as.double(tmp.mu),
		residuals = as.double(y - tmp.mu), iter = out$outiter, family
		 = family, rank = p, errorcode = 0, sealp = sqrt(tmpsig$
		var.alp), qls=tmpsig$qls, pan.aic=tmpsig$pan.aic, 
                sealp.OK=sealp.OK, varnames=varnames, n=n, nclus=nclus,
		corstruct.tag=corstruct, maxni=maxni,
                wcormat=wcormat, wcovmat=wcovmat, 
		Call=Call)
	# provisional error
        if (is.na(tmpsig$ua)) warning("U.alpha apparently unsolved, check working model adequacy")
	else if(abs(tmpsig$ua) > ctl$Ua.tol) warning(paste(
			"Alp est func not solved, results suspect, Ua=", round(
			tmpsig$ua, 5), sep = ""))
        M2LG.given = m2LG(final.out, as.double(y), x, id, cor.met, invlink=family$linkinv,
             hetfac=family$variance)
#yags <- function(formula, id, 
#        cor.met=NULL, family=gaussian(),
#        corstruct="independence", control=yags.control(), 
#        weights=NULL, betainit=NULL, alphainit=NULL, data=list(), subset=NULL)
    if (allcrit) {
        if (is.null(lhetfam)) stop("must define lhetfam with allcrit=TRUE")
        if (is.null(qhetfam)) stop("must define qhetfam with allcrit=TRUE")
        if (!is.null(subset)) stop("can't use subset with allcrit=TRUE -- please create basic data frame")
        ndata= cbind(data, id=id, cor.met=cor.met, weights=weights)
        if (is.null(icritalp)) icritalp=.5
        cat("hom...")
        indmod.hom = yags(formula, id, cor.met, family, corstruct="independence", control=control,
          weights=weights, betainit=betainit, alphainit=icritalp, data=ndata)
        excmod.hom = yags(formula, id, cor.met, family, corstruct="exchangeable", control=control,
          weights=weights, betainit=betainit, alphainit=icritalp, data=ndata)
        ar1mod.hom = yags(formula, id, cor.met, family, corstruct=critar1tag, control=control,
          weights=weights, betainit=betainit, alphainit=icritalp, data=ndata)
        cat("lin...")
        indmod.lin = yags(formula, id, cor.met, family=lhetfam, corstruct="independence", control=control,
          weights=weights, betainit=betainit, alphainit=icritalp, data=ndata)
        excmod.lin = yags(formula, id, cor.met, family=lhetfam, corstruct="exchangeable", control=control,
          weights=weights, betainit=betainit, alphainit=icritalp, data=ndata)
        ar1mod.lin = yags(formula, id, cor.met, family=lhetfam, corstruct=critar1tag, control=control,
          weights=weights, betainit=betainit, alphainit=icritalp, data=ndata)
        cat("qua...")
        indmod.qua = yags(formula, id, cor.met, family=qhetfam, corstruct="independence", control=control,
          weights=weights, betainit=betainit, alphainit=alphainit, data=ndata)
        excmod.qua = yags(formula, id, cor.met, family=qhetfam, corstruct="exchangeable", control=control,
          weights=weights, betainit=betainit, alphainit=icritalp, data=ndata)
        ar1mod.qua = yags(formula, id, cor.met, family=qhetfam, corstruct=critar1tag, control=control,
          weights=weights, betainit=betainit, alphainit=icritalp, data=ndata)
        cat("...\n")
        allm = lapply(list(indmod.hom, excmod.hom, ar1mod.hom,
                            indmod.lin, excmod.lin, ar1mod.lin,
                            indmod.qua, excmod.qua, ar1mod.qua), function(x)x@m2LG)
        alld = lapply(list(indmod.hom, excmod.hom, ar1mod.hom,
                            indmod.lin, excmod.lin, ar1mod.lin,
                            indmod.qua, excmod.qua, ar1mod.qua), del1)
        alld2 = lapply(list(indmod.hom, excmod.hom, ar1mod.hom,
                            indmod.lin, excmod.lin, ar1mod.lin,
                            indmod.qua, excmod.qua, ar1mod.qua), del2)
        allc = lapply(list(indmod.hom, excmod.hom, ar1mod.hom,
                            indmod.lin, excmod.lin, ar1mod.lin,
                            indmod.qua, excmod.qua, ar1mod.qua), coef)
        ansm = unlist(allm)
        ansd = unlist(alld)
        ansd2 = unlist(alld2)
        names(ansd2) = names(allc) = names(ansd) = names(ansm) = c("ind.hom", "exch.hom", "ar1.hom", 
                           "ind.lin", "exch.lin", "ar1.lin", 
                           "ind.qua", "exch.qua", "ar1.qua")
        final.out@m2LG = c(ansm, given=M2LG.given)
        final.out@del1 = c(ansd, given=as.numeric(del1(final.out)))
        final.out@del2 = c(ansd2, given=as.numeric(del2(final.out)))
        final.out@critcoefs = allc
        }
        else {
               final.out@m2LG = M2LG.given
               final.out@del1 = del1(final.out)
             }
        
#m2LG = function(gmod,response,x,id,tim,invlink=function(x)x,hetfac=function(m)1) {
	final.out
}

#
# some simulation assistance
#

csmat <- function(n=4, rho=.5)
 {
 o <- matrix( rho, n, n )
 diag(o) <- 1
 o
 }

ar1mat <- function(n=4, rho=.5)
 {
 o <- as.matrix(dist(1:n))
 rho^o
 }

mvnsamp <- function(n=100, m=rep(0,4), v=diag(4))
 {
 p <- length(m)
 if (nrow(v) != p) stop("v and m nonconforming")
 if (ncol(v) != p) stop("v and m nonconforming")
 id <- sort(rep(1:n,p))
 tim <- rep(1:p,n)
 x <- mvrnorm(n,m,v)
 data.frame(y=as.numeric(x),id=id,tim=tim)
 }

