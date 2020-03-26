setClass("simout", representation(df="data.frame", trueCor="character", varfun4sim="function", mvgenfun="character", empcor="matrix", respmat="matrix"))

hetSim = function(NCLUST=50, CLUSTSIZE=4, RHO=.5, alpinit=.25, trueCor="ind",
  working = c("exchangeable", "ar1", "independence"),
  varfun4sim = function(x) rep(1,length(x)), mvgen = mvrnorm, logscale=FALSE) {
    require(yags)
    require(MASS)
    require(nlme)
    XL = list()
    
    TX = function(x)x
    if (logscale) TX=log
    for (i in 1:NCLUST) {
     XL[[i]] = rep(NA, CLUSTSIZE)
     for (j in 1:CLUSTSIZE) {
      XL[[i]][j] = TX(runif(1,j,j+1))
     }
    }
    
    X = unlist(XL)
    MU = t(matrix(X, nr=4))
    R = diag(CLUSTSIZE)
    if (trueCor == "exch") {
        R[] = RHO
        diag(R) = 1
        }
    else if (trueCor == "ar1") {
        R = as.matrix(dist(1:CLUSTSIZE))
        R = RHO^R
    }
    else if (trueCor != "ind") {
        stop("true must be ind, exch, or ar1")
    }
    
    Y = matrix(NA, nc=CLUSTSIZE, nr=nrow(MU))
    if (logscale) TX=exp
    for (i in 1:nrow(MU)) {
     A = sqrt(diag(varfun4sim(TX(MU[i,]))))
     Y[i,] = mvgen(1, TX(MU[i,]), A %*% R %*% A)
    }
    empcor = cor(Y)
    
    met = as.numeric(rep(1:CLUSTSIZE, NCLUST))
    df = data.frame(y=as.numeric(t(Y)), x=X, id = rep(1:NCLUST, each=CLUSTSIZE), met=met)
    new("simout", df=df, trueCor=trueCor, varfun4sim=varfun4sim, mvgenfun=
       deparse(substitute(mvgen)), empcor=empcor, respmat=Y)
}

setMethod("show", "simout", function(object) {
 cat("simout instance, dim = ", dim(object@df), "\n")
 cat("mvgenfun: ", object@mvgenfun, "\n")
 cat("Corr:\n")
 print(object@trueCor)
 cat("var:\n")
 print(object@varfun4sim)
})

getCriteria = function(simout, 
       working = c("exchangeable", "ar1", "independence"), 
       yagsfam = gaussian(), alpinit=.5, yagsformula=y~x,
       glsformula=y~x, glsVarFunc=varIdent(form=~1), ... ) {
    df = simout@df
    trueCor = simout@trueCor
    themet <<- df$met
    doy = function(z) yags(yagsformula, id=id, data=df, corstr=z, alphainit=alpinit,
         cor.met=themet, family=yagsfam, ...)
    ans = lapply(working, doy)
    names(ans) = working
    exchGLS = gls(glsformula, corr=corCompSymm(alpinit, form=~1|id), weights=
      glsVarFunc, data=df, method="ML")
    ans[["exchGLS"]] = exchGLS
    ar1GLS = gls(glsformula, corr=corAR1(alpinit, form=~1|id), weights=
      glsVarFunc, data=df, method="ML")
    ans[["ar1GLS"]] = ar1GLS
    identGLS = gls(glsformula, corr=corIdent(form=~1|id), weights=
      glsVarFunc, data=df, method="ML")
    ans[["identGLS"]] = identGLS
    ans[["data"]] = df
    ans[["true"]] = simout@trueCor
    ans[["varfun4sim"]] = simout@varfun4sim
    ans
}

hetsim.control = function(NCLUST=50, CLUSTSIZE=4, RHO=.5, alpinit=.25, trueCor="ind",
  working = c("exchangeable", "ar1", "independence"),
  varfun4sim = function(x) rep(1,length(x)), mvgen = mvrnorm, logscale=FALSE) {
  list(NCLUST=NCLUST, CLUSTSIZE=CLUSTSIZE, RHO=RHO,
    alpinit=alpinit, trueCor=trueCor, working=working, varfun4sim=varfun4sim,
    mvgen=mvgen, logscale=logscale)
}

getCriteria.control = function(working=c("exchangeable", "ar1", "independence"),
       yagsfam = gaussian(), alpinit=.5, yagsformula=y~x,
       glsformula=y~x, glsVarFunc=varIdent(form=~1), ...) {
     list(working=working, yagsfam=yagsfam, alpinit=alpinit, yagsformula=yagsformula,
       glsformula=glsformula, glsVarFunc=glsVarFunc, ...)
}

setClass("hetEval", contains="list")
setMethod("show", "hetEval", function(object) {
 cat("hetEval object.  There were ", length(object[[1]]), "runs.\n")
 cat("components of runs:\n")
 print(names(object[[1]][[1]]))
})
    
hetEvalOLD = function(NSIM=100, 
   hetSimParms=hetsim.control(varfun4sim = function(x)1:4, mvgen=mvrnorm),
   naiveCritParms = getCriteria.control(yagsfam=gaussian(), yagsformula=y~x,
     glsformula=y~x),
   oracleCritParms = getCriteria.control(yagsfam=quasi(var=mu), yagsformula=y~x,
     glsformula=y~x, glsVarFunc=varFunc(~x)),
     doer = get("%do%")) {
 require(foreach)
 outnaive = list()
 outoracle = list()
 runif(1)
 assign("%doer%", doer)
 foreach (i = 1:NSIM) %doer% {
   data = do.call("hetSim", hetSimParms)
   nCritParms = c(simout=data, naiveCritParms)
   outnaive[[i]] = do.call("getCriteria", nCritParms)
   oCritParms = c(simout=data, oracleCritParms)
   outoracle[[i]] = do.call("getCriteria", oCritParms)
 }
 new("hetEval", list(naive=outnaive, oracle=outoracle))
}

hetEval = function(NSIM=100, 
   hetSimParms=hetsim.control(varfun4sim = function(x)1:4, mvgen=mvrnorm),
   naiveCritParms = getCriteria.control(yagsfam=gaussian(), yagsformula=y~x,
     glsformula=y~x),
   oracleCritParms = getCriteria.control(yagsfam=quasi(var=mu), yagsformula=y~x,
     glsformula=y~x, glsVarFunc=varFunc(~x)), applier=lapply) {
 outboth = applier(1:NSIM, function(x) {
   data = do.call("hetSim", hetSimParms)
   nCritParms = c(simout=data, naiveCritParms)
   ansn = do.call("getCriteria", nCritParms) 
   oCritParms = c(simout=data, oracleCritParms)
   anso = do.call("getCriteria", oCritParms)
   gc()
   list(naive=ansn, oracle=anso)
})
 outnaive = lapply(outboth, function(x)x[[1]])
 outoracle = lapply(outboth, function(x)x[[2]])
 new("hetEval", list(naive=outnaive, oracle=outoracle))
}

 
MSEslo = function(x,trueval=1,mod="exchangeable") {
    if (mod == "exch") mod = "exchangeable"
    ll = lapply(x, SEslo, trueval, mod); 
     bad = sapply(ll, function(x) any(is.na(x)))
     if (any(bad)) warning(paste(sum(bad), "NAs in MSEslo"))
    sapply(ll, mean,na.rm=TRUE)}

SEslo = function(z,trueVal=1,mod="exchangeable") {
    sapply(z, function(w)(coef(w[[mod]])[2]-trueVal)^2) }


del1 = function (x) {
    tmp = sapply(c("rj1", "rj2"), function(z) slot(yags.adeqReport(x),
        z))
    tmp[2] - 2 * tmp[1] + 1
}

del2 = function (x) {
    tmp = sapply(c("rj1", "rj2"), function(z) slot(yags.adeqReport(x),
        z))
    sum(log(tmp)^2)
}



allcrit = function (x) 
{
    #z = try(c(Lg.ind = x$iden$log, Lg.exch = x$exchG$log, Lg.ar1 = x$ar1G$log, 
    z = try(c(Lg.ind = x$indepen@m2LG/-2., Lg.exch = x$exchange@m2LG/-2., 
                   Lg.ar1 = x[["ar1"]]@m2LG/-2.,
        del1.ind = del1(x$indepen), del1.exch = del1(x$exchange), 
        del1.ar1 = del1(x[["ar1"]]), pan.ind = x$indepen@pan.aic, 
        pan.exch = x$exchange@pan.aic, pan.ar1 = x[["ar1"]]@pan.aic, 
        del2.ind = del2(x$indepen), del2.exch = del2(x$exchange), 
        del2.ar1 = del2(x[["ar1"]])))
    if (inherits(z, "try-error")) {
        warning("some model's criteria were unavailable possibly due to nonconvergence")
        return(rep(NA, 12))
        }
    else return(z)
}


picks = function(x) {
  tru = x[["true"]]
  cr = allcrit(x)
  lpic = which.max(cr[1:3])
  d1pic = which.min(cr[4:6])
  qpic = which.min(cr[7:9])
  d2pic = which.min(cr[10:12])
  pn = c(names(lpic), names(d1pic), names(qpic), names(d2pic))
  ok = rep(FALSE, 4)
  names(ok) = c("Lg", "del1", "qic", "del2")
  hit = grep(tru, pn)
  if (length(hit)>0) ok[hit] = TRUE
  ok
}


allcrit2 = function (x, ind) 
{
    xx = x
    x = xx$naive[[ind]]
    naicrit.nai = c(Lg.ind.nai = x$iden$log, Lg.exch.nai = x$exchG$log, 
        Lg.ar1.nai = x$ar1G$log, del1.ind.nai = del1(x$indepen), 
        del1.exch.nai = del1(x$exchange), del1.ar1.nai = del1(x[["ar1"]]), 
        pan.ind.nai = x$indepen@pan.aic, pan.exch.nai = x$exchange@pan.aic, 
        pan.ar1.nai = x[["ar1"]]@pan.aic, del2.ind.nai = del2(x$indepen), 
        del2.exch.nai = del2(x$exchange), del2.ar1.nai = del2(x[["ar1"]]))
    x = xx$oracle[[ind]]
    naicrit.ora = c(Lg.ind.ora = x$iden$log, Lg.exch.ora = x$exchG$log, 
        Lg.ar1.ora = x$ar1G$log, del1.ind.ora = del1(x$indepen), 
        del1.exch.ora = del1(x$exchange), del1.ar1.ora = del1(x[["ar1"]]), 
        pan.ind.ora = x$indepen@pan.aic, pan.exch.ora = x$exchange@pan.aic, 
        pan.ar1.ora = x[["ar1"]]@pan.aic, del2.ind.ora = del2(x$indepen), 
        del2.exch.ora = del2(x$exchange), del2.ar1.ora = del2(x[["ar1"]]))
    c(naicrit.nai, naicrit.ora)
}

picks2 = function (x, ind, hettag="ora") 
{
    tru = x[[1]][[ind]][["true"]]
    cr = allcrit2(x, ind)
    lpic = which.max(cr[c(1:3, 13:15)])
    d1pic = which.min(cr[c(4:6, 16:18)])
    qpic = which.min(cr[c(7:9, 19:21)])
    d2pic = which.min(cr[c(10:12, 22:24)])
    pn = c(names(lpic), names(d1pic), names(qpic), names(d2pic))
    ok = oraneed = rep(FALSE, 4)
    names(oraneed) = names(ok) = c("Lg", "del1", "qic", "del2")
    hit = grep(tru, pn)
    orahit = grep(hettag, pn)
    if (length(hit) > 0) 
        ok[hit] = TRUE
    if (length(orahit) > 0) 
        oraneed[orahit] = TRUE
    apply(data.frame(critpicks = ok, usesOra=oraneed),1,all)
}

summarize = function (outcomp, tag = c("ora", "nai")[1]) 
{
    tmp = sapply(1:length(outcomp[[1]]), function(x) picks2(outcomp, 
        x, tag))
    apply(tmp, 1, mean)
}

getPicks = function(struc, varmod="f", corstr="exch") {
 rr = lapply(struc, sapply, allcrit)
 rownames(rr[[1]]) = paste("n", rownames(rr[[1]]), sep=".")
 rownames(rr[[2]]) = paste("f", rownames(rr[[2]]), sep=".")
 arr = rbind(rr[[1]], rr[[2]])
 L = arr[c(1:3, 13:15), ]
 D1 = arr[c(4:6, 16:18), ]
 P = arr[c(4:6, 16:18)+3, ]
 D2 = arr[c(4:6, 16:18)+6, ]
 LPICKS = rownames(L)[apply(L,2,which.max)]
 D1PICKS = rownames(D1)[apply(D1,2,which.min)]
 D2PICKS = rownames(D2)[apply(D2,2,which.min)]
 PPICKS = rownames(P)[apply(P,2,which.min)]
   probRight = function(vec) length(grep(corstr, grep(paste("^",
        varmod, sep = ""), vec, value = TRUE)))/length(vec)

 Lprob = probRight(LPICKS)
 D1prob = probRight(D1PICKS)
 D2prob = probRight(D2PICKS)
 Pprob = probRight(PPICKS)
 c(Lprob=Lprob, D1prob=D1prob, D2prob=D2prob, Pprob=Pprob)
}

