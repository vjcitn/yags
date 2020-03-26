adap = function (struc, varmod = "f", corstr = "exch") 
{
    rr = lapply(struc, sapply, allcrit)
    rownames(rr[[1]]) = paste("n", rownames(rr[[1]]), sep = ".")
    rownames(rr[[2]]) = paste("f", rownames(rr[[2]]), sep = ".")
    arr = rbind(rr[[1]], rr[[2]])
    L = arr[c(1:3, 13:15), ]
    D1 = arr[c(4:6, 16:18), ]
    P = arr[c(4:6, 16:18) + 3, ]
    D2 = arr[c(4:6, 16:18) + 6, ]
    LPICKS = rownames(L)[apply(L, 2, which.max)]
    D1PICKS = rownames(D1)[apply(D1, 2, which.min)]
    D2PICKS = rownames(D2)[apply(D2, 2, which.min)]
    PPICKS = rownames(P)[apply(P, 2, which.min)]
    allco = function (x) sapply(x, function(x) c(ex = coef(x$exchange)[2], ar1 = coef(x$ar1)[2], 
         ind = coef(x$indep)[2]))
    plist = list(L = LPICKS, D1 = D1PICKS, D2 = D2PICKS, P = PPICKS, 
        CON = allco(struc$naive), COF = allco(struc$oracle))
    ainds = sapply(plist[1:4], function(x) grep("ar1", x))
    einds = sapply(plist[1:4], function(x) grep("exch", x))
    iinds = sapply(plist[1:4], function(x) grep("ind", x))
    useN = lapply(plist[1:4], function(x) grep("^n",x))  # sapply -> lapply Mar 30 2011
    useF = sapply(plist[1:4], function(x) grep("^f",x))
    ncoef = ncol(plist$CON)
    Ladap = rep(NA, ncoef)
    D1adap = rep(NA, ncoef)
    D2adap = rep(NA, ncoef)
    Padap = rep(NA, ncoef)
# following updates a vector of adaptively chosen coefficients
# if run on all cormods gives a full vector of chosen coeffs
    pullco = function( adaco, crittype = "L", cormod="ar1" ) {
      if (cormod == "ar1") vindset = ainds
      else if (cormod == "ex") vindset = einds
      else if (cormod == "ind") vindset = iinds
      if (length(vindset[[crittype]]) > 0) {
          doN = intersect(vindset[[crittype]], useN[[crittype]])
          if (length(doN > 0))
                  adaco[ doN ] = plist$CON[cormod, doN ]
          doF = intersect(vindset[[crittype]], useF[[crittype]])
          if (length(doF > 0))
                  adaco[ doF  ] = plist$COF[cormod, doF ]
          }
      adaco
    }
    Ladap = pullco(Ladap, "L", "ar1")
    Ladap = pullco(Ladap, "L", "ex")
    Ladap = pullco(Ladap, "L", "ind")
    D1adap = pullco(D1adap, "D1", "ar1")
    D1adap = pullco(D1adap, "D1", "ex")
    D1adap = pullco(D1adap, "D1", "ind")
    D2adap = pullco(D2adap, "D2", "ar1")
    D2adap = pullco(D2adap, "D2", "ex")
    D2adap = pullco(D2adap, "D2", "ind")
    Padap = pullco(Padap, "P", "ar1")
    Padap = pullco(Padap, "P", "ex")
    Padap = pullco(Padap, "P", "ind")
    list(Ladap=Ladap, D1adap=D1adap, D2adap=D2adap, Padap=Padap, plist=plist)
}
