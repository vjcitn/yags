rar1poi = function(n=100, muvec = c(2,3,4), cov=.7^as.matrix(dist(1:3))) {
        ni = length(muvec)
        if (ni < 2) stop("must have at least 2 elements in muvec")
        if (nrow(cov) != ni) stop("muvec and cov not conformant")
        rho = cov[1,2]
        out <- matrix(NA, nr = n, nc = length(muvec))
        out[, 1] <- rpois(n, muvec[1])
        for(j in 2:length(muvec))
                out[, j] <- rpois(n, muvec[j] + (rho * muvec[j])/muvec[j - 1] *
                        (out[, j - 1] - muvec[j - 1]))
        out
}
