'qicyags'<-function (obj, response, indep = NULL, family = gaussian, ql = c("Gamma", 
    "binomial", "gaussian", "inverse.gaussian", "poisson"))
{


# This script calculates the QIC(R) according to the book by Hardin and Hilbe (The Stata macro).

    y <- response
    mu <- obj@fitted.values

#Calculate the sum of quasi-log-likelihood, i.e., the first term of QIC(R)

    if (ql[1] == "Gamma") 
        q.mu <- (-1) * sum(y/mu + log(mu))
    if (ql[1] == "binomial") 
        q.mu <- sum(y * log(mu/(1 - mu)) + log(1 - mu))
    if (ql[1] == "gaussian") 
        q.mu <- (-1/2) * sum((y - mu)^2)
    if (ql[1] == "inverse.gaussian") 
        q.mu <- sum(-y/(2 * mu^2) + 1/mu)
    if (ql[1] == "poisson") 
        q.mu <- sum(y * log(mu) - mu)

#Calculate the second term of QIC(R) 

    if (qicr == T) {
        A.I.par <- (indep@naive.parmvar)/as.numeric(indep@phi)
        V.R.par <- obj@robust.parmvar
        if (length(A.I.par) == 1) 
            qic.r <- -2 * q.mu + 2 * ((as.numeric(A.I.par)^-1) * as.numeric(V.R.par))
        if (length(A.I.par) > 1) 
            qic.r <- -2 * q.mu + 2 * sum(diag((solve(A.I.par)) %*% V.R.par))
    }
    qic.r
}



