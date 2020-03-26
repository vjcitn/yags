\name{seizlong}
\alias{seizlong}
\docType{data}
\title{
data frame on seizure counts arising from a study of progabide
}
\description{
data frame on seizure counts arising from a study of progabide
}
\usage{data(seizlong)}
\format{
  A data frame with 236 observations on the following 5 variables.
  \describe{
    \item{\code{resp}}{a numeric vector of 4 successive two-week seizure counts per patient (N=59)}
    \item{\code{lbas}}{a numeric vector of log of the baseline two-week count}
    \item{\code{lage}}{a numeric vector log patient age}
    \item{\code{tx}}{a numeric vector 1 if progabide 0 otherwise}
    \item{\code{id}}{a numeric vector patient identifier}
  }
}
\details{
obtained from Table 2 of Thall and Vail 1990
}
\source{
Thall PF, Vail SC.  Some Covariance Models for Longitudinal Count Data
with Overdispersion.  Biometrics 1990;46:657-671.
}
%\references{
%%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(seizlong)
v4 = rep(c(0,0,0,1), 59)
f1 = yags(resp~lbas*tx+lage+v4, id=id, data=seizlong, family=poisson(), corstr="exchangeable", alphainit=.2)
f1
with(seizlong, m2LG(f1, resp, cbind(1,lbas,tx,lage,v4,lbas*tx), id=id, tim=id, invlink=function(x)exp(x), hetfac=function(m) m))
f1@del1
f2 = yags(resp~lbas*tx+lage+v4, id=id, data=seizlong, family=quasi(link="log", var="mu^2"), corstr="exchangeable", alphainit=.2)
f2
with(seizlong, m2LG(f2, resp, cbind(1,lbas,tx,lage,v4,lbas*tx), id=id, tim=id, invlink=function(x)exp(x), hetfac=function(m) m^2))
f2@del1
f3 = yags(resp~lbas*tx+lage+v4, id=id, data=seizlong, family=poisson(), corstr="ar1", cor.met=rep(1:4,59), alphainit=.2)
with(seizlong, m2LG(f3, resp, cbind(1,lbas,tx,lage,v4,lbas*tx), id=id, tim=rep(1:4,59), invlink=function(x)exp(x), hetfac=function(m) m))
f3@del1
f4 = yags(resp~lbas*tx+lage+v4, id=id, data=seizlong, family=quasi(link="log", var="mu^2"), corstr="ar1", cor.met=rep(1:4,59), alphainit=.2)
with(seizlong, m2LG(f4, resp, cbind(1,lbas,tx,lage,v4,lbas*tx), id=id, 
       tim=rep(1:4,59), invlink=function(x)exp(x), hetfac=function(m) m^2))
f4@del1
}
\keyword{datasets}
