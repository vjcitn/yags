\name{m2LG}
\alias{m2LG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
compute a gaussian log likelihood criterion for a GEE fit
}
\description{
compute a gaussian log likelihood criterion for a GEE fit
}
\usage{
m2LG(gmod, response, x, id, tim, invlink = function(x) x, hetfac = function(m) rep(1, length(m)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gmod}{
a gee or yagsResult instance
}
  \item{response}{
the response data vector for the instance
}
  \item{x}{
the design matrix (will typically need an intercept column) for the instance
}
  \item{id}{
the cluster discriminator
}
  \item{tim}{
the time coordinates when relevant for time series models
}
  \item{invlink}{
transformation of Xbeta to generate mean response
}
  \item{hetfac}{
function of estimated mean to multiply variance function for heteroskedasticity
}
}
\details{
Note that yags now computes the m2LG value as a matter of course
for several working structures.  This code is simply called
after yags has converged before returning the yagsResult instance.
}
\value{
a scalar
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
VJ Carey
}
%\note{
%%%  ~~further notes~~
%}

\examples{
set.seed(2340)
y = rnorm(200,10)
x = runif(200)
id = as.numeric(rep(1:50, each=4))
require(gee)
require(yags)
f1 = gee(y~x, id=id, corstr="exchangeable")
f2 = yags(y~x, id=id, corstr="exchangeable", alphainit=.2)
require(nlme)
f3 = gls(y~x, cor=corCompSymm(.2, form=~1|id), method="ML")
m2LG(f1, y, cbind(1,x), id, id ) # tim variable irrelevant
m2LG(f2, y, cbind(1,x), id, id ) # tim variable irrelevant
-2*f3$logLik # for ML estimate of correlation
f4 = gls(y~x, cor=corCompSymm(f2@alpha, form=~1|id,fixed=TRUE), method="ML")
-2*f4$logLik
# relate m2LG to deviance in simple case
f5 = yags(y~x, id=id, corstr="independence")
summary(f5)
f5@m2LG
f5@phi*(f5@m2LG-200*log(2*pi*f5@phi))
f6 = glm(y~x)
f6
# generate heteroskedastic data and check
y2 = rnorm(200,25,(1+2*x)^2)
f7 = yags(y2~x, id=id, corstr="independence", family=quasi(var="mu^2"), betainit=c(25,.03))
f7@m2LG
f7@phi*(f7@m2LG-200*log(2*pi*f7@phi))
f8 = glm(y2~x, family=quasi(var="mu^2"))
summary(f8)
}
\keyword{ models }
