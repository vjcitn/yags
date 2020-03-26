\name{yags}
\alias{yags}
\alias{trylib}
\alias{tryset}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Yet Another GEE Solver}
\description{
Solves GEEs in an extensible way, with a C++ engine.
$Header: /udd/stvjc/VCROOT/yags/man/yags.Rd,v 5.1 2007/12/11 16:57:19 stvjc Exp $
}
\usage{
yags(formula, id, cor.met=NULL, family=gaussian(), corstruct="independence", 
    control=yags.control(), weights=NULL, betainit=NULL, alphainit=NULL, 
    data=list(), subset=NULL, allcrit=FALSE, lhetfam= NULL, qhetfam = NULL,
    icritalp = NULL, critar1tag = "ar1")
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{ GLM formula }
  \item{id}{ numeric vector of cluster discriminators; data assumed
sorted by this vector }
  \item{cor.met}{ correlation metameter: vector of observation
times for longitudinal data, matrix of coordinates
for other designs.  If unstructured working correlation,
this must have nonnegative integer values with zero origin. }
  \item{family}{ GLM family -- can use
default (gaussian -- linear link, constant variance),
poisson() (log link, variance = mean),
Gamma() (reciprocal link, variance = mu*mu),
binomial (logit link, binomial variance), 
quasi(link=log, var=\eqn{mu^2}), quasi(link = "identity", var=mu), 
or quasi(link=:identity:, var=\eqn{mu^2})}
  \item{corstruct}{ string describing working correlation model.
For Wang and Carey 2004 JASA article structure, use "UQ.fom"

Other options are "independence", "exchangeable", "UJ.fom", 
        "ar1", "unstructured".
 }
  \item{control}{ list of control parameters, see yags.control() }
  \item{weights}{ vector of weights }
  \item{betainit}{ initial value of regression parameters }
  \item{alphainit}{ initial value of working correlation model parameters }
  \item{data}{ data source for model fit }
  \item{subset}{ expression selecting a subset for fitting }
  \item{allcrit}{ logical -- if TRUE, compute all criteria for candidate models }
  \item{lhetfam}{ the glm family object defining the candidate model expressing
 linear heteroskedasticity }
  \item{qhetfam}{ the glm family object defining the candidate model expressing
 quadratic heteroskedasticity }
  \item{icritalp}{ numeric initial value of alpha parameter for candidate model fits
for criteria computation }
  \item{critar1tag}{ string defining the corstruct to be used for candidate model fits }
}

\value{
see \code{\link{yags.object}}
}
\references{ Liang KY, Zeger SL.  Biometrika 1986 }
\author{ VJ Carey, stvjc@gauss.med.harvard.edu }

\seealso{ \code{\link{glm}} }

\examples{
#
# trivial demo
#
library(methods)
data(stackloss)
Y1 <- yags(stack.loss~Air.Flow,id=1:21, data=stackloss)
#
# spruce examples for gaussian family and various corstructs
#
data(Spruce)
SPRind <- yags(y~Time+treated,id=Subject,data=Spruce)
SPRexc <- yags(y~Time+treated,id=Subject,data=Spruce,
		corstr="exchangeable", alphainit=0.)
# unstructured:
SPRuns <- yags(y~Time+treated,id=Subject,data=Spruce,
		corstr="unstructured", alphainit=rep(.1,45),
                cor.met=as.double(rep(0:9,79)))
# U_J of Wang and Carey JASA 2004
SPRUJ <- yags(y~Time+treated,id=Subject,data=Spruce,
		corstr="UJ.fom", alphainit=.1,
                cor.met=as.double(rep(0:9,79)))
#
# some criteria of adequacy over ranges of models
SPRUJ <- yags(y~Time+treated,id=Subject,data=Spruce,
                corstr="ar1", alphainit=.1,
                cor.met=as.double(rep(0:9,79)), allcrit=TRUE,
                  lhetfam=quasi(variance="mu"), qhetfam=quasi(variance="mu^2"))
sort(SPRUJ@m2LG)
sort(SPRUJ@del1)
}

\keyword{ models }
\keyword{ multivariate }

