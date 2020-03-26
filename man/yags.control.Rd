\name{yags.control}
\alias{yags.control}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ control iteration, etc. for yags() }
\description{
yags can involve a nested iteration procedure with
solutions to an auxiliary estimating equation Ua embedded
in the standard NR.  the iteration process is controlled
by parameters set in this function.
}
\usage{
yags.control(maxiter=15, tol=1e-04, verbose=FALSE, Ua.maxit=20, Ua.tol=0.001, 
      Ua.gridlo=0.1, Ua.gridhi=0.9, Ua.gridnpts=10, Ua.secantdel=0.01, 
      fixscale=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{maxiter}{ bound on number of NR iterations }
  \item{tol}{ convergence criterion for NR }
  \item{verbose}{ should extra info on iterations be emitted? }
  \item{Ua.maxit}{ max number iterations for auxiliary estimating equation solver }
  \item{Ua.tol}{ convergence criterion for aux. solver }
  \item{Ua.gridlo}{ aux solver starts with a grid search: lower starting point }
  \item{Ua.gridhi}{ upper bound to starting grid }
  \item{Ua.gridnpts}{ number of nodes in starting grid }
  \item{Ua.secantdel}{ delta for secant solution }
  \item{fixscale}{ is scale parameter fixed? }
}
%\details{
%}
%\value{
%}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }

%\seealso{  }

%\examples{
%}

\keyword{ models }

