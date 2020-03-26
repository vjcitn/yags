\name{rmvlnorm}
\alias{rmvlnorm}
\title{ simulated correlated lognormal deviates, according to
appendix of Wang and Carey 2004 JASA }
\description{ simulated correlated lognormal deviates, according to
appendix of Wang and Carey 2004 JASA }
\usage{
rmvlnorm(n = 100, lmu = c(2, 3, 4), lcov = diag(c(2, 3, 4),3,3))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ number of p-vectors to be returned }
  \item{lmu}{ p-vector of response mean }
  \item{lcov}{ pxp matrix of response covariance }
}
%\details{
%}
\value{
an n x p matrix of marginally lognormal responses with
specified mean and covariance
}
\references{ Wang YG, Carey VJ, Unbiased estimating equations ...,
JASA Sept 2004, v99, no 467, 845-853 }
\author{ VJ Carey }
%\note{ ~~further notes~~ 

% ~Make other sections like Warning with \section{Warning }{....} ~
%}
%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
mm = matrix(.75, 3, 3)
diag(mm) = 1
var( rmvlnorm( 1000, c(2,3,4), mm ) )
}
\keyword{ models }
