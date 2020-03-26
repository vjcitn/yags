\name{Spruce}
\alias{Spruce}
\non_function{}
\title{ Spruce data from oswald }
\usage{data(Spruce)}
\description{
The \code{Spruce} data frame has 790 rows and 5 columns.
Two strata: control, and treated with ozone. 
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{a numeric vector}
    \item{Time}{a numeric vector}
    \item{Subject}{a numeric vector}
    \item{treated}{a numeric vector}
    \item{trtfac}{a factor with levels control, treated}
  }
}
\details{
Compare to V\&R sitka data.
}
\source{
See Diggle, Liang and Zeger 1994, and the oswald distribution at 
statlib.
}
%\references{
%}
\examples{
data(Spruce)
## maybe str(Spruce) ; plot(Spruce) ...
}
\keyword{datasets}

