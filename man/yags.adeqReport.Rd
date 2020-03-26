\name{yags.adeqReport}
\alias{yags.adeqReport}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ report on a yags fit }
\description{
reporting function that focuses on measures of GEE solution adequacy;
a component to be used by print.yagsResult; not intended to
be used directly
}
\usage{
yags.adeqReport(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class yagsResult (but not checked to be such) }
  \item{\dots}{ not currently used }
}
\details{
This is an interim solution to a problem of providing a
good class design for model reporting utilities.
}
\value{
a string including adequacy statistics including
the Rotnitzky-Jewell criterion and Wei Pan's AIC analogue.
}
%\references{  }
\author{ VJ Carey, stvjc@gauss.med.harvard.edu }
%\note{ }

%\seealso{ }

%\examples{ }
\keyword{ models }% at least one, from doc/KEYWORDS

