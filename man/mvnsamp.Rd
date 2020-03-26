\name{mvnsamp}
\alias{mvnsamp}
\title{ create a balanced sample of MVN outcomes for input to gee-like programs }
\description{
create a balanced sample of MVN outcomes for input to gee-like programs 
}
\usage{
mvnsamp(n=100, m=rep(0, 4), v=diag(4))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n}{ number of clusters }
  \item{m}{ number of entries/cluster }
  \item{v}{ covariance }
}
%\details{
%}
\value{
a data.frame with columns y, id and tim
}

\examples{
mydf <- mvnsamp()
with(mydf, yags(y~tim,id=id,corstr="exchangeable", alphainit=.1))
}
\keyword{multivariate}

