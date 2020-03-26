\name{yags.object}
\alias{yags.object}
\non_function{}
\title{ object returned by yags() }
\usage{none}
\description{
yags returns a list of the following structure:
}
\format{
\describe{
\item{coefficients}{
a numeric p-vector: regression coefficients
}
\item{naive.parmvar}{
the pxp model-based parameter covariance matrix
}
\item{robust.parmvar}{
the pxp model-robust parameter covariance matrix
}
\item{alpha}{
the q-vector of working correlation parameters
}
\item{phi}{
the scale parameter estimate
}
\item{linear.predictors}{
a numeric vector
}
\item{fitted.values}{
a numeric vector
}
\item{residuals}{
a numeric vector
}
\item{iter}{
a numeric vector
}
\item{family}{
a function
}
\item{rank}{
a numeric vector
}
\item{errorcode}{
an indicator 
}
\item{sealp}{
a numeric vector of standard errors
}
\item{qls}{
the QLS (Chaganty/Shults) criterion value
}
\item{pan.aic}{
the AIC for GEE of Wei Pan (Bcs 2000)
}
\item{sealp.OK}{
a logical vector
}
}
}
%\details{
%}
\keyword{ models }
%\keyword{GEE}

