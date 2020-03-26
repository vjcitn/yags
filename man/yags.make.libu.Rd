\name{yags.make.libu}
\alias{yags.make.libu}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ create libu.so in current working directory }
\description{
THIS IS ONLY A POTENTIAL FUNCTIONALITY -- IF YOU NEED TO USE IT CONTACT
VJ CAREY --

creates a trivial C++ file to allow dynamic loading to
succeed and as a template for C++ extensions to the
working correlation model repertory of yags


}
\usage{
# NOT FUNCTIONING UNLESS EXPERIMENTAL_USERCOR defined in C++
# IGNORE UNLESS YOU REALLY KNOW THE SOURCE CODE
yags.make.libu(lib.loc=.lib.loc, cxx="g++ -c")
}
\arguments{
  \item{lib.loc}{ path to library where yags is installed, to
allow access to include files }
  \item{cxx}{ string encoding C++ compiler invocation }
}
\details{
A small C++ program is written to the file user\_wcor.cc,
compiled according to the value of \code{cxx}, and
transformed to \code{libu.so} by \code{R SHLIB}.
The program can be modified and relinked to specify
calculations of working correlation models selectable through
corstruc="user" in yags.
}
\value{
Nothing returned, all side-effects.  You must typically
exit R, reset LD\_LIBRARY\_PATH, and restart in order
for yags to work.
}
%\references{  }
\author{ VJ Carey, stvjc@channing.harvard.edu }
\note{ This approach is dependent on a fairly obsolete approach (MC++)
to matrix computations in C++.  The next version is likely
to be based on the .Call interface, and high-performance user-written extensions
will follow R rather than MC++ protocols. }


%\seealso{  }

\examples{
\dontrun{
yags.make.libu()
}
}
\keyword{ models }% at least one, from doc/KEYWORDS
%\keyword{ C++ }% at least one, from doc/KEYWORDS

