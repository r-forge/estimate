% $Id$

\name{ceil}
\alias{ceil}
\title{MATLAB ceil function}
\description{
  Rounds to the nearest integer.
}
\usage{
ceil(x)
}
\arguments{
  \item{x}{numeric to be rounded}
}
\details{
  Simply invokes \code{ceiling} for those more used to C library API name.
}
\value{
  Returns numeric vector containing smallest integers not less than the
  corresponding elements of argument \code{x}.
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\seealso{
  \code{\link{fix}},
  \code{\link[base]{Round}}
}
\examples{
ceil(c(0.9, 1.3, 2.4))
}
\keyword{arith}

