% $Id$

\name{logspace}
\alias{logspace}
\title{MATLAB logspace function}
\description{
  Generate logarithmically spaced vectors.
}
\usage{
logspace(a, b, n=50)
}
\arguments{
  \item{a}{numeric scalar specifying exponent for starting point}
  \item{b}{numeric scalar specifying exponent for ending point}
  \item{n}{numeric scalar specifying number of points to be generated}
}
\details{
  Useful for creating frequency vectors, it is a logarithmic equivalent
  of \code{linspace}.
}
\value{
  Returns vector containing containing \code{n} points logarithmically
  spaced between decades \eqn{10^\code{a}} and \eqn{10^\code{b}}. For
  \eqn{n < 2}, \code{b} is returned.
}
\seealso{
  \code{\link{linspace}}
}
\examples{
logspace(1, pi, 36)
}
\author{
  P. Roebuck \email{proebuck@mdanderson.org}
}
\keyword{manip}

