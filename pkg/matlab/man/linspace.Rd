% $Id$

\name{linspace}
\alias{linspace}
\title{MATLAB linspace function}
\description{
  Generate linearly spaced vectors.
}
\usage{
linspace(a, b, n=100)
}
\arguments{
  \item{a}{numeric scalar specifying starting point}
  \item{b}{numeric scalar specifying ending point}
  \item{n}{numeric scalar specifying number of points to be generated}
}
\details{
  Similar to colon operator but gives direct control over the number of points.
  Note also that although MATLAB doesn't specifically document this, the number
  of points generated is actually \code{floor(n)}.
}
\value{
  Returns vector containing containing \code{n} points linearly spaced
  between \code{a} and \code{b} inclusive. If \eqn{\code{n} < 2}, the
  result will be the ending point \code{b}.
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\seealso{
  \code{\link{logspace}}
}
\examples{
linspace(1, 10, 4)
}
\keyword{manip}

