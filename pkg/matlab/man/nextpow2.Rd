% $Id$

\name{nextpow2}
\alias{nextpow2}
\title{MATLAB nextpow2 function}
\description{
  Smallest power of 2 greater than or equal to its argument.
}
\usage{
nextpow2(x)
}
\arguments{
  \item{x}{numeric or complex value(s).}
}
\details{
  Computes the smallest power of two that is greater than or equal to the
  absolute value of \code{x}. (That is, \code{p} that satisfies
  \eqn{2^p \ge abs(x)}).
  For negative or complex values, the absolute value will be taken.
}
\value{
  Returns numeric result containing integer \code{p} as described above.
  Nonscalar input returns an element-by-element result (of same size/dimensions
  as its input).
}
\author{
  H. Borchers \email{hwborchers@googlemail.com},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\seealso{
  \code{\link{pow2}}
}
\examples{
nextpow2(10)                   # 4
nextpow2(1:10)                 # 0 1 2 2 3 3 3 3 4 4
nextpow2(-2^10)                # 10
nextpow2(.Machine$double.eps)  # -52
nextpow2(c(0.5, 0.25, 0.125))  # -1 -2 -3
}
\keyword{arith}

