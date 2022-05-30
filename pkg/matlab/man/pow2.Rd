% $Id$

\name{pow2}
\alias{pow2}
\title{MATLAB pow2 function}
\description{
  Power with base 2.
}
\usage{
pow2(f, e)
}
\arguments{
  \item{f}{numeric vector of factors}
  \item{e}{numeric vector of exponents for base 2}
}
\details{
  Computes the expression \code{f * 2^e} for corresponding elements of \code{f}
  and \code{e}. If \code{e} is missing, it sets \code{e} to \code{f}
  and \code{f} to 1.
  Imaginary parts of complex values are ignored unless \code{e} is missing.
}
\value{
  Returns numeric vector constructed as described above.
}
\author{
  H. Borchers \email{hwborchers@googlemail.com},
  P. Roebuck \email{proebuck1701@gmail.com}
}
\seealso{
  \code{\link{nextpow2}}
}
\examples{
pow2(c(0, 1, 2, 3))               # 1 2 4 8
pow2(c(0, -1, 2, 3), c(0,1,-2,3)) # 0.0 -2.0  0.5 24.0
pow2(1i)                          # 0.7692389+0.6389613i

# For IEEE arithmetic...
pow2(1/2, 1)                      # 1
pow2(pi/4, 2)                     # pi
pow2(-3/4, 2)                     # -3
pow2(1/2, -51)                    # .Machine$double.eps
pow2(1/2, -1021)                  # .Machine$double.xmin
}
\keyword{arith}

