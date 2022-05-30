% $Id$

\name{fix}
\alias{fix}
\title{MATLAB fix function}
\description{
  Rounds toward zero.
}
\usage{
fix(A)
}
\arguments{
  \item{A}{numeric to be rounded}
}
\details{
  Simply invokes \code{\link[base:Round]{trunc}}.
}
\value{
  Returns vector containing integers by truncating the
  corresponding values of argument \code{A} toward zero.
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\seealso{ 
  \code{\link{ceil}},
  \code{\link[base]{Round}}
}   
\examples{
fix(c(1.3, 2.5, 3.7))
}
\keyword{arith}

