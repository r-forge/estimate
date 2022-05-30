% $Id$

\name{std}
\alias{std}
\title{MATLAB std function}
\description{
  Computes the standard deviation of the values of \code{x}.
}
\usage{
std(x, flag=0)
}
\arguments{
  \item{x}{numeric vector or matrix}
  \item{flag}{numeric scalar. If \code{0}, selects unbiased algorithm.
  If \code{1}, selects biased algorithm (currently unsupported).}
}
\details{
  Simply invokes \code{\link[stats]{sd}}.
}
\value{
  Return value depends on argument \code{x}. If vector, returns the
  standard deviation. If matrix, returns vector containing the
  standard deviation of each column.
}
\examples{
std(1:2) ^ 2
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\keyword{univar}

