% $Id$

\name{rosser}
\alias{rosser}
\title{MATLAB rosser function}
\description{
  Create the Rosser matrix, a classic symmetric eigenvalue test problem.
}
\usage{
rosser()
}
\details{
  The returned matrix has the following features:
  \itemize{
    \item{a double eigenvalue}
    \item{three nearly equal eigenvalues}
    \item{dominant eigenvalues of opposite sign}
    \item{a zero eigenvalue}
    \item{a small, nonzero eigenvalue}
  }
}
\value{
  Returns an \code{8}-by-\code{8} matrix with integer elements.
}
\examples{
rosser()
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\keyword{array}

