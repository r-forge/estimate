% $Id$

\name{find}
\alias{find}
\title{MATLAB find function}
\description{
  Finds indices of elements.
}
\usage{
find(x)
}
\arguments{
  \item{x}{expression to evaluate}
}
\details{
  If expression is not logical, finds indices of nonzero elements of
  argument \code{x}.
}
\value{
  Returns indices of corresponding elements matching the expression \code{x}.
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\examples{
find(-3:3 >= 0)
find(c(0, 1, 0, 2, 3))
}
\keyword{logic}
\keyword{attribute}

