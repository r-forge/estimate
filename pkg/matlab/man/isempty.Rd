% $Id$

\name{isempty}
\alias{isempty}
\title{MATLAB isempty function}
\description{
  Determine if object is empty.
}
\usage{
isempty(A)
}
\arguments{
  \item{A}{object to evaluate}
}
\details{
  An empty object has at least one dimension of size zero.
}
\value{
  Returns \code{TRUE} if \code{x} is an empty object; otherwise, \code{FALSE}.
}
\examples{
isempty(1:3)  # FALSE
isempty(array(NA, c(2, 0, 2)))	# TRUE
}
\author{
  P. Roebuck \email{proebuck@mdanderson.org}
}
\keyword{logic}
\keyword{attribute}

