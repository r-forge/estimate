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
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\examples{
isempty(1:3)  # FALSE
isempty(array(NA, c(2, 0, 2)))	# TRUE
}
\keyword{logic}
\keyword{attribute}

