% $Id$

\name{strcmp}
\alias{strcmp}
\title{MATLAB strcmp function}
\description{
  Compare strings.
}
\usage{
strcmp(S, T)
}
\arguments{
  \item{S, T}{character vectors to evaluate}
}
\details{
  Comparisons are case-sensitive and any leading and trailing blanks in
  either of the strings are explicitly included in the comparison.
}
\value{
  Returns \code{TRUE} if \code{S} is identical to \code{T};
  otherwise, \code{FALSE}.
}
\note{
  Value returned is the opposite of the C language convention.
}
\examples{
strcmp("foo", "bar")  # FALSE
strcmp(c("yes", "no"), c("yes", "no"))	# TRUE
}
\author{
  P. Roebuck \email{proebuck@mdanderson.org}
}
\keyword{character}

