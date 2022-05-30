% $Id$

\name{pathsep}
\alias{pathsep}
\title{MATLAB pathsep function}
\description{
  Returns the character that separates directory names in a list such
  as the \env{PATH} environment variable.
}
\usage{
pathsep
}
\details{
  Variable that contains the value of \code{.Platform$path.sep}.
}
\value{
  Returns character representing this platform's path separator.
}
\note{
  Implemented as an R variable rather than a function such that it more
  closely resembles normal MATLAB usage. 
}
\seealso{
  \code{\link{filesep}}
}
\author{
  P. Roebuck \email{proebuck@mdanderson.org}
}
\keyword{file}

