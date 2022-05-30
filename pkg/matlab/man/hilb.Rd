% $Id$

\name{hilb}
\alias{hilb}
\title{MATLAB hilb function}
\description{
  Create a Hilbert matrix.
}
\usage{
hilb(n)
}
\arguments{
  \item{n}{numeric scalar specifying dimensions for the result}
}
\details{
  The Hilbert matrix is a notable example of a poorly conditioned matrix.
  Its elements are \preformatted{
    H[i, j] = 1 / (i + j - 1)}.
}
\value{
  Returns an \code{n}-by-\code{n} matrix constructed as described above.
}
\examples{
hilb(3)
}
\author{
  H. Borchers \email{hwborchers@googlemail.com},
  P. Roebuck \email{proebuck@mdanderson.org}
}
\keyword{array}

