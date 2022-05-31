% $Id$

\name{daubcoeff}
\alias{.getDaubechiesCoefficients}
\title{Daubechies Coefficients}
\description{
  Returns the Daubechies coefficients for given number of coefficients.
}
\usage{
.getDaubechiesCoefficients(N)
}
\arguments{
  \item{N}{numeric scalar specifying number of coefficients desired}
}
\value{
  Returns a vector containing coefficients.
}
\details{
  Supported values of \code{N} are [4|6|8|10|12|14|16|18|20].
}
\references{
  \url{https://cas.ensmp.fr/~chaplais/Wavetour_presentation/Tables/Daubechies/}
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\keyword{internal}

