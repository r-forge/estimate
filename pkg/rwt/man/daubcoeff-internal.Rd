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
\details{
  Supported values of \code{N} are [4|6|8|10|12|14|16|18|20].
}
\value{
  Returns a vector containing coefficients.
}
\references{
  \url{https://cas.mines-paristech.fr/~chaplais/Wavetour_presentation/Tables/Daubechies/Daubechies_coefficients.csv}
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\keyword{internal}

