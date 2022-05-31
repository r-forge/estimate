% $Id$

\name{threshold}
\alias{hardTh}
\alias{softTh}
\title{Threshold Input Signal}
\description{
  Thresholds the input signal \code{y} with the threshold value \code{thld}.
}
\usage{
hardTh(y, thld)
softTh(y, thld)
}
\arguments{
  \item{y}{1D or 2D signal to be thresholded}
  \item{thld}{numeric threshold value to be applied}
}
\value{
  \item{x}{Thresholded output}
}
\references{
  D.L. Donoho\cr
  \emph{De-noising via Soft-Thresholding}\cr
  Tech. Rept. Statistics, Stanford (1992)
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\examples{
sig <- makesig(SIGNAL.WERNER.SORROWS, 8)
thld <- 1
x <- rwt:::hardTh(sig$x, thld)
}
\keyword{interface}

