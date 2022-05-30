% $Id$

\name{jet.colors}
\alias{jet.colors}
\title{MATLAB jet function}
\description{
  Creates a vector of \code{n} colors beginning with dark blue, ranging through
  shades of blue, cyan, green, yellow and red, and ending with dark red.
}
\usage{
jet.colors(n)
}
\arguments{
  \item{n}{numeric scalar specifying number of colors to be in the palette}
}
\value{
  Returns vector of \code{n} color names. This can be used either to create
  a user-defined color palette for subsequent graphics, a \code{col=}
  specification in graphics functions, or in \code{par}.
}
\seealso{
  \code{\link[grDevices]{palette}},
  \code{\link[graphics]{par}},
  \code{\link[grDevices]{rgb}}
}
\examples{
require(graphics)
x <- 1:16
pie(x, col=jet.colors(length(x)))
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\keyword{color}

