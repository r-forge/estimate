% $Id$

\name{multiline.plot.colors}
\alias{multiline.plot.colors}
\title{MATLAB multiline plot colors}
\description{
  Creates a vector of colors equivalent to MATLAB's default colors to use
  for multiline plots.
}
\usage{
multiline.plot.colors()
}
\details{
  This is equivalent to the MATLAB command \preformatted{

    get(gca, 'ColorOrder')

  }
}
\value{
  Returns vector of color names. This can be used either to create
  a user-defined color palette for subsequent graphics, a \code{col=}
  specification in graphics functions, or in \code{par}.
}
\note{
  Method should be considered experimental and will most likely be removed
  and replaced with similar functionality in the near future.
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\seealso{
  \code{\link[grDevices]{palette}},
  \code{\link[graphics]{par}},
  \code{\link[grDevices]{rgb}}
}
\examples{
require(graphics)
x <- matrix(1:16, nrow=2, byrow=TRUE)
matplot(x, type="l", col=multiline.plot.colors())
}
\keyword{color}

