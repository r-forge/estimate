% $Id$

\name{mrdwt}
\alias{mrdwt}
\title{Redundant Discrete Wavelet Transform}
\description{
  Computes the redundant discrete wavelet transform \code{y} for
  input signal \code{x} using the scaling filter \code{h}. Redundant
  means here that the sub-sampling after each stage is omitted.
}
\usage{
mrdwt(x, h, L)
}
\arguments{
  \item{x}{Finite 1D or 2D signal (implicitly periodized)}
  \item{h}{Scaling filter to be applied}
  \item{L}{Number of levels in wavelet decomposition.
In the case of a 1D signal, length(x) must be divisible by \eqn{2^L};
in the case of a 2D signal, the row and the column dimension must be
divisible by \eqn{2^L}. If no argument is specified, a full DWT is returned
for maximal possible L.
  }
}
\value{
  Returns a list with components:
    \item{yl}{Lowpass component}
    \item{yh}{Highpass components}
    \item{L}{Number of levels in wavelet decomposition}
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\examples{
sig <- makesig(SIGNAL.LEOPOLD, 8)
h <- daubcqf(4)
L <- 1
ret.mrdwt <- mrdwt(sig$x, h$h.0, L)
}
\keyword{interface}

