% $Id$

\name{mirdwt}
\alias{mirdwt}
\title{Inverse Redundant Discrete Wavelet Transform}
\description{
  Computes the inverse redundant discrete wavelet transform \code{x} for
  input signal \code{y} using the scaling filter \code{h}. (Redundant means
  here that the sub-sampling after each stage of the forward transform has
  been omitted.)
}
\usage{
mirdwt(yl, yh, h, L)
}
\arguments{
  \item{yl}{Lowpass component}
  \item{yh}{Highpass components}
  \item{h}{Scaling filter to be applied}
  \item{L}{Number of levels in wavelet decomposition.
In the case of a 1D signal, length(yl) must be divisible by \eqn{2^L};
in the case of a 2D signal, the row and the column dimension must be
divisible by \eqn{2^L}. If no argument is specified, a full DWT is returned
for maximal possible L.
  }
}
\value{
  Returns a list with components:
    \item{x}{Finite length 1D or 2D signal}
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
ret.mirdwt <- mirdwt(ret.mrdwt$yl, ret.mrdwt$yh, h$h.0, ret.mrdwt$L)
}
\keyword{interface}

