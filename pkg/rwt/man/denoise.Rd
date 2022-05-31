% $Id$

\name{denoise}
\alias{denoise}
\alias{denoise.dwt}
\alias{denoise.udwt}
\alias{DWT.TRANSFORM.TYPE}
\alias{UDWT.TRANSFORM.TYPE}
\alias{MAD.VARIANCE.ESTIMATOR}
\alias{STD.VARIANCE.ESTIMATOR}
\alias{SOFT.THRESHOLD.TYPE}
\alias{HARD.THRESHOLD.TYPE}
\alias{MAX.DECOMPOSITION}
\alias{CALC.THRESHOLD.TO.USE}
\alias{DEFAULT.DWT.THRESHOLD.MULTIPLIER}
\alias{DEFAULT.UDWT.THRESHOLD.MULTIPLIER}
\alias{default.dwt.option}
\alias{default.udwt.option}
\title{Wavelet-based Denoising}
\description{
  Denoise the signal \code{x} using the 2-band wavelet system described
  by the filter \code{h} using either the traditional discrete wavelet
  transform (DWT) or the linear shift invariant discrete wavelet
  transform (also known as the undecimated DWT (UDWT)).
}
\usage{
denoise(x, h, type, option)
denoise.dwt(x, h, option = default.dwt.option)
denoise.udwt(x, h, option = default.udwt.option)
}
\arguments{
  \item{x}{1D or 2D signal to be denoised}
  \item{h}{numeric scalar specifying scaling filter to be applied}
  \item{type}{type of transform. Valid values are:
  \tabular{l}{
    \code{DWT.TRANSFORM.TYPE}\cr
    \code{UDWT.TRANSFORM.TYPE}
  }}
  \item{option}{list containing desired transformation settings}
}
\details{
  The transformation settings in the \code{option} list are:

  \describe{
    \item{threshold.low.pass.part:}{
      logical scalar. If \code{TRUE}, threshold the low-pass component.}
    \item{threshold.multiplier:}{\code{thld = c*MAD(noise_estimate)} }
    \item{variance.estimator:}{Valid values are:
      \tabular{ll}{
        MAD.VARIANCE.ESTIMATOR \tab Mean absolute deviation\cr
        STD.VARIANCE.ESTIMATOR \tab Classical numerical std estimate
      }
    }
    \item{threshold.type:}{Valid values are:
      \tabular{ll}{
        SOFT.THRESHOLD.TYPE \tab Soft thresholding\cr
        HARD.THRESHOLD.TYPE \tab Hard thresholding
      }
    }
    \item{num.decompression.levels:}{number of levels in wavelet
      decomposition. Setting this to \code{MAX.DECOMPOSITION} will allow
      maximal decomposition.
    }
    \item{threshold:}{actual threshold to use. Setting this to anything but
      \code{CALC.THRESHOLD.TO.USE} will disable the \code{variance.estimator}
      setting.
    }
  }
}
\value{
  Returns a list with components:
    \item{xd}{estimate of noise free signal}
    \item{xn}{estimated noise signal (\code{x}-\code{xd})}
    \item{option}{list of actual parameters used. It is configured
      the same way as the input option list with an additional
      element - \code{option[[7]] = type}.
    }
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\note{
Both \code{denoise.dwt} and \code{denoise.udwt} are convenience routines
that call the \code{denoise} routine with appropriate default arguments.
}
\examples{
sig <- makesig(SIGNAL.DOPPLER)
h <- daubcqf(6)
ret.dwt <- denoise.dwt(sig$x, h$h.0)
}
\keyword{interface}

