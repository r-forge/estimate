% $Id$

\name{padarray}
\alias{padarray}
\alias{padarray,array,numeric,missing,missing-method}
\alias{padarray,array,numeric,ANY,character-method}
\alias{padarray,array,numeric,character,character-method}
\alias{padarray,vector,numeric,ANY,ANY-method}
\title{MATLAB padarray function}
\description{
  Pad array.
}
\usage{
padarray(A, padsize, padval=0, direction=c("both", "pre", "post"))
}
\arguments{
  \item{A}{vector, matrix, or array to be padded}
  \item{padsize}{integer vector specifying both amount of padding and
    the dimension along which to add it}
  \item{padval}{scalar value specifying pad value, which defaults to 0. \cr
    Instead, it may specify the method used to determine pad values. \cr
    Valid values for the method are:
    \tabular{ll}{
      \code{"circular"}  \tab pad with circular repetition of elements within
                              the dimension \cr
      \code{"replicate"} \tab pad by repeating border elements of array \cr
      \code{"symmetric"} \tab pad array with mirror reflections of itself \cr
    }
  }
  \item{direction}{character string specifying direction to apply padding. \cr
    Valid values are:
    \tabular{ll}{
      \code{"both"} \tab pad before first element and after last array
                         element along each dimension \cr
      \code{"pre"}  \tab pad after last array element along each dimension \cr
      \code{"post"} \tab pad before first array element along each dimension \cr
    }
  }
}
\details{
  This is an S4 generic function.
}
\value{
  Return value is the same type as argument \code{A} with requested padding.
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\examples{
padarray(1:4, c(0, 2))				# 0 0 [1 2 3 4] 0 0
padarray(1:4, c(0, 2), -1)			# -1 -1 [1 2 3 4] -1 -1
padarray(1:4, c(0, 2), -1, "post")		# [1 2 3 4] -1 -1
padarray(1:4, c(0, 3), "symmetric", "pre")	# 3 2 1 [1 2 3 4]
padarray(letters[1:5], c(0, 3), "replicate")	# a a a [a b c d e] e e e
padarray(letters[1:5], c(0, 3), "circular", "post")	# [a b c d e] a b c
}
\keyword{array}

