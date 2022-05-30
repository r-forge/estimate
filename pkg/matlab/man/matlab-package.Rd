% $Id$

\name{matlab-package}
\alias{matlab-package}
\docType{package}
\title{MATLAB Emulation Functions}
\description{
  Wrapper functions and variables used to replicate MATLAB
  function calls as best possible to simplify porting.
}
\details{
  \tabular{ll}{
    Package: \tab matlab\cr
    Type:    \tab Package\cr
    Version: \tab \Sexpr[results=text]{packageDescription("matlab")$Version}\cr
    Date:    \tab \Sexpr[results=text]{packageDescription("matlab")$Date}\cr
    License: \tab \Sexpr[results=text]{packageDescription("matlab")$License}\cr
  }
  They are no more complete than absolutely necessary and are quite
  possibly broken for fringe cases.

  For a complete list of functions, use \code{library(help="matlab")}.\cr
  For a high-level summary of the changes for each revision, use
    \code{file.show(system.file("NEWS", package="matlab"))}.
}
\note{
  In certain cases, these may not correspond exactly with MATLAB API
  as sometimes it just wasn't possible.
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\keyword{package}

