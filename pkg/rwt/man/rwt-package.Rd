% $Id$

\name{rwt-package}
\alias{rwt-package}
\docType{package}
\title{
  \Sexpr[stage=build,results=hide]{pd <- utils::packageDescription("rwt")}
  \Sexpr[stage=build,results=text]{pd$Title}
}
\description{
  \Sexpr[stage=build,results=text]{pd$Description}
}
\details{
  \tabular{ll}{
    Package: \tab \Sexpr[stage=build,results=text]{pd$Package}\cr
    Type:    \tab package\cr
    Version: \tab \Sexpr[stage=build,results=text]{pd$Version}\cr
    Date:    \tab \Sexpr[stage=build,results=text]{pd$Date}\cr
    License: \tab \Sexpr[stage=build,results=text]{pd$License}\cr
  }

  For a complete list of functions, use \code{library(help="rwt")}.\cr
  For a high-level summary of the changes for each revision, use
    \code{file.show(system.file("NEWS", package="rwt"))}.
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\keyword{package}

