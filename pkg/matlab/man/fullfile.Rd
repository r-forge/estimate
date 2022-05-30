% $Id$

\name{fullfile}
\alias{fullfile}
\title{MATLAB fullfile function}
\description{
  Contructs path to a file from components in platform-independent manner
}
\usage{
fullfile(\dots)
}
\arguments{
  \item{\dots}{character strings representing path components}
}
\details{
  Builds a full filename from the directories and filename specified. This
  is conceptually equivalent to \preformatted{

    paste(dir1, dir2, dir3, filename, sep=filesep)

  }
  with care taken to handle cases when directories begin or end with a
  separator.
}
\value{
  Returns character vector of arguments concatenated term-by-term and
  separated by file separator if all arguments have a positive length;
  otherwise, an empty character vector.
}
\seealso{
  \code{\link{fileparts}},
  \code{\link{filesep}}
}
\examples{
fullfile("", "etc", "profile")	# /etc/profile
}
\author{
  P. Roebuck \email{proebuck1701@gmail.com}
}
\keyword{file}

