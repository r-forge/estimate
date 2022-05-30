###
### $Id$
###
### Provides the number of dimensions.
###


##-----------------------------------------------------------------------------
ndims <- function(A) {
    length(matlab::size(A))
}

