###
### $Id$
###
### Provides the number of dimensions.
###


##-----------------------------------------------------------------------------
ndims <- function(A) {
    return(length(matlab::size(A)))
}

