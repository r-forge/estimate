###
### $Id$
###
### Provides number of elements.
###


##-----------------------------------------------------------------------------
numel <- function(A, varargin) {
    if (!missing(varargin)) {
        stop("not implemented")         # need example
    }

    prod(matlab::size(A))
}

