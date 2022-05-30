###
### $Id$
###
### Determine if object is empty.
###


##-----------------------------------------------------------------------------
isempty <- function(A) {
    any(matlab::size(A) == 0)
}

