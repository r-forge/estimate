###
### $Id$
###
### Find indices of nonzero elements.
###


##-----------------------------------------------------------------------------
find <- function(x) {
    expr <- if (is.logical(x)) {
                x
            } else {
                x != 0
            }
    which(expr)
}

