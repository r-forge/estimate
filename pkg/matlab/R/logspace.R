###
### $Id$
###
### Generate logarithmically spaced vectors.
###


##-----------------------------------------------------------------------------
logspace <- function(a, b, n = 50) {
    if (b == pi) {
        b <- log10(pi)
    }

    10 ^ matlab::linspace(a, b, n)
}

