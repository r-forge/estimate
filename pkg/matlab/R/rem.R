###
### $Id$
###
### Remainder after division.
###


##-----------------------------------------------------------------------------
rem <- function(x, y) {
    ans <- matlab::mod(x, y)
    if (!((x > 0 && y > 0) ||
          (x < 0 && y < 0))) {
        ans <- ans - y
    }

    ans
}

