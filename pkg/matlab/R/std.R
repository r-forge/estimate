###
### $Id$
###
### Standard deviation.
###


##-----------------------------------------------------------------------------
std <- function(x, flag = 0) {
    if (flag != 0) {
        stop("biased standard deviation not implemented")
    }

    return(sd(x))
}

