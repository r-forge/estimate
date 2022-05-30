###
### $Id$
###
### Remainder after division.
###


##-----------------------------------------------------------------------------
rem <- function(x, y) {
    ## Logically, it's this code but R doesn't work like MATLAB
    #ans <- matlab::mod(x, y)
    #if (!((x > 0 && y > 0) ||
    #      (x < 0 && y < 0))) {
    #    ans <- ans - y
    #}

    flip <- ifelse(y < 0, yes = TRUE, no = FALSE)
    if (flip) {
        x <- -x
        y <- -y
    }
    ans <- abs(x) %% y
    ans <- ifelse(x < 0, yes = -ans, no = ans)
    if (flip) {
        ans <- -ans
    }

    ans
}

