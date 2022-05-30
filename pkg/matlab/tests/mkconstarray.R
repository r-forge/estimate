###
### $Id$
###


##-----------------------------------------------------------------------------
test.mkconstarray <- function(input, expected) {
    output <- do.call(getFromNamespace("mkconstarray", "matlab"), input)
    identical(output, expected)
}

test.mkconstarray(list(class.type = "double", value = pi, size = 4), rep(pi, 4))

