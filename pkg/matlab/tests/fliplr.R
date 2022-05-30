###
### $Id$
###


##-----------------------------------------------------------------------------
test.fliplr <- function(input, expected) {
    output <- do.call(getFromNamespace("fliplr", "matlab"), input)
    identical(output, expected)
}

X.mat <- matrix(1:6, 2, 3, byrow = TRUE)
fliplr.expected.X.mat <- matrix(c(3:1, 6:4), 2, 3, byrow = TRUE)

test.fliplr(list(X.mat), fliplr.expected.X.mat)

X.vec <- seq(1, 9, by = 2)
fliplr.expected.X.vec <- seq(9, 1, by = -2)

test.fliplr(list(X.vec), fliplr.expected.X.vec)

