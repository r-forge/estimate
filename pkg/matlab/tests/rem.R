###
### $Id$
###


##-----------------------------------------------------------------------------
test.rem <- function(input, expected) {
    output <- do.call(getFromNamespace("rem", "matlab"), input)
    identical(output, expected)
}


## Remainder After Division of Scalar
rem.expected.div.scalar <- 3
test.rem(list(x = 23, y = 5), rem.expected.div.scalar)

## Remainder After Division of Vector
X.vec <- 1:5
rem.expected.div.vec <- c(1, 2, 0, 1, 2)
test.rem(list(x = X.vec, y = 3), rem.expected.div.vec)

## Remainder After Division of Vector by Zero
X.vec <- 1:5
rem.expected.div.vec.by0 <- rep(NaN, length(X.vec))
test.rem(list(x = X.vec, y = 0), rem.expected.div.vec.by0)

## Remainder After Division for Positive and Negative Values
## Note that nonzero results have the same sign as the dividend.
X.posneg.vec <- c(-4, -1, 7, 9)
rem.expected.div.posneg.vec <- c(-1, -1, 1, 0)
test.rem(list(x = X.posneg.vec, y = 3), rem.expected.div.posneg.vec)

## Remainder After Division for Floating-Point Values
X.theta <- c(0.0, 3.5, 5.9, 6.2, 9.0, 4 * pi)
b <- 2 * pi;
expected.div.fp.vec <- c(0, 3.5, 5.9, 6.2, 2.716815, 0)
test.rem(list(x = X.theta, y = b), expected.div.fp.vec)

## Remainder After Division of Matrix
X.mat <- matrix(1:9, nrow = 3, byrow = TRUE)
rem.expected.X.mat <- matrix(c(3:1, 6:4, 9:7), nrow = 3, byrow = TRUE)
rem.expected.X.mat.Y0 <- matrix(rep(NaN, length(X.mat)), nrow = nrow(X.mat))
rem.expected.X.mat.Y1 <- matrix(rep(1, length(X.mat)), nrow = nrow(X.mat))
rem.expected.X.mat.Y2 <- matrix(rep(c(1, 0), 5)[1:length(X.mat)],
                                nrow = nrow(X.mat))
rem.expected.X.mat.Y3 <- matrix(rep(c(1, 2, 0), nrow(X.mat)),
                                nrow = nrow(X.mat),
                                byrow = TRUE)

test.rem(list(x = X.mat, y = 0), rem.expected.X.mat.Y0)
test.rem(list(x = X.mat, y = 1), rem.expected.X.mat.Y1)
test.rem(list(x = X.mat, y = 2), rem.expected.X.mat.Y2)
test.rem(list(x = X.mat, y = 3), rem.expected.X.mat)


## rem & mod give same results with X, Y having same sign
test.rem(list(x = 5, y = 3), matlab::mod(5, 3))
test.rem(list(x = -5, y = -3), matlab::mod(-5, -3))

## alternate formula used when X, Y having different signs
test.rem(list(x = 5, y = -3), (matlab::mod(5, -3) - -3))
test.rem(list(x = -5, y = 3), (matlab::mod(-5, 3) - 3))

