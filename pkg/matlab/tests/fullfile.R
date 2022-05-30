###
### $Id$
###


##-----------------------------------------------------------------------------
test.fullfile <- function(input, expected) {
    output <- do.call(getFromNamespace("fullfile", "matlab"), input)
    identical(output, expected)
}

fullfile.expected <- file.path(path.expand("~"), "somedir", "foo.txt")

test.fullfile(list(dir    = path.expand("~"),
                   subdir = "somedir",
                   file   = "foo.txt"), fullfile.expected)

