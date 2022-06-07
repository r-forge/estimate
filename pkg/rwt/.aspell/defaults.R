###
### Package spell checking configuration file (R >= 3.0.0)
###
### Performed as part of R CMD check.
### Assumes GNU Aspell and appropriate language dictionary is installed.
###
### $Id$
###

lang <- "en"
pkgname <- "rwt"           # replace value with name of package


##
## Settings used by "aspell-utils" package
##

## Dictionaries to use
aspell_dictionaries_R <- (function() {
    ## Recreate 'utils:::aspell_dictionaries_R' (aka "en_stats")
    R_dict_path <- file.path(R.home("share"), "dictionaries")
    tools::file_path_sans_ext(list.files(R_dict_path))
})()
aspell_dictionaries_pkg <- sprintf("%s_%s", lang, pkgname)
names(aspell_dictionaries_pkg) <- pkgname

## Generic settings for all files
genericParams <- list(encoding = "UTF-8",
                      language = lang,
                      dictionaries = c(aspell_dictionaries_R,
                                       aspell_dictionaries_pkg))

## Settings used by "aspell-utils" methods for each type of package file
R_files <- C_files <- vignettes <- description <- genericParams
Rd_files <- append(genericParams,
                   list(drop = c("\\author",
                                 "\\references",
                                 "\\source")))
rm(lang, pkgname, genericParams)

