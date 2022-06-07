###
### Package spell checking (R >= 3.0.0)
###
### Performed by hand.
### Assumes GNU Aspell and appropriate language dictionary is installed.
###
### $Id$
###

##
## Basic setup
## 
## R> library(utils)
## R> source("defaults.R")
## R> pkgname <- "rwt"           # replace value with name of package
##
## Process installed package
## R> pkg.path <- system.file(package=pkgname)
##        -- or --
## R> pkg.path <- file.path(..., pkgname)  # Modify as needed
##

##
## Spellcheck the meat of the package as follows:
## 
## R> utils:::aspell_package_description(pkg.path)  # not exported
## R> aspell_package_Rd_files(pkg.path)
## R> aspell_package_vignettes(pkg.path)
## R> aspell_package_R_files(pkg.path)
## R> aspell_package_C_files(pkg.path)  # template file 'po/PACKAGE.pot' only
##
## To get just the list of questioned words on manpages...
##
## R> asa <- aspell_package_Rd_files(pkg.path)
## R> unique(sort(asa$Original))
##

