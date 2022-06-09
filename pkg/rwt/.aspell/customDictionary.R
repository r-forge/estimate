###
### Package spell checking (R >= 3.0.0)
###
### Creation of custom dictionary.
### Assumes GNU Aspell and appropriate language dictionary is installed.
### Also assumes write permission to the package installation directory.
###
### $Id$
###

##
## Custom dictionaries
##

##
## Basic setup
##
## R> source('defaults.R')
## R> pkgname <- "rwt"           # replace value with name of package
##
## Process installed package
## R> pkg.path <- system.file(package=pkgname)
##        -- or --
## R> pkg.path <- file.path(..., pkgname)  # Modify as needed
##
## R> dict.file <- sprintf("%s.rds", aspell_dictionaries_pkg[1])
## R> dict.path <- file.path(pkg.path, ".aspell", dict.file)
##

##
## Create new personal dictionary
##
## R> words <- c("foo", "bar", "baz")  # Specify initial words
## R> saveRDS(sort(words), dict.path)

words <- c(
  ## Acronyms
    "DWT",              # Discrete Wavelet Transform
    "UDWT",             # Undecimated Discrete Wavelet Transform
  ## Parameters
    "Highpass",         # mirdwt.Rd:17:13, mrdwt.Rd:27:15
    "Lowpass",          # mirdwt.Rd:16:13, mrdwt.Rd:26:15
    "yl",               # mirdwt.Rd:20:36
  ## People
    "Daubechies",       # Ingrid Daubechies
    "Donoho",           # David Donoho
    "Heisenberg",       # Werner Heisenberg
    "Johnstone",        # Iain Johnstone
    "Kronecker",        # Leopold Kronecker
  ## Terms
    "deconvolution",
    "denoise",
    "denoised",
    "denoising",
    "periodize",
    "periodized",
    "undecimated",
  ## Software
    "WaveLab",          # MATLAB toolbox
    toupper(pkgname),   # Rice Wavelet Toolkit
    pkgname             # R package
)
#saveRDS(sort(words), dict.path)

##
## Add new words to existing dictionary
##
## R> words <- c("qux", "quux")        # Specify words to augment existing
## R> words <- c(readRDS(dict.path), words)
## R> saveRDS(sort(words), dict.path)
##

