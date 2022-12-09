#This package and function will create a repository with all R and package versions used in this R project. In another computer it will retrieve the one created in mine and allow a fully reproducible R code.
library(checkpoint)
library(withr)
library(pkgcache)
library(pkgdepends)
#getRversion() '4.0.5'
checkpoint("2022-09-18",checkpoint_location=getwd(),r_version = '4.0.5' )
#check they work
getOption("repos")
.libPaths()
View(installed.packages(.libPaths()[1])[, "Package"])
#END-----


