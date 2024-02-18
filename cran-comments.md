
Bug fixes related to automatic coercion of matrices to vectors
Workaround for change to roxygen that removed package-alias
Removed maptools from Suggests field of DESCRIPTION

## Test environments

* local: ubuntu 20.04, R 4.3.2
* Win-builder
* r-hub: Ubuntu Linux 20.04.1 LTS, R-release, GCC
* r-hub: Fedora Linux, R-devel, clang, gfortran


## R CMD check results

There were no ERRORs or WARNINGs

There were 3 NOTEs:

installed size is 13.7Mb
   sub-directories of 1Mb or more:
      data   1.4Mb
      libs  11.9Mb

Found the following (possibly) invalid URLs:
   URL: http://www.jstor.org/stable/2345768
      From: DESCRIPTION
      Status: 403
      Message: Forbidden

Found the following (possibly) invalid URLs:
   URL: https://www.ncei.noaa.gov/products/jason-satellite-products
      From: man/jason3.Rd
      Status: Error
      Message: Empty reply from server

Both urls are valid. Possibly jstor rejects requests that are not from browsers.

