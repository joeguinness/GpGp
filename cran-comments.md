
## Test environments

* local ubuntu 18.10, R 4.0.2
* Debian Linux, R-devel, GCC ASAN/UBSAN
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Ubuntu Linux 16.04 LTS, R-release, GCC
* Fedora Linux, R-devel, clang, gfortran


## R CMD check results

There were no ERRORs or WARNINGs

There were 2 NOTEs:

* note about package size. I have made
  an effort to reduce the package size.

    > checking installed package size ... NOTE
        installed size is 10.6Mb
    	sub-directories of 1Mb or more:
    	  data    1.4Mb
    	  libs    8.8Mb
    
* note about a URL in documentation. This url
  works when I check it in my browser.

    > checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Joseph Guinness '
    Found the following (possibly) invalid URLs:
    URL: http://www.jstor.org/stable/2345768
    From: DESCRIPTION
    Status: 403
    Message: Forbidden


## Downstream dependencies

There is one downstream dependency (GPvecchia), but it
uses only 1 function that is not affected by this
package update.


