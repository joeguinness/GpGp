
## Test environments

* local: ubuntu 18.10, R 4.0.2
* r-hub: Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* r-hub: macOS 10.13.6 High Sierra, R-release, CRAN's setup
* r-hub: 


## R CMD check results

There were no ERRORs or WARNINGs

There was 1 NOTEs:

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


