
The major purpose of this release is to fix an error
on CRAN checks for Oracle Developer Studio. We have also
parallelized the computation of the inverse cholesky
factor, and included safeguards against large smoothness
parameters in matern_isotropic.


## Test environments

* local: ubuntu 18.10, R 4.0.3
* r-hub: Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* r-hub: macOS 10.13.6 High Sierra, R-release, CRAN's setup
* r-hub: Oracle Solaris 10, x86, 32 bit, R-release, Oracle Developer Studio 12.6


## R CMD check results

There were no ERRORs or WARNINGs

There was 1 NOTEs:

    N  checking installed package size ...
         installed size is 10.8Mb
         sub-directories of 1Mb or more:
           data   1.4Mb
           libs   9.1Mb


## Downstream dependencies

There is one downstream dependency (GPvecchia), but it
uses only 1 function that is not affected by this
package update.


