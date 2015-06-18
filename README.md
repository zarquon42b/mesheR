mesheR [![Travis Build Status](https://travis-ci.org/zarquon42b/mesheR.png?branch=master)](https://travis-ci.org/zarquon42b/mesheR)

====
__mesheR__ is an R-package providing methods for mesh manipulations in R. Its main feature are advanced surface registration algorithms.

##### Warning: #####
This code is still under development, the API may change without further notice, so see commit log in case your old code does not work any more.


##### Installation of the R-package "mesheR" (latest development code) using *devtools*: ####

0. Make sure, you already have [Morpho](http://sourceforge.net/p/morpho-rpackage/wiki/Installation_Morpho/) and [Rvcg](http://sourceforge.net/p/morpho-rpackage/wiki/Installation_Rvcg/) installed.
1. install *devtools* from within R (Ubuntu/Debian users will have to install *libcurl4-gnutls-dev* beforehand):

		install.packages("devtools")


2. Install build environment
    * **Windows:** Install latest version of *[Rtools](http://cran.r-project.org/bin/windows/Rtools)*.
During installation of *Rtools* make sure to install the *toolchain*, and to select *"Edit the system path"* (and confirming the installers suggestions).
    * **OSX:** Install *[XCODE](https://developer.apple.com/xcode/)*


3. In R run the command:

	
		require(devtools)
		install_github("zarquon42b/mesheR")
