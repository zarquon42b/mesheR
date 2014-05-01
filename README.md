mesheR
====
__mesheR__ is an R-package providing methods for mesh manipulations in R. Its main feature are advanced surface registration algorithms.


##### Installation of the R package "mesheR": ####
   0. Make sure to work with the latest version of R and install dependencies  
     1. Make sure, you already have [Morpho](http://sourceforge.net/p/morpho-rpackage/wiki/Installation_Morpho/) and [Rvcg](http://sourceforge.net/p/morpho-rpackage/wiki/Installation_Rvcg/) installed.
     2. install additional dependencies from CRAN (type the following commands into your R terminal):: 
               
               setRepositories(ind=1:2)
               install.packages(c("spam","colorspace","rhdf5"))


   1. Download the version suitable for your OS Either the compiled package (for Windows and OS X) or the source tarball (Linux) from [here](https://sourceforge.net/projects/morpho-rpackage/files/mesheR/).

   2. Installation command from within R: 
   
        install.packages("Path_to_downloaded_package_mesheR[Version_OS]",repos=NULL)

   3. check if the package can be loaded:
        
        load package: library(mesheR)

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
	setRepositories(ind=1:2)	
	install_github("zarquon42b/mesheR")
