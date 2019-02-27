# cxbtools
Tools to deal with (unresolved) point sources from the Cosmic X-ray Background (CXB) in soft X-ray spectra.

The CXB tools are written to help determining the contribution from unresolved X-ray sources to the 
background spectrum in X-ray spectra. They were especially written for the EPIC instruments aboard 
XMM-Newton. Please cite [Mernier et al. (2015)](http://dx.doi.org/10.1051/0004-6361/201425282) 
and [Zenodo 2575495](http://doi.org/10.5281/zenodo.2575495) if you use this software for a scientific publication.

See the reference or doc/cxb.pdf for an explanation of the method used in this software.

## Install

This installation requires the CMake package to be installed on your system. On Linux, this
can be installed through the package repository of your distribution. For MacOS, install packages
can be downloaded from: [CMake.org](https://cmake.org).

In addtion, the code depends on the GNU Scientific Library (GSL). This library can also be installed 
through the package manager of your distribution. On MacOS, this package needs to be compiled 
separately. For more information: [GSL website](https://www.gnu.org/software/gsl/). 

For the cxbtools compilation, download or clone the source code from github:

```
  user@unix:~> git clone https://github.com/jdeplaa/cxbtools.git
```

Then change to the cxbtools directory and run:

```
  user@unix:~/cxbtools> cmake .
  user@unix:~/cxbtools> make 
```

If you want to install cxbtools in a system location, you can OPTIONALLY do:

```
  user@unix:~/cxbtools> sudo make install
```

The excutables for the tools described in the next section can be found in the source code directory.
See the CMake documentation for more install and compilation options.

## Using CXBTools

There are three tools included:

### cxbups
CXBups is a program to calculate the remaining Cosmic X-ray Background
flux after the exclusion of point sources up to a certain limit.

  Usage:
    ./cxbups flux-limit area

  Example:
    ./cxbups 3.E-15 0.0549

The flux-limit is in erg/cm^2/s
The area is in square degrees.

### cxbrnd
CXBrnd is a program to calculate the remaining Cosmic X-ray Background
flux after the exclusion of point sources up to a certain limit through
Monte Carlo simulations.

  Usage:
    ./cxbrnd flux-limit area niter

  Example:
    ./cxbrnd 3.E-15 0.0549 1000

The flux-limit is in erg/cm^2/s
The area is in square degrees.

### cxbopt
CXBopt is a program to calculate optimal source extraction radius
and flux cut for a certain annulus.

  Usage:
    ./cxbopt source-cnts backg-counts region

  Example:
    ./cxbopt 10000. 15000. 0.0549

The region area is in square degrees.


