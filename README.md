# Photospline

Detector response to a high-energy physics process is often estimated by Monte
Carlo simulation. For purposes of data analysis, the results of this simulation
are typically stored in large multi-dimensional histograms, which can quickly
become unwieldy in terms of size or numerically problematic due to unfilled
bins or interpolation artifacts. Photospline is a library that uses the
[penalized spline technique](http://dx.doi.org/10.1214/ss/1038425655) to
efficiently compute, store, and evaluate B-spline representations of such
tables.

## Installing

To build the core photospline libraries, you will need:

* CMake >= 3.1
* A C++11-compliant compiler
* cfitsio

This will allow you to evaluate splines, but not much else. The fitting library further requires:

* SuiteSparse
* BLAS and LAPACK libraries (we recommend OpenBLAS)

To build and install in /usr/local, run `cmake . -DCMAKE_INSTALL_PREFIX=/usr/local` followed by `make install` in the source directory.
