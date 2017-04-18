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

* [CMake](https://cmake.org) >= 3.1
* A C++11-compliant compiler (e.g. [gcc](https://gcc.gnu.org) >= 4.8.1 or [clang](https://clang.llvm.org) >= 3.3)
* [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html)

This will allow you to evaluate splines, but not much else. The fitting library further requires:

* [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html)
* BLAS and LAPACK libraries (we recommend [OpenBLAS](http://www.openblas.net) on Linux; the built-in Accelerate.framework is fine on OS X)

To build and install in /usr/local, run `cmake . -DCMAKE_INSTALL_PREFIX=/usr/local`, followed by `make install`, in the source directory.

## Documentation

A tutorial on fitting spline surfaces with photospline is given in the [sphinx
documentation](docs/source/index.rst). If you have
[sphinx](http://www.sphinx-doc.org/en/stable/),
[doxygen](http://www.doxygen.org/), and
[breate](http://breathe.readthedocs.io/en/latest/) (a sphinx-doxygen bridge),
you can build a pretty-looking HTML version, complete with API documentation,
with `make html`. The output will appear in the build directory under
`docs/index.html`.

## Using photospline in downstream projects

Photospline exports its build configuration for use in downstream projects
that are built with CMake. The files are stored in `CMAKE_INSTALL_PREFIX/share/photospline/cmake`,
so as long as the install prefix is in CMake's prefix path (which you can
specify with [`CMAKE_PREFIX_PATH`](https://cmake.org/cmake/help/v3.0/command/find_package.html)), compiling a C++ executable against photospline is as simple as

    find_package(photospline REQUIRED)
    
    add_executable(foo foo.cxx)
    target_link_libraries(foo photospline)
,

or for a C executable,
    
    add_executable(cfoo foo.c)
    target_link_libraries(cfoo cphotospline)
.

This causes the targets `foo` and `cfoo` to be built with the include paths,
compiler options, and link libraries specified when the photospline libraries
were built.

