
.. index:: Photospline
.. _Photospline:

Photospline
===========

Detector response to a high-energy physics process is often estimated by
Monte Carlo simulation. For purposes of data analysis, the results of this
simulation are typically stored in large multi-dimensional histograms,
which can quickly become unwieldy in terms of size or numerically
problematic due to unfilled bins or interpolation artifacts. Photospline is
a library that uses the `penalized spline technique`_ to efficiently
compute, store, and evaluate B-spline representations of such tables.

.. _`penalized spline technique`: http://dx.doi.org/10.1214/ss/1038425655

Fitting tutorial
================

.. highlightlang:: python

A basic fitting script starts with a few imports::
	
	import numpy
	from photospline import spglam as glam
	from photospline import splinefitstable
	from photospline.glam.glam import grideval
	from photospline.glam.bspline import bspline

Some of these are strictly required.

* ``numpy`` - We need this in order to pass n-dimensional arrays around
* ``glam`` - We will use the sparse-matrix implementation of the 
  least-squares fit. Replace this line with 
  ``from photospline.glam import glam`` in order to use the 
  inefficient (but readable) Python reference implementation.
* ``splinefitstable`` - We will use
  :py:func:`photospline.splinefitstable.write` to store the fit
  results in a FITS file.

We will also need some extra functions in order to visualize the result.

* ``grideval`` - We will use this function to evaluate the spline surface
  on a coordinate grid.
* ``bspline`` - We will also plot the individual basis functions to get a
  better idea of how the spline surface is constructed.

For this example, we will make up some random, one-dimensional data drawn
from a cosine overlaid on a polynomial::
	
	numpts = 500
	x1 = numpy.sort(numpy.random.uniform(-4,25,size=numpts))
	z = numpy.random.poisson(numpy.cos(x1)**2 + (x1-3.0)**2 + 10)

We can improve the quality of the fit by giving less weight to data points
that are likely to be distorted by statistical fluctuations. Here we'll use
a minimum weight of 1, and weight the high-occupancy points up::
	
	w = 1. + z

Now we need to choose a set of B-spline basis functions. We will use a
grid of order-2 B-splines that extend beyond the the data points, relying
on the regularization term for extrapolation::
	
	order = 2
	knots = [numpy.linspace(-8,35,30)]
	smooth = 3.14159e3

This generalizes easily to multi-dimensional data; we would simply add an
extra entry to ``knots`` for each extra dimension.

To actually run the fit, we call :py:func:`photospline.spglam.fit`  ::
	
	>>> result = glam.fit(z,w,[x1],knots,order,smooth)
	Calculating penalty matrix...
	Calculating spline basis...
	Reticulating splines...
		Convolving bases...
			Convolving dimension 0
		Flattening residuals matrix...
	Transforming fit array...
	Computing least square solution...
	Analyze[27]: 0.000049 s
	Factorize[27]: 0.000027 s
	Solve[27]: 0.000022 s
	Done: cleaning up

Now we can save the fit result for later use with :c:func:`ndsplineeval`
or :py:func:`photospline.glam.glam.grideval` ::
	
	splinefitstable.write(result, 'splinefit-1d.fits')

To see the result, we can plot it with `matplotlib`_::
	
	import pylab
	# Plot individual basis splines 
	xfine = numpy.linspace(knots[0][0], knots[0][-1], 10001)
	splines = [numpy.array([bspline(knots[0], x, n, order) for x in xfine]) for n in range(0,len(knots[0])-2-1)]
	for c, n in colorize(range(len(splines))):
		pylab.plot(xfine, result.coefficients[n]*splines[n], color=c)
	# Plot the spline surface (sum of all the basis functions)
	pylab.plot(xfine, glam.grideval(result, [xfine]), label='Spline fit', color='k')
	pylab.legend(loc='upper left')
	pylab.show()

.. _matplotlib: http://matplotlib.org/

The result is shown below:

.. figure:: splinefit_1d.png
	:width: 50%
	
	Example 1-dimensional spline fit. The spline surface (black 
	line) is the sum of basis functions (colored lines) weighted with
	coefficients determined by a linear least-squares fit to the data
	(blue dots). In the region to the right where there are no data,
	the order-2 penalty term produces a straight line.

Python library reference
========================

The interface to the fitting library is entirely in Python, using Numpy
arrays to as containers for the data, spline coefficients, knot vectors, 
etc. The spline coefficients determined in the fit are stored in FITS 
files that can be loaded and evaluated using the bundled C library.

.. note::
	The reference Python fitting implementation 
	:func:`photospline.glam.glam.fit` uses
	dense matrix operations from numpy, and can become 
	extremely memory-hungry in large numbers of dimensions.
		
	The C implementation :func:`photospline.spglam.fit`
	implements the same functions as the Python version, but
	uses sparse matrix operations from SuiteSparse_. This is
	both orders of magnitude faster and uses 100x less memory
	for typical problems. You should use spglam for large-scale
	fitting and fall back to glam if necessary for debugging.
	
.. _SuiteSparse: http://www.cise.ufl.edu/research/sparse/

.. autofunction:: photospline.glam.glam.fit

.. function:: photospline.spglam.fit(z, w, coords, knots, order, smooth=1, periods=None, penalties=None, monodim=None)
	
	A drop-in replacement for 
	:py:func:`photospline.glam.glam.fit`.

.. autofunction:: photospline.splinefitstable.write

.. autofunction:: photospline.splinefitstable.read

.. autofunction:: photospline.glam.glam.grideval

.. function:: photospline.spglam.grideval(table, coords)
	
	A drop-in replacement for 
	:py:func:`photospline.glam.glam.grideval`.
	
	.. note:: The ``bases`` argument is not supported in this version.

.. autofunction:: photospline.glam.glam.bspline

C library reference
===================

.. highlightlang:: c

The photospline C library contains a collection of functions for 
efficiently evaluating the tensor-product B-spline surfaces produced by
:py:func`photospline.spglam.fit` and written to disk as FITS
files by :py:func:`photospline.splinefitstable.write`.

Spline table I/O and manipulation
---------------------------------

.. c:type:: struct splinetable
	
	A data structure that holds the coefficients, knot grids,
	and associated metadata needed to evaluate a spline surface::
		
		struct splinetable {
			int ndim;
			int *order;
			double **knots;
			long *nknots;
			double **extents;
			double *periods;
			float *coefficients;
			long *naxes;
			unsigned long *strides;  
			int naux;
			char ***aux;
		};

.. c:function:: int readsplinefitstable(const char *path, struct splinetable *table)
	
	Read a spline table from a FITS file on disk.
	
	:param path: the filesystem path to read from
	:param splinetable: a pointer to a pre-allocated :c:type:`splinetable`
	:returns: 0 upon success

.. c:function:: int readsplinefitstable_mem(struct splinetable_buffer *buffer, struct splinetable *table)
	
	Read a spline table from a FITS file in memory.
	
	:param buffer: The memory buffer to read from. ``data`` must point
	               to the beginning of the memory area where the FITS
		       file is stored, and ``size`` should give its size.
		       ``mem_alloc`` and ``mem_realloc`` must be set to
		       valid addresses, but will not be called.
	
	Example::
		
		struct splinetable table;
		struct splinetable_buffer buf;
		
		readsplinefitstable("foo.fits", &table);
		buf.mem_alloc = &malloc;
		buf.mem_realloc = &realloc;
		writesplinefitstable_mem(&buf, &table);
		splinetable_free(&table);
		
		readsplinefitstable_mem(&buf, &table);
		free(buf.data);
		/* do stuff with table */
		splinetable_free(&table);

.. c:function:: int writesplinefitstable(const char *path, const struct splinetable *table)
	
	Write a spline table to a FITS file on disk.

.. c:function:: int writesplinefitstable_mem(struct splinetable_buffer *buffer, const struct splinetable *table)
	
	Write a spline table to a FITS file in memory.
	
	:param buffer: the memory buffer to write to. Memory will be
	               allocated internally via ``mem_alloc`` and
	               ``mem_realloc``, which should point to
	               :c:func:`malloc` and :c:func:`realloc` or
	               functional equivalents.
	
	Example::
		
		struct splinetable table;
		struct splinetable_buffer buf;
		
		readsplinefitstable("foo.fits", &table);
		buf.mem_alloc = &malloc;
		buf.mem_realloc = &realloc;
		writesplinefitstable_mem(&buf, &table);
		splinetable_free(&table);
		/* do stuff with buf */
		free(buf.data);

.. c:type: struct splinetable_buffer

	A structure describing a memory area and how to resize it::
		struct splinetable_buffer {
			void *data;
			size_t size;
			void *(*mem_alloc)(size_t newsize);
			void *(*mem_realloc)(void *p, size_t newsize);
		};

.. c:function:: char * splinetable_get_key(struct splinetable *table, const char *key)

	Get a pointer to the raw string representation of an auxiliary
	field stored in the FITS header.

	:param key: A null-terminated string containing the key.
	:returns: A pointer to the beginning of the value string, or NULL
	          if the key doesn't exist.

.. c:function:: int splinetable_read_key(struct splinetable *table, splinetable_dtype type, const char *key, void *result)
	
	Parse and store the value of an auxiliary field stored in the FITS
	header.
	
	:param type: The type of value stored.
	:param key: A null-terminated string containing the key.
	:param result: Where to store the result
	:returns: 0 upon success

.. c:type:: splinetable_dtype

	Types that can be parsed by :c:func:`splinetable_read_key`::
		
		typedef enum {
			SPLINETABLE_INT,
			SPLINETABLE_DOUBLE
		} splinetable_dtype;

.. c:function: int splinetable_convolve(struct splinetable *table, const int dim, const double *knots, size_t n_knots)
	
	Convolve a table with the spline defined on a set of knots
	along a given dimension and store the spline expansion of the
	convolved surface in the table. This will raise the order of the
	splines in the given dimension by (n_knots - 1), i.e. convolving
	with an order-0 spline (a box function, defined on two knots) will
	raise the order of the spline surface by 1.
	
	:param dim: The dimension along which to apply the convolution
	:param knots: The knots defining the convolution kernel spline.
	:param n_knots: The number of knots. The order of the kernel spline
	                will be n_knots-2.
	:returns: 0 upon success

.. c:function:: void splinetable_free(struct splinetable *table)
	
	Free memory allocated by :c:func:`readsplinefitstable`.

Spline evaluation
-----------------

.. c:function:: int tablesearchcenters(const struct splinetable *table, const double *x, int *centers)
	
	Find the index of the order-0 spline in each dimension that has
	support at the corresponding coordinate.
	
	:param x: an array of coordinates, one for each dimension
	:param centers: the array where the index corresponding to each
	                coordinate will be stored.
	:returns: 0 upon success. A non-zero return value indicates that
	          one or more of the coordinates is outside the
	          support of the spline basis in its dimension, in which
	          case the value of the spline surface is 0 by 
	          construction.


.. c:function:: double ndsplineeval(const struct splinetable *table, const double *x, const int *centers, int derivatives)
	
	Evaluate the spline surface or its derivatives.
	
	:param x: an array of coordinates, one for each dimension
	:param centers: an array of central-spline indices in each
	                dimension, filled in a previous call to
	                :c:func:`tablesearchcenters`
	:param derivatives: a bitmask indicating the type of basis to use
	                    in each dimension. If the bit corresponding to
	                    a dimension is set, the basis in that
	                    dimension will consist of the derivatives of
	                    the usual B-spline basis, and the return value
	                    will be the gradient of the surface in that
	                    dimension.
	
	For example, passing an unset bitmask evaluates the surface::
		
		struct splinetable table;
		int centers[table.ndim];
		double x[table.ndim];
		double v;
		
		if (tablesearchcenters(&table, &x, &centers) != 0)
			v = ndsplineeval(&table, &x, &centers, 0);
		else
			v = 0;

	while setting individual bits evalulates the elements of the gradient::
		
		double gradient[table.ndim];
		for (i = 0; i < table.ndim; i++)
			gradient[i] = ndsplineeval(&table, &x, &centers, (1 << i));

.. c:function:: ndsplineeval_gradient(const struct splinetable *table, const double *x, const int *centers, double *evaluates)
	
	Evaluate the spline surface and all of its derivatives, using
	SIMD operations to efficiently multiply the same coefficient matrix
	into multiple bases. 
	
	:param x: an array of coordinates, one for each dimension
	:param centers: an array of central-spline indices in each
	                dimension, filled in a previous call to
	                :c:func:`tablesearchcenters`
	:param evaluates: an array of size at least ndim+1 where the value
	                  of the surface and the elements of the gradient
	                  will be stored
	
	On most platforms, each group of 4 bases (e.g. the value and first
	3 elements of the gradient) takes the nearly the same number of
	operations as a single call to :c:func:`ndsplineeval`. As a
	result, this version is much faster than sequential calls to
	:c:func:`ndsplineeval` in applications like maximum-likelihood
	fitting where both the value and entire gradient are required.
