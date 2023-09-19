#!/usr/bin/env python

import sys, os
sys.path.append(os.getcwd())
import numpy
import photospline

def pad_knots(knots, order=2):
	"""
	Pad knots out for full support at the boundaries
	"""
	pre = knots[0] - (knots[1]-knots[0])*numpy.arange(order, 0, -1)
	post = knots[-1] + (knots[-1]-knots[-2])*numpy.arange(1, order+1)
	return numpy.concatenate((pre, knots, post))

# reproducibility considered good
numpy.random.seed(42)

z, edges = numpy.histogramdd(numpy.random.multivariate_normal([0,0,0], numpy.eye(3), size=1000),
    bins=[numpy.linspace(-3, 3, 101)]*3)
z = z.cumsum(axis=2).astype(float)
w = numpy.ones(z.shape)

# evaluate CDF at right-hand bin edges
centers = [0.5*(edges[i][:-1] + edges[i][1:]) for i in range(2)] + [edges[-1][1:]]
knots = [pad_knots(numpy.linspace(-3, 3, 5))]*3
order = [2,2,3]
smooth = 1

data, ws = photospline.ndsparse.from_data(z, w)
spline = photospline.glam_fit(data, ws, centers, knots, order, [smooth]*3, [2]*3, monodim=2, verbose=False)
y = spline.grideval(centers)

residual = ((z-y)**2).sum()
numpy.testing.assert_almost_equal(residual, 50791.31, 2)
