#ifndef PHOTOSPLINE_CORE_BSPLINE_H
#define PHOTOSPLINE_CORE_BSPLINE_H

#include <cstddef>

namespace photospline{

/*
 * Compute the value of the ith nth-order basis spline of a set
 * defined by knots at the point x.
 */

double bspline(const double* knots, double x, int i, int n);
double bspline_deriv(const double* knots, double x, int i, int n, unsigned order);

/*
 * A brain-dead reimplementation of de Boor's BSPLVB, which generates
 * the values of the non-zero B-splines at x from the bottom up without
 * unnecessarily recalculating terms. 
 * 
 * NB: for bsplvb_simple(), bspline_nonzero(), and bspline_deriv_nonzero(),
 * `left' must be the index of the nearest fully-supported knot
 * span for splines of order n, i.e. n <= left <= nknots-n-2. For bsplvb(),
 * `left' must be the index of the nonzero 0-th order spline, i.e.
 * knots[left] <= x < knots[left+1].
 *
 * See Chapter X in: 
 * 
 * Carl de Boor. A Practical Guide to Splines, volume 27 of Applied
 *     Mathematical Sciences. Springer-Verlag, 1978.
 */

void bsplvb_simple(const double* knots, const unsigned nknots,
    double x, int left, int jhigh, float* biatx);
void bsplvb(const double* knots, const double x, const int left, const int jlow,
    const int jhigh, float* biatx,
    double* delta_l, double* delta_r);
void bspline_nonzero(const double* knots, const unsigned nknots,
    const double x, int left, const int n, float* values, float* derivs);
void bspline_deriv_nonzero(const double* knots, const unsigned nknots,
    const double x, const int left, const int n, float* biatx);

/*
 * Evaluates the results of a full spline basis given a set of knots,
 * a position, an order, and a central spline for the position (or -1).
 * The central spline should be the index of the 0th order basis spline
 * that is non-zero at the position x.
 */

double splineeval(const double* knots, const double* weights, int nknots, double x,
    int order, int center);
	
} //namespace photospline

#endif /* PHOTOSPLINE_CORE_BSPLINE_H */

