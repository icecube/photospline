#ifndef PHOTOSPLINE_CORE_BSPLINE_H
#define PHOTOSPLINE_CORE_BSPLINE_H

#include <cstddef>
#include <cassert>

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

template <typename Float=float>
void bsplvb_simple(const double* knots, const unsigned nknots,
    double x, int left, int degree, Float* biatx)
{
	int i, j;
	double saved, term;
	double delta_l[degree], delta_r[degree];
	
	biatx[0] = 1.0;
	
	/*
	 * Handle the (rare) cases where x is outside the full
	 * support of the spline surface.
	 */
	if (left == degree-1)
		while (left >= 0 && x < knots[left])
			left--;
	else if (left == int(nknots)-degree-1)
		while (left < int(nknots)-1 && x > knots[left+1])
			left++;	
	
	/* 
	 * NB: if left < degree-1 or left > nknots-degree-1,
	 * the following loop will dereference addresses ouside
	 * of knots[0:nknots]. While terms involving invalid knot
	 * indices will be discarded, it is important that `knots'
	 * have (maxdegree-1)*sizeof(double) bytes of padding
	 * before and after its valid range to prevent segfaults
	 * (see parsefitstable()).
	 */
	for (j = 0; j < degree-1; j++) {
		delta_r[j] = knots[left+j+1] - x;
		delta_l[j] = x - knots[left-j];
		
		saved = 0.0;
		
		for (i = 0; i < j+1; i++) {
			term = biatx[i] / (delta_r[i] + delta_l[j-i]);
			biatx[i] = saved + delta_r[i]*term;
			saved = delta_l[j-i]*term;
		}
		
		biatx[j+1] = saved;
	}
	
	/* 
	 * If left < (spline order), only the first (left+1)
	 * splines are valid; the remainder are utter nonsense.
	 */
	if ((i = degree-1-left) > 0) {
		for (j = 0; j < left+1; j++)
			biatx[j] = biatx[j+i]; /* Move valid splines over. */
		for ( ; j < degree; j++)
			biatx[j] = 0.0; /* The rest are zero by construction. */
	} else if ((i = left+degree+1-nknots) > 0) {
		for (j = degree-1; j > i-1; j--)
			biatx[j] = biatx[j-i];
		for ( ; j >= 0; j--)
			biatx[j] = 0.0;
	}
}

template <typename Float=float>
void bsplvb(const double* knots, const double x, const int left, const int jlow,
    const int jhigh, Float* biatx,
    double* delta_l, double* delta_r)
{
	int i, j;
	double saved, term;

	if (jlow == 0)
		biatx[0] = 1.0;
		
	for (j = jlow; j < jhigh-1; j++) {
		delta_r[j] = knots[left+j+1] - x;
		delta_l[j] = x - knots[left-j];
		
		saved = 0.0;
		
		for (i = 0; i < j+1; i++) {
			term = biatx[i] / (delta_r[i] + delta_l[j-i]);
			biatx[i] = saved + delta_r[i]*term;
			saved = delta_l[j-i]*term;
		}
		
		biatx[j+1] = saved;
	}
}

template <typename Float=float>
void bspline_nonzero(const double* knots, const unsigned nknots,
    const double x, int left, const int n, Float* values, Float* derivs)
{
	/* Special case for constant splines */
	if (n == 0) {
		values[0] = 1;
		derivs[0] = 0;
		return;
	}
	
	/*
	 * Handle the (rare) cases where x is outside the full
	 * support of the spline surface.
	 */
	assert(left >= n && left <= int(nknots)-n-2);
	if (left == n)
		while (left >= 0 && x < knots[left])
			left--;
	else if (left == int(nknots)-n-2)
		while (left < int(nknots)-1 && x > knots[left+1])
			left++;
	
	double delta_r[n+1], delta_l[n+1];
	
	/* Get the non-zero n-1th order B-splines at x */
	bsplvb(knots, x, left, 0, n, values, delta_r, delta_l);
	
	/* 
	 * Now, form the derivatives of the nth order B-splines from
	 * linear combinations of the lower-order splines.
	 *
	 * NB: bspline_deriv_nonzero() uses double-precision
	 *     temporaries, so we do the same here to ensure that
	 *     the results are identical.
	 */
	
	/* 
	 * On the last supported segment of the ith nth order spline,
	 * only the i+1th n-1th order spline is nonzero.
	 */
	double temp = values[0];
	derivs[0] =  - n*temp / ((knots[left+1] - knots[left+1-n]));
	/* On the middle segments, both the ith and i+1th splines contribute. */
	int i, j;
	for (i = 1; i < n; i++) {
		double a = n*temp/((knots[left+i] - knots[left+i-n]));
		temp = values[i];
		derivs[i] = a - n*temp/(knots[left+i+1] - knots[left+i+1-n]);
	}
	/*
	 * On the first supported segment of the i+nth nth order spline,
	 * only the ith n-1th order spline is nonzero.
	 */
	derivs[n] = n*temp/((knots[left+n] - knots[left]));
	
	/* Now, continue to the non-zero nth order B-splines at x */
	bsplvb(knots, x, left, n-1, n+1, values, delta_r, delta_l);
	
	/* Rearrange for partially-supported points. */
	if ((i = n-left) > 0) {
		for (j = 0; j < left+1; j++) {
			values[j] = values[j+i]; /* Move valid splines over. */
			derivs[j] = derivs[j+i];
		}
		for ( ; j < n+1; j++)
			values[j] = derivs[j] = 0.0;
	} else if ((i = left+n+2-nknots) > 0) {
		for (j = n; j > i-1; j--) {
			values[j] = values[j-i];
			derivs[j] = derivs[j-i];
		}
		for ( ; j >= 0; j--)
			values[j] = derivs[j] = 0.0;
	}
}

template <typename Float=float>
void bspline_deriv_nonzero(const double* knots, const unsigned nknots,
    const double x, int left, const int n, Float* biatx)
{
	int i, j;
	double temp, a;
	double delta_l[n], delta_r[n];
	
	/* Special case for constant splines */
	if (n == 0)
		return;
	
	/*
	 * Handle the (rare) cases where x is outside the full
	 * support of the spline surface.
	 */
	if (left == n)
		while (left >= 0 && x < knots[left])
			left--;
	else if (left == int(nknots)-n-2)
		while (left < int(nknots)-1 && x > knots[left+1])
			left++;
	
	/* Get the non-zero n-1th order B-splines at x */
	bsplvb(knots, x, left, 0 /* jlow */, n /* jhigh */,
	    biatx, delta_l, delta_r);
	
	/* 
	 * Now, form the derivatives of the nth order B-splines from
	 * linear combinations of the lower-order splines.
	 */
	
	/* 
	 * On the last supported segment of the ith nth order spline,
	 * only the i+1th n-1th order spline is nonzero.
	 */
	temp = biatx[0];
	biatx[0] =  - n*temp / ((knots[left+1] - knots[left+1-n]));
	
	/* On the middle segments, both the ith and i+1th splines contribute. */
	for (i = 1; i < n; i++) {
		a = n*temp/((knots[left+i] - knots[left+i-n]));
		temp = biatx[i];
		biatx[i] = a - n*temp/(knots[left+i+1] - knots[left+i+1-n]);
	}
	/*
	 * On the first supported segment of the i+nth nth order spline,
	 * only the ith n-1th order spline is nonzero.
	 */
	biatx[n] = n*temp/((knots[left+n] - knots[left]));

	/* Rearrange for partially-supported points. */
	if ((i = n-left) > 0) {
		for (j = 0; j < left+1; j++)
			biatx[j] = biatx[j+i]; /* Move valid splines over. */
		for ( ; j < n+1; j++)
			biatx[j] = 0.0; /* The rest are zero by construction. */
	} else if ((i = left+n+2-nknots) > 0) {
		for (j = n; j > i-1; j--)
			biatx[j] = biatx[j-i];
		for ( ; j >= 0; j--)
			biatx[j] = 0.0;
	}
}

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

