/*
 * bspline_multi.c: Provides efficient routines using vector intrinsics to
 *    simultaneously compute the value and gradient of an N-dimensional
 *    B-spline surface.
 */

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cassert>

#include "photospline/bspline.h"
#include "photospline/detail/simd.h"

namespace photospline{

void
bspline_nonzero(const double* knots, const unsigned nknots,
    const double x, int left, const int n,
    float* values, float* derivs)
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
	assert(left >= n && left <= nknots-n-2);
	if (left == n)
		while (left >= 0 && x < knots[left])
			left--;
	else if (left == nknots-n-2)
		while (left < nknots-1 && x > knots[left+1])
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

} //namespace photospline
