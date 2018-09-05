/*
 * bspline.c: Routines for calculating values of B-splines and their
 *  derivatives, as well as efficient computation of the values of
 *  N-dimensional B-spline surfaces.
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>

//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_blas.h>

#include <photospline/bspline.h>

namespace photospline{

/*
 * Compute the value of the ith nth-order basis spline of a set
 * defined by knots at the point x.
 *
 * This is implemented using the De Boor algorithm, as outlined on
 * Wikipedia.
 */

double
bspline(const double *knots, double x, int i, int n)
{
	double result;

	if (n == 0) {
		/*
		 * Special case the 0th order case, where B-Splines
		 * are constant functions from one knot to the next.
		 */

		if (x >= knots[i] && x < knots[i+1])
			return 1.0;
		else
			return 0.0;
	}

	result = (x - knots[i])*bspline(knots, x, i, n-1) /
	    (knots[i+n] - knots[i]);
	result += (knots[i+n+1] - x)*bspline(knots, x, i+1, n-1) /
	    (knots[i+n+1] - knots[i+1]);

	return result;
}

double
bspline_deriv(const double *knots, double x, int i, int n, unsigned order)
{
	double result;

	if (n == 0) {
		/*
		 * Special case the 0th order case, where B-Splines
		 * are constant functions from one knot to the next.
		 */

		return 0.0;
	}

	if (order <= 1) {
		result = n * bspline(knots, x, i, n-1) / (knots[i+n] - knots[i]);
		result -= n * bspline(knots, x, i+1, n-1) / (knots[i+n+1] - knots[i+1]);
	} else {
		result = n * bspline_deriv(knots, x, i, n-1, order-1) / (knots[i+n] - knots[i]);
		result -= n * bspline_deriv(knots, x, i+1, n-1, order-1) / (knots[i+n+1] - knots[i+1]);
	}
	
	return result;
}

/*
 * Evaluates the results of a full spline basis given a set of knots,
 * a position, an order, and a central spline for the position (or -1).
 * The central spline should be the index of the 0th order basis spline
 * that is non-zero at the position x.
 */

double
splineeval(const double *knots, const double *weights, int nknots, double x, int order,
    int center)
{
	double work = 0.0;
	int i;

	if (center < 0) {
		/* XXX: should be a binary search */
		for (center = 0; center+1 < nknots; center++) {
			if (x > knots[center] && x < knots[center+1])
				break;
		}
	
		if (center+1 >= nknots)
			return 0.0;
	}

	i = center - order;
	if (i < 0)
		i = 0;

	while (i < nknots-order-1 && i <= center) {
		work += weights[i]*bspline(knots, x, i, order);
		i++;
	}

	return work;
}

// What is this?
//double
//ndsplineeval_linalg(const struct splinetable *table, const double *x,
//    const int *centers, int derivatives)
//{
//	int totalcoeff, n;
//	int coeffstrides[table->ndim];
//	gsl_matrix_float *basis1, *basis2, *basis_elem;
//
//	coeffstrides[table->ndim - 1] = totalcoeff = 1;
//        for (n = table->ndim-1; n >= 0; n--) {
//                totalcoeff *= (table->order[n] + 1);
//                if (n > 0)
//                        coeffstrides[n-1] = totalcoeff;
//        }
//
//	float basis1_data[totalcoeff], basis2_data[totalcoeff],
//	    elem_data[maxorder(table->order, table->ndim) + 1];
//	gsl_matrix_float b1, b2, be;
//	basis1 = &b1; basis2 = &b2; basis_elem = &be;
//	basis1->data = basis1_data;
//	basis2->data = basis2_data;
//	basis_elem->data = elem_data;
//
//	/*
//	 * Form outer product basis1 = basis2 x basis_elem, filling basis_elem
//	 * every time with the non-zero basis functions on each axis and
//	 * swapping basis1 and basis2 via tmp_basis.
//	 */
//	basis2->size1 = 1;
//	basis2->size2 = 1;
//	basis2->data[0] = 1.0;
//
//	for (n = table->ndim-1; n >= 0; n--) {
//		gsl_matrix_float *tmp_basis;
//		if (derivatives & (1 << n)) {
//			bspline_deriv_nonzero(table->knots[n], 
//			    table->nknots[n], x[n], centers[n],
//			    table->order[n], basis_elem->data);
//		} else {
//			bsplvb_simple(table->knots[n], table->nknots[n],
//			    x[n], centers[n], table->order[n] + 1,
//			    basis_elem->data);
//		}
//
//		basis_elem->size1 = table->order[n] + 1;
//		basis_elem->size2 = 1;
//
//		basis1->size2 = basis2->size2;
//		basis1->size1 = basis_elem->size1;
//		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
//		    basis1->size1, basis1->size2, basis_elem->size2, 1.0,
//		    basis_elem->data, 1, basis2->data, basis2->size2, 0,
//		    basis1->data, basis1->size2);
//		basis1->size2 = basis1->size1 * basis1->size2;
//		basis1->size1 = 1;
//
//		tmp_basis = basis1;
//		basis1 = basis2;
//		basis2 = tmp_basis;
//	}
//
//	/* Now basis1 is free, so fill it with the spline coefficients */
//	int i, tablepos;
//	int decomposedposition[table->ndim];
//	tablepos = 0;
//	for (n = 0; n < table->ndim; n++) {
//		decomposedposition[n] = 0;
//		tablepos += (centers[n] - table->order[n])*table->strides[n];
//	}
//
//	for (i = 0; i < table->order[table->ndim-1] + 1; i++)
//		basis1->data[i] = table->coefficients[tablepos + i];
//	for (n = 1; n < coeffstrides[0] /* number of chunks */; n++) {
//		tablepos += table->strides[table->ndim-2];
//		decomposedposition[table->ndim-2]++;
//		/* Carry to higher dimensions */
//		for (i = table->ndim-2; __builtin_expect(i > 0 && 
//		    decomposedposition[i] > table->order[i], 0); i--) {
//			decomposedposition[i-1]++;
//			tablepos += (table->strides[i-1]
//			    - decomposedposition[i]*table->strides[i]);
//			decomposedposition[i] = 0;
//		}
//
//		for (i = 0; i < table->order[table->ndim-1] + 1; i++)
//			basis1->data[n*(table->order[table->ndim-1] + 1) + i] =
//			    table->coefficients[tablepos + i];
//	}
//
//	/* Take the dot product */
//	__builtin_prefetch(basis1->data);
//	__builtin_prefetch(basis2->data);
//	return cblas_sdot(totalcoeff, basis1->data, 1, basis2->data, 1);
//}
	
} //namespace photospline
