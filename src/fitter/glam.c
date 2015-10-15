
#include "cholesky_solve.h"
#include "photospline/detail/splineutil.h"
#include "photospline/detail/glam.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

static cholmod_sparse* flatten_ndarray_to_sparse(struct ndsparse *array,
    size_t nrow, size_t ncol, cholmod_common* c);
cholmod_sparse* calc_penalty(uint64_t* nsplines, double *knots, uint32_t ndim, uint32_t i,
    uint32_t order, uint32_t porder, int mono, cholmod_common* c);

void print_ndsparse_py(struct ndsparse *a);

#define max(a,b) ((a > b) ? a : b)

int
glamfit_complex(const struct ndsparse* data, const double* weights, const double* const* coords,
    uint32_t ndim, const uint64_t* nknots, const double* const* knots,
	const uint64_t* naxes, float* out_coefficients, const uint32_t* order,
	cholmod_sparse* penalty, uint32_t monodim, int verbose, cholmod_common* c)
{
	cholmod_sparse** bases, ** boxedbases;
	cholmod_dense*coefficients, * Rdens;
	struct ndsparse F, R;
	cholmod_sparse* Fmat, * Rmat, * fitmat;
	double scale1[2] = {1.0, 0.0}, scale2[2] = {1.0, 0.0};
	size_t sidelen;
	uint64_t* nsplines;
	long i, j, k;
	int err = 0;
	
	assert(data->ndim>0);
	nsplines = calloc(data->ndim,sizeof(uint64_t));

	/*
	 * We will start by precomputing some things before we don't need to
	 * worry about memory pressure on temporaries.
	 */

	/* Figure out the total number of spline nodes */
	sidelen = 1;
	for (i = 0; i < ndim; i++) {
		nsplines[i] = nknots[i] - order[i] - 1;
		sidelen *= nsplines[i];
	}

	/*
	 * Compute the spline basis matrices, and pre-box them.
	 */

	if (verbose)
		printf("Calculating spline basis...\n");

	bases = calloc(ndim, sizeof(cholmod_sparse*));
	boxedbases = calloc(ndim, sizeof(cholmod_sparse*));

	for (i = 0; i < ndim; i++) {
		bases[i] = bsplinebasis(knots[i], nknots[i],
		    coords[i], data->ranges[i], order[i], c);
		if (monodim == i) {
			/* If this dimension is the monotonic one, convert
			 * the basis to T-Splines */
			cholmod_sparse *oldbasis, *tril;

			oldbasis = bases[i];
			tril = cholmod_tril(nsplines[i], c);
			bases[i] = cholmod_l_ssmult(oldbasis, tril, 0, 1, 0, c);
			cholmod_l_free_sparse(&oldbasis, c);
			cholmod_l_free_sparse(&tril, c);
		}

		boxedbases[i] = box(bases[i], bases[i], c);
	}
	if (c->status != CHOLMOD_OK) {
		printf("Basis calculation failed\n");
		return(-1);
	}

	if (verbose)
		printf("Reticulating splines...\n");

	/*
	 * Initialize F and R. 
	 * F = weights
	 * R = weights * data
	 */
	R.rows = F.rows = data->rows;
	R.ndim = F.ndim = data->ndim;
	F.x = malloc(data->rows * sizeof(double));
	F.i = malloc(2*data->ndim * sizeof(int *));
	F.ranges = malloc(2*data->ndim * sizeof(unsigned));
	R.x = malloc(data->rows * sizeof(double));
	R.i = malloc(data->ndim * sizeof(int *));
	R.ranges = malloc(data->ndim * sizeof(unsigned));

	for (i = 0; i < data->ndim; i++) {
		F.i[i] = malloc(data->rows * sizeof(int));
		R.i[i] = malloc(data->rows * sizeof(int));
	}

	memcpy(R.x, weights, data->rows * sizeof(double));
	memcpy(F.x, weights, data->rows * sizeof(double));

	for (i = 0; i < data->ndim; i++) {
		memcpy(R.i[i], data->i[i], data->rows * sizeof(int));
		memcpy(F.i[i], data->i[i], data->rows * sizeof(int));
		R.ranges[i] = F.ranges[i] = data->ranges[i];
	}

	for (i = 0; i < data->rows; i++)
		R.x[i] *= data->x[i];

	/*
	 * Convolve F and R with the basis matrices
	 */

	if (verbose)
		printf("\tConvolving bases...\n");

	for (i = 0; i < data->ndim; i++) {
		if (verbose)
			printf("\t\tConvolving dimension %ld\n",i);

		err=slicemultiply(&F, boxedbases[i], i, c);
		if (err != 0) {
			printf("slicemultiply (F) failed\n");
			return(1);
		}
		err=slicemultiply(&R, bases[i], i, c);
		if (err != 0) {
			printf("slicemultiply (R) failed\n");
			return(1);
		}
	}

	/* Now flatten R into a matrix */

	if (verbose)
		printf("\tFlattening residuals matrix...\n");

	Rmat = flatten_ndarray_to_sparse(&R, sidelen, 1, c);

	for (i = 0; i < R.ndim; i++)
		free(R.i[i]);
	free(R.x); free(R.i); free(R.ranges);

	/* XXX: We reshape, transpose, and then flatten F, which is
	 * potentially memory hungry. This can probably be done in one
	 * step. */

	/*
	 * F now has the number of splines squared as each axis
	 * dimension. We now want to double the dimensionality of F
	 * so that it is n1xn1xn2xn2x... instead of (n1xn1)x(n2xn2)x...
	 */

	if (verbose)
		printf("Transforming fit array...\n");

	/* Fill the ranges array in-place by starting at the back, and
	 * make use of 3/2 = 2/2 = 1 to get the pairs. While here,
	 * rearrange the index columns using the same logic. */
	F.ndim *= 2;
	for (i = F.ndim-1; i >= 0; i--) {
		F.ranges[i] = sqrt(F.ranges[i/2]);
		if (i % 2 == 0)
			F.i[i] = F.i[i/2];
		else
			F.i[i] = malloc(F.rows * sizeof(int));
	}

	/* Now figure out each point's new coordinates */
	for (i = 0; i < F.rows; i++) {
		for (j = 0; j < F.ndim; j += 2) {
			F.i[j+1][i] = F.i[j][i] % F.ranges[j];
			F.i[j][i] = F.i[j][i] / F.ranges[j];
		}
	}

	/* Reorder dimensions of F so that the even-numbered axes come first */
	{
		unsigned int** oldi;
		unsigned* oldranges;

		oldi = F.i;
		oldranges = F.ranges;
		assert(F.ndim>0);
		F.i = malloc(F.ndim * sizeof(int *));
		F.ranges = malloc(F.ndim * sizeof(unsigned));
		for (i = 0; i < F.ndim; i++) {
			if (i % 2 == 0) {
				F.i[i/2] = oldi[i];
				F.ranges[i/2] = oldranges[i];
			} else {
				F.i[F.ndim/2 + i/2] = oldi[i];
				F.ranges[F.ndim/2 + i/2] = oldranges[i];
			}
		}
		free(oldi);
		free(oldranges);
	}

	/* Now flatten F */

	Fmat = flatten_ndarray_to_sparse(&F, sidelen, sidelen, c);
	for (i = 0; i < F.ndim; i++)
		free(F.i[i]);
	free(F.x); free(F.i); free(F.ranges);

	/* XXX: optimization possibilities ended */

	scale2[0] = 1.0;
	fitmat = cholmod_l_add(Fmat, penalty, scale1, scale2, 1, 0, c);
	
	cholmod_l_free_sparse(&Fmat, c);/* we don't need Fmat anymore */

	Rdens = cholmod_l_sparse_to_dense(Rmat, c);
	cholmod_l_free_sparse(&Rmat, c); /* nor Rmat */

	/*
	 * Now, we can solve the linear system 
	 */

	if (verbose)
	    printf("Computing least square solution...\n");

	if (monodim != PHOTOSPLINE_GLAM_NO_MONODIM) {
		coefficients = nnls_normal_block3(fitmat, Rdens,
		    verbose, c);
	} else {
		/* XXX: clamped to one iteration */
		coefficients = cholesky_solve(fitmat, Rdens, c,
		    verbose, 0);
	}

	cholmod_l_free_sparse(&fitmat, c);
	
	/* Clean up detritus */

	if (verbose)
		printf("Done: cleaning up\n");

	cholmod_l_free_sparse(&Fmat, c);
	cholmod_l_free_dense(&Rdens, c);

	for (i = 0; i < ndim; i++) {
		cholmod_l_free_sparse(&bases[i], c);
		cholmod_l_free_sparse(&boxedbases[i], c);
	}
	free(bases);
	free(boxedbases);
	free(nsplines);

	/* Copy out the coefficients */
	if (coefficients == NULL) {
		printf("Solution FAILED\n");
		//out->coefficients = NULL;
		//out->naxes = NULL;
		return(1);
	}

	//out->coefficients = malloc(coefficients->nrow * coefficients->ncol *
	//    sizeof(float));
	for (i = 0; i < coefficients->nrow * coefficients->ncol; i++)
		out_coefficients[i] = ((double *)(coefficients->x))[i];
	//out->naxes = malloc(out->ndim * sizeof(long));
	//for (i = 0; i < out->ndim; i++)
	//	out->naxes[i] = out->nknots[i] - order[i] - 1;

	/* Free our last matrix */
	cholmod_l_free_dense(&coefficients, c);

	/*
	 * If we had a monotonic dimension, it was fit with t-splines. These
	 * must be converted back to b-spline coefficients before returning.
	 *
	 * The equivalent set of b-spline coefficients is the set where
	 * each b-spline coefficient is the sum of all preceding t-spline ones.
	 */
	if (monodim != PHOTOSPLINE_GLAM_NO_MONODIM) {
		long stride1, stride2;
		stride1 = stride2 = 1;

		for (i = 0; i < ndim; i++) {
			if (i < monodim)	
				stride1 *= naxes[i];
			else if (i > monodim)
				stride2 *= naxes[i];
		}

		for (i = 0; i < stride1; i++) {
			for (j = 1; j < naxes[monodim]; j++) {
				for (k = 0; k < stride2; k++) {
				  out_coefficients[i*stride2*naxes[monodim] +
					j*stride2 + k] += out_coefficients[
					i*stride2*naxes[monodim] + (j-1)*stride2 + k];
				}
			}
		}
	}
	return(0);
}

cholmod_sparse*
add_penalty_term(uint64_t* nsplines, double* knots, uint32_t ndim, uint32_t dim, uint32_t order,
   uint32_t porder, double scale, int mono, cholmod_sparse* penalty,
   cholmod_common* c)
{
	cholmod_sparse* penalty_tmp, * penalty_chunk;
	double scale1[2] = {1.0, 0.0}; double scale2[2] = {1.0, 0.0};

	if (scale == 0.0)
		return (penalty);		

	penalty_chunk = calc_penalty(nsplines, knots, ndim, dim, order,
	    porder, mono, c);
	penalty_tmp = penalty;

	/* Add each chunk to the big matrix, scaling by smooth */
	scale2[0] = scale;
	penalty = cholmod_l_add(penalty, penalty_chunk, scale1, scale2,
	    1, 0, c);

	cholmod_l_free_sparse(&penalty_tmp, c);
	cholmod_l_free_sparse(&penalty_chunk, c);
	return (penalty);
}

static cholmod_sparse*
flatten_ndarray_to_sparse(struct ndsparse *array, size_t nrow, size_t ncol,
    cholmod_common* c)
{
	assert(array->ndim>0);
	cholmod_triplet* trip;
	cholmod_sparse* sparse;
	long moduli[array->ndim];
	long i, j, k;

	trip = cholmod_l_allocate_triplet(nrow, ncol, array->rows, 0,
	    CHOLMOD_REAL, c);

	moduli[array->ndim-1] = 1;
	for (i = array->ndim-2; i >= 0; i--)
		moduli[i] = moduli[i+1]*array->ranges[i+1];

	for (i = 0; i < array->rows; i++) {
		k = 0;
		for (j = 0; j < array->ndim; j++)
			k += array->i[j][i]*moduli[j];

		((long *)(trip->j))[i] = k % ncol;
		((long *)(trip->i))[i] = k / ncol;
		((double *)(trip->x))[i] = array->x[i];
	}
	trip->nnz = array->rows;

	sparse = cholmod_l_triplet_to_sparse(trip, trip->nnz, c);
	cholmod_l_free_triplet(&trip, c);

	return (sparse);
}

static void
divided_diffs(int order, int porder, int j, double* knots, double* out)
{
	double a[order], b[order];
	double delta;
	int i;

	/*
	 * Recursively calculate divided differences. Returns a coefficient
	 * array like [1, -2, 1] used to approximate the nth derivative.
	 *
	 * From C. DeBoor, "A Practical Guide to Splines", equation X.16
	 *  2001 edition, page 117
	 */

	/*
	 * Special case order 1 (first derivative by finite differences)
	 * to terminate the recursion.
	 */
	if (porder == 0) {
		out[0] = 1.0;
		return;
	}

	/*
	 * Get each of the (n-1)th derivatives.
	 */

	divided_diffs(order, porder - 1, j + 1, knots, a);
	divided_diffs(order, porder - 1, j, knots, b);

	/*
	 * Get the denominator
	 */

	delta = (knots[j+order+1] - knots[j+porder])/
	    (double)(order - (porder - 1));

	/*
	 * Now subtract them, divide by delta, and return
	 */

	out[0] = -b[0]/delta;
	out[porder] = a[porder-1]/delta;
	for (i = 1; i < porder; i++)
		out[i] = (a[i-1] - b[i])/delta;
}

cholmod_sparse*
calc_penalty(uint64_t* nsplines, double* knots, uint32_t ndim, uint32_t dim, uint32_t order,
    uint32_t porder, int mono, cholmod_common* c)
{
	cholmod_sparse* finitediff, * fd_trans, * DtD, * result;
	cholmod_sparse* tmp, * tmp2;
	cholmod_triplet* trip;
	double divd[porder + 1];
	long i, row, col;

	/* First, we will compute the finite difference matrix,
	 * which looks like this for order 2:
	 * 
	 * 1 -2  1  0 0 ...
	 * 0  1 -2  1 0 ...
	 * 0  0  1 -2 1 ...
	 */

	trip = cholmod_l_allocate_triplet(nsplines[dim] - porder, nsplines[dim],
	    (nsplines[dim] - porder)*(porder+1), 0, CHOLMOD_REAL, c);

	for (row = 0; row < nsplines[dim] - porder; row++) {
		divided_diffs(order, porder, row, knots, divd);

		for (col = row; col < row + porder + 1; col++) {
			((long *)(trip->i))[trip->nnz] = row;
			((long *)(trip->j))[trip->nnz] = col;
			((double *)(trip->x))[trip->nnz] = divd[col-row];

			trip->nnz++;
		}
	}

	finitediff = cholmod_l_triplet_to_sparse(trip, trip->nnz, c);
	cholmod_l_free_triplet(&trip, c);

	if (mono) {
		/* If this dimension is the monotonic one, convert
		 * the basis to T-Splines */
		cholmod_sparse *old, *tril;

		old = finitediff;
		tril = cholmod_tril(nsplines[dim], c);
		finitediff = cholmod_l_ssmult(old, tril, 0, 1, 0, c);
		cholmod_l_free_sparse(&old, c);
		cholmod_l_free_sparse(&tril, c);
	}

	/*
	 * Now we want DtD, which is the transpose of finitediff
	 * multiplied by finitediff
	 */

	fd_trans = cholmod_l_transpose(finitediff, 1, c);
	DtD = cholmod_l_ssmult(fd_trans, finitediff,
	    1 /* DtD is symmetric */, 1, 0, c);
	cholmod_l_free_sparse(&finitediff, c);
	cholmod_l_free_sparse(&fd_trans, c);

	/* Next take kronecker products to form the full P */

	tmp = NULL;
	result = NULL;
	for (i = 0; i < ndim; i++) {
		tmp2 = (i == dim) ? DtD : cholmod_l_speye(
		    nsplines[i], nsplines[i],
		    CHOLMOD_REAL, c);
		tmp2->stype = 1; /* The identity matrix is always symmetric. */

		if (result == NULL) {
			result = tmp2;
			continue;
		}

		tmp = kronecker_product(result, tmp2, c);
		cholmod_l_free_sparse(&result, c);
		cholmod_l_free_sparse(&tmp2, c);

		result = tmp;
	}

	return result;
}

