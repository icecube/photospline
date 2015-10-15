#include <cholmod.h>
#include <math.h>
#include <assert.h>

#include "photospline/detail/splineutil.h"

int
ndsparse_allocate(struct ndsparse* s, size_t rows, size_t ndim)
{
	s->rows = rows;
	s->ndim = ndim;
	s->i = calloc(ndim, sizeof(int*));
	if (!s->i)
		return (1);
	s->ranges = calloc(ndim, sizeof(int));
	if (!s->ranges) {
		free(s->i);
		return (1);
	}
	s->x = calloc(rows, sizeof(double));
	if (!s->x) {
		free(s->i);
		free(s->ranges);
		return (1);
	}
	for (size_t i=0; i<ndim; i++) {
		s->i[i] = calloc(rows, sizeof(int));
		if (!s->i[i]) {
			free(s->x);
			free(s->ranges);
			for (size_t j=0; j<i; j++)
				free(s->i[j]);
			free(s->i);
			return(1);
		}
	}
	return (0);
}

void
ndsparse_free(struct ndsparse* s)
{
	free(s->x);
	free(s->ranges);
	for (size_t i=0; i<s->ndim; i++)
		free(s->i[i]);
	free(s->i);
}

cholmod_sparse *
cholmod_tril(int dim, cholmod_common *c)
{
	cholmod_dense *dtril;
	cholmod_sparse *stril;
	int row, col;

	dtril = cholmod_l_zeros(dim, dim, CHOLMOD_REAL, c);

	/* CHOLMOD dense matrices are in column-major order */
	for (col = 0; col < dim; col++) for (row = col; row < dim; row++)
		((double *)(dtril->x))[col*dim + row] = 1;

	stril = cholmod_l_dense_to_sparse(dtril, 1, c);
	cholmod_l_free_dense(&dtril, c);

	return (stril);
}

static double
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

cholmod_sparse*
bsplinebasis(const double* knots, size_t nknots, const double* x, size_t npts, int order,
    cholmod_common* c)
{
	cholmod_dense* basis;
	cholmod_sparse* sbasis;
	size_t nsplines;
	int row,col,k;

	/*
	 * This matrix has npts rows and nknots-order-1 columns, with order+1
	 * maximum non-zero elements per row (basis functions have limited
	 * extent).
	 *
	 * XXX: We could do this more efficiently by filling a triplet matrix,
	 * but this has little effect on total execution time.
	 */

	nsplines = nknots-order-1;
	basis = cholmod_l_allocate_dense(npts, nsplines, npts, CHOLMOD_REAL, c);

	/* CHOLMOD dense matrices are in column-major order */
	k = 0;
	for (col = 0; col < nsplines; col++) 
		for (row = 0; row < npts; row++, k++) 
			((double *)(basis->x))[k] = bspline(knots, x[row],
			    col, order);

	sbasis = cholmod_l_dense_to_sparse(basis, 1, c);
	cholmod_l_free_dense(&basis, c);

	return (sbasis);
}

/*
 * slicemultiply() multiplies an n-dimensional array by a sparse matrix along
 * the axis specified by dim, returning the result in the original array. It
 * is equivalent to glam.rho() in the Python implementation, and the
 * function 'rho' defined by Eilers and Currie.
 * 
 * If you want to understand how it works, look at the Python implementation.
 */

int
slicemultiply(struct ndsparse *a, cholmod_sparse *b, int dim,
    cholmod_common *c)
{
	cholmod_triplet *section;
	cholmod_sparse *ssection;
	int cols, i, j, k, stride;

	/* Check that the dimensions match */
	if (b->nrow != a->ranges[dim])
		return -1;

	/*
	 * Construct a matrix with a->ranges[dim] rows based on flattening a 
	 * transposed version of a.
	 *
	 * Since we are just shuffling things about, the max number of
	 * non-zero elements is just the number of elements in a.
	 */

	cols = 1;
	for (i = 0; i < a->ndim; i++)
		if (i != dim) cols *= a->ranges[i];

	section = cholmod_l_allocate_triplet(a->ranges[dim], cols, a->rows, 0,
	    CHOLMOD_REAL,c);
	section->nnz = a->rows;

	/* We rotate the array so that dim is at the front, then flatten. */

	assert(a->ndim>0);
	for (i = 0; i < a->rows; i++) {
		/* The flattened row number is straightforward: it is
		 * the dim'th coordinate. */
	
		((long *)(section->i))[i] = a->i[dim][i];

		/* The fastest running index is the dimension we are looking at.
		 * The next fastest is that the one that came after that, in a
		 * modular fashion. For non-primary (dim) indices, build up
		 * an aggregate range step.
		 */

		stride = 1;
		((long *)(section->j))[i] = 0;
		for (k = dim + a->ndim - 1; k > dim; k--) {
			((long *)(section->j))[i] += stride*a->i[k % a->ndim][i];
			stride *= a->ranges[k % a->ndim];
		}
	
		/* Copy the data */
		((double *)(section->x))[i] = a->x[i];
	}

	/* Now obtain the sparse representation */
	ssection = cholmod_l_triplet_to_sparse(section, 0, c);
	if (c->status != CHOLMOD_OK) {
		return(-1);
	}

	cholmod_l_free_triplet(&section, c);

	/* Compute bT . ssection, putting the result in ssection */
	{
		cholmod_sparse *bt, *bta;
		bt = cholmod_l_transpose(b,1,c);
		bta = cholmod_l_ssmult(bt,ssection,0 /*assymmetric */,
		    1 /* compute values */, 0 /* don't bother sorting */,c);

		cholmod_l_free_sparse(&bt,c);
		cholmod_l_free_sparse(&ssection,c);

		ssection = bta;
	}

	/* Now we need to undo the rotation and flattening from above */

	section = cholmod_l_sparse_to_triplet(ssection,c);
	cholmod_l_free_sparse(&ssection,c);

	/* Set up the nd-array again, bearing in mind that it need not have
	 * the same number of non-zero elements as before, and that the range
	 * along dimension dim is also now the number of columns in b */

	for (i = 0; i < a->ndim; i++)
		a->i[i] = realloc(a->i[i],sizeof(int)*section->nnz);
	a->x = realloc(a->x,sizeof(double)*section->nnz);
	a->rows = section->nnz;
	a->ranges[dim] = b->ncol;

	/* We unflatten the array and rotate so that the front is at dim */

	for (i = 0; i < a->rows; i++) {
		/* The flattened column number is straightforward: it is
		 * the dim'th coordinate. */
	
		a->i[dim][i] = ((long *)(section->i))[i];

		/* The fastest running index is the dimension we are looking at.
		 * The next fastest is that the one that came after that, in a
		 * modular fashion. For non-primary (dim) indices, build up
		 * an aggregate range step.
		 */

		stride = 1;
		for (k = dim + a->ndim - 1; k > dim; k--)
			stride *= a->ranges[k % a->ndim];
		j = ((long *)(section->j))[i];
		for (k = dim+1; k < dim + a->ndim; k++) {
			/* See how many of this dimension's strides we fit */
			stride /= a->ranges[k % a->ndim];
			a->i[k % a->ndim][i] = j/stride;
			/* Collect the remainder */
			j = j % stride;
		}

		/* Copy the data */
		a->x[i] = ((double *)(section->x))[i];
	}

	cholmod_l_free_triplet(&section,c);

	return 0;
}

/* Computes the kronecker product of sparse matrices a and b */

cholmod_sparse *
kronecker_product(cholmod_sparse *a, cholmod_sparse *b, cholmod_common *c)
{
	cholmod_sparse *final;
	cholmod_triplet *ta, *tb, *tf;
	int i, j;

	ta = cholmod_l_sparse_to_triplet(a,c);
	tb = cholmod_l_sparse_to_triplet(b,c);
	tf = cholmod_l_allocate_triplet(ta->nrow*tb->nrow, ta->ncol*tb->ncol,
	    ta->nnz*tb->nnz, (ta->stype == tb->stype) ? ta->stype : 0,
	    CHOLMOD_REAL, c);
	tf->nnz = ta->nnz*tb->nnz;

	for (i = 0; i < ta->nnz; i++) {
		for (j = 0; j < tb->nnz; j++) {
			/* New data is x_a*x_b */
			((double *)(tf->x))[i*tb->nnz + j] = ((double *)(ta->x))[i]
			    * ((double *)(tb->x))[j];
			/* New column index is col_a*ncol_b + col_b */
			((long *)(tf->j))[i*tb->nnz + j] = ((long *)(ta->j))[i]*tb->ncol
			    + ((long *)(tb->j))[j];
			/* New row index is row_a*nrow_b + row_b */
			((long *)(tf->i))[i*tb->nnz + j] = ((long *)(ta->i))[i]*tb->nrow
			    + ((long *)(tb->i))[j];
		}
	}
		
	final = cholmod_l_triplet_to_sparse(tf, 0, c);
	
	cholmod_l_free_triplet(&ta,c);
	cholmod_l_free_triplet(&tb,c);
	cholmod_l_free_triplet(&tf,c);

	return final;
}

/*
 * box is implemented as follows in python:
 * def box(A,B):
 * 	ea = numpy.ones((1,A.shape[1]),float)
 * 	eb = numpy.ones((1,B.shape[1]),float)
 * 	return numpy.matrix(numpy.asarray(numpy.kron(A, eb)) * \
 *          numpy.asarray(numpy.kron(ea, B)))
 *
 * What it does is form a new matrix with the same number of rows as A and B,
 * where each element is the row is the product of an element of A and the
 * corresponding row of B. For example:
 *
 * [ A_00*B_00 A_00*B_01 A_00*B_02 A_01*B_00 A_01*B_01 ... ]
 * [ A_10*B_10 A_10*B_11 ... ]
 */


cholmod_sparse *
box(cholmod_sparse *a, cholmod_sparse *b, cholmod_common *c)
{
	cholmod_sparse *final;
	cholmod_triplet *ta, *tb, *tf;
	int i, rownum;

	struct rowitem {
		long col;
		double val;

		struct rowitem *next;
	};

	struct rowitem **brows, **curbrows;
	struct rowitem *item;

	/* We need the number of rows in both matrices to be equal */
	if (a->nrow != b->nrow)
		return NULL;
	
	ta = cholmod_l_sparse_to_triplet(a,c);
	tb = cholmod_l_sparse_to_triplet(b,c);

	/*
	 * The number of non-zero entries cannot be larger than the number
	 * of nonzero entries in a times the length of a row in b
	 */
	tf = cholmod_l_allocate_triplet(ta->nrow, ta->ncol*tb->ncol,
	    ta->nnz*tb->ncol, 0, CHOLMOD_REAL, c);
	tf->nnz = 0;

	/* Decompose b into its (sparse) rows */
	brows = calloc(tb->nrow,sizeof(struct rowitem *));
	curbrows = calloc(tb->nrow,sizeof(struct rowitem *)); /* fast lookup for list end*/

	for (i = 0; i < tb->nnz; i++) {
		rownum = ((long *)(tb->i))[i];
		item = malloc(sizeof(struct rowitem));
		item->next = NULL;
		item->col = ((long *)(tb->j))[i];
		item->val = ((double *)(tb->x))[i];

		if (brows[rownum] == NULL)
			brows[rownum] = item;
		else
			curbrows[rownum]->next = item;
		
		curbrows[rownum] = item;
	}

	/* Free curbrows now that we are done with it */
	free(curbrows);

	for (i = 0; i < ta->nnz; i++) {
		rownum = ((long *)(ta->i))[i];

		for (item = brows[rownum]; item != NULL; item = item->next) {
			((long *)(tf->i))[tf->nnz] = rownum;
			((long *)(tf->j))[tf->nnz] = ((long *)(ta->j))[i]*tb->ncol +
			    item->col;
			((double *)(tf->x))[tf->nnz] = ((double *)(ta->x))[i] * item->val;

			tf->nnz++;
		}
	}

	/* Free all the row items */
	for (i = 0; i < tb->nrow; i++) {
		struct rowitem *tmp;
		item = brows[i];

		while (item != NULL) {
			tmp = item;
			item = item->next;
			free(tmp);
		}
	}
	free(brows);

	final = cholmod_l_triplet_to_sparse(tf, 0, c);

	/* clean up */
	cholmod_l_free_triplet(&ta, c);
	cholmod_l_free_triplet(&tb, c);
	cholmod_l_free_triplet(&tf, c);

	return (final);
}
		

void
print_sparse(cholmod_sparse *a, cholmod_common *c)
{
	cholmod_dense *dense;
	int i,j;
	double val;

	dense = cholmod_l_sparse_to_dense(a,c);
	printf("----\n");
	for (i = 0; i < dense->nrow; i++) {
		for (j = 0; j < dense->ncol; j++) {
			val = ((double *)(dense->x))[j*dense->nrow + i];
			if (fabs(val) <= 1e-12) {
				printf("    -    ");
			} else {
				printf(" %- .1e",val);
			}
		}
		printf("\n");
	}
	printf("----\n");

	cholmod_l_free_dense(&dense,c);
}

void
print_factor(cholmod_factor *L, cholmod_common *c)
{
	cholmod_sparse *Lsparse;
	cholmod_factor *Lcopy;
	cholmod_dense *Ldense;
	int i, j;
	double val;
	Lcopy = cholmod_l_copy_factor(L, c);
	Lsparse = cholmod_l_factor_to_sparse(Lcopy, c);
	Ldense = cholmod_l_sparse_to_dense(Lsparse, c);

	printf("----\n");
	for (i = 0; i < Ldense->nrow; i++) {
		for (j = 0; j < i+1; j++) {
			val = ((double*)(Ldense->x))[j*(Ldense->nrow)+i];
			if (i >= j) {
				if (val == 0) printf("    -    ");
				else if (val == 1) printf("    %d    ",
				    (int)val);
				else printf(" % -.1e",val);
			} else {
				printf("         ");
			}
		}	
		printf("\n");
	}
	printf("----\n");

	cholmod_l_free_factor(&Lcopy, c);
	cholmod_l_free_dense(&Ldense, c);
	cholmod_l_free_sparse(&Lsparse, c);
}

