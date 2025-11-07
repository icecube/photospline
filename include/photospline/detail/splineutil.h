
#ifndef PHOTOSPLINE_CFITTER_SPLINEUTIL_H
#define PHOTOSPLINE_CFITTER_SPLINEUTIL_H

#include <cholmod.h>

#ifdef __cplusplus
extern "C" {
#endif

struct ndsparse {
	/*
	 * This is an ntuple, and is similar to the CHOLMOD triplet
	 * type.
	 */

	size_t rows;
	size_t ndim;
	
	/* i contains coordinate indices, x the data values */
	double* x;
	unsigned int** i;
	unsigned int* ranges; /* Index ranges for each column */
};

/* returns 0 on success, nonzero otherwise */
int ndsparse_allocate(struct ndsparse* s, size_t rows, size_t ndim);
void ndsparse_free(struct ndsparse* s);

cholmod_sparse *bsplinebasis(const double *knots, size_t nknots, const double *x,
    size_t npts, int order, cholmod_common *c);

int slicemultiply(struct ndsparse *a, cholmod_sparse *b, int dim,
    cholmod_common *c);
cholmod_sparse *kronecker_product(cholmod_sparse *a, cholmod_sparse *b,
    cholmod_common *c);
cholmod_sparse *box(cholmod_sparse *a, cholmod_sparse *b, cholmod_common *c);
cholmod_sparse *cholmod_tril(int dim, cholmod_common *c);

cholmod_dense *nnls_lawson_hanson(cholmod_sparse *A, cholmod_dense *y,
    double tolerance, int min_iterations, int max_iterations, unsigned npos,
    int normaleq, int verbose, cholmod_common *c);
cholmod_dense *nnls_normal_block(cholmod_sparse *AtA, cholmod_dense *Atb,
   int verbose, cholmod_common *c);

cholmod_dense *nnls_normal_block_updown(cholmod_sparse *AtA, cholmod_dense *Atb,
   int verbose, cholmod_common *c);

cholmod_dense *
nnls_normal_block3(cholmod_sparse *AtA, cholmod_dense *Atb, int verbose,
   cholmod_common *c);

void print_sparse(cholmod_sparse *a, cholmod_common *c);

void print_factor(cholmod_factor *L, cholmod_common *c);

#ifdef __cplusplus
} //extern "C"
#endif

#endif /* PHOTOSPLINE_CFITTER_SPLINEUTIL_H */

