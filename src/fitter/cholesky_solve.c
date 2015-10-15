/* Use non-portable extensions */
#define _GNU_SOURCE 1

#include <stdio.h>
#include <unistd.h>
#include <math.h>

#ifdef __FreeBSD__
#include <sys/param.h>
#include <sys/cpuset.h>

typedef cpuset_t cpu_set_t;
#endif

#include <time.h>
#include <float.h>
#include <stdbool.h>
#include <assert.h>
#include <string.h>
#include <pthread.h>
#include <sched.h>

#include <cholmod.h>
#include "photospline/detail/splineutil.h"
#include "cholesky_solve.h"

static int
intcmp(const void *xa, const void *xb);

cholmod_dense *
cholesky_solve(cholmod_sparse *AtA, cholmod_dense *Atb, cholmod_common *c,
     int verbose, int n_resolves)
{
	int i,j,nvar;
	double ones[2] = {1., 0}, mones[2] = {-1., 0};
	double sum_delta_b;
	cholmod_dense *x, *delta_x, *delta_Atb;
	cholmod_factor *L;
	clock_t t0, t1;

	/* set symmetry */
	AtA->stype = 1;

	/* allocate workspace */
	nvar = AtA->nrow;
	if (n_resolves > 0) {
		delta_Atb = cholmod_l_allocate_dense(nvar, 1, nvar,
                            CHOLMOD_REAL, c);
		delta_x = cholmod_l_allocate_dense(nvar, 1, nvar,
                            CHOLMOD_REAL, c);
	}
	
	t0 = clock();

	/* compute the nonzero pattern of the Cholesky factorization */	
	L = cholmod_l_analyze(AtA, c);

	if (verbose) {
		t1 = clock();
		printf("Analyze[%d]: %f s\n", nvar,
		    (double)(t1-t0)/(CLOCKS_PER_SEC));
		t0 = clock();
	}

	/* compute the numerical values of the Cholesky factorization */	
	cholmod_l_factorize(AtA, L, c);

	if (verbose) {
		t1 = clock();
		printf("Factorize[%d]: %f s\n",nvar,(double)(t1-t0)/(CLOCKS_PER_SEC));
		t0 = clock();
	}

	/* solve the system AtA*x = Atb */
	x = cholmod_l_solve(CHOLMOD_A, L, Atb, c);

	if (verbose) {
		t1 = clock();
		printf("Solve[%d]: %f s\n",nvar,(double)(t1-t0)/(CLOCKS_PER_SEC));
		t0 = clock();
	}

	/* refine the solution if needed */
	for (j = 0; j < n_resolves; j++) {
		sum_delta_b = 0.;

		if (verbose)
			t0 = clock();

		/* copy the right-hand side */
		for (i = 0; i < nvar; i++)
		    ((double *)(delta_Atb->x))[i] =
		        ((double *)(Atb->x))[i];

		/* delta_b = A*(x + delta_x) - b */
		cholmod_l_sdmult(AtA, 0, ones, mones, x, delta_Atb, c);

		/* solve A*delta_x = delta_b */
		delta_x = cholmod_l_solve(CHOLMOD_A, L, delta_Atb, c);

		/* x' = x - delta_x */
		for (i = 0; i < nvar; i++) {
			sum_delta_b += (((double *)(delta_Atb->x))[i])*(((double *)(delta_Atb->x))[i]);
			((double *)(x->x))[i] -= ((double *)(delta_x->x))[i];
		}

		if (verbose) {
			t1 = clock();
			printf("reSolve %d: %f s (db: %e)\n",
		    	    j,(double)(t1-t0)/(CLOCKS_PER_SEC),sqrt(sum_delta_b));
		}

		cholmod_l_free_dense(&delta_x, c);
	}

	cholmod_l_free_factor(&L, c);
	if (n_resolves > 0) cholmod_l_free_dense(&delta_Atb, c);

	return (x);
}

cholmod_sparse* get_column(cholmod_sparse *A, long k, 
    long *iPerm, long *Fset, long nF, cholmod_common *c)
{
	cholmod_sparse *R;
	long *Ap, *Ai, *Anz, *Rp, *Ri;
	double *Ax, *Rx;
	long row;
	int i, j, nz, A_col_nz;

	/* n-by-1, rows not sorted, packed, stype 0 (unsymmetric) */
	R = cholmod_l_allocate_sparse(A->nrow, 1, nF, 
	    false, true, 0, CHOLMOD_REAL, c);

	Ap = (long*)(A->p);   /* pointers to the start of each column */
	Ai = (long*)(A->i);   /* row indices */
	Ax = (double*)(A->x); /* numerical values */
	Anz = (long*)(A->nz); /* length of each column */

	/* get the number of non-zero entries
	 * in column k
	 */
	if (A->packed) A_col_nz = Ap[k+1]-Ap[k];
	else A_col_nz = Anz[k];

	Rp = (long*)(R->p);
	Ri = (long*)(R->i);
	Rx = (double*)(R->x);

	nz = 0;
	/* copy rows from A[F,k] to R */
	for (i = 0; i < nF; i++) {
	/* XXX: O(nF*A_col_nz) brute force search */
		for (j = 0; j < A_col_nz; j++) {
			if (Ai[Ap[k]+j] == Fset[i]) {
				Ri[nz] = Ai[Ap[k]+j];
				Rx[nz] = Ax[Ap[k]+j];
				nz++;	
			}
		} 
	}

	/* symbolically permute the order of the rows */
	if (iPerm != NULL) {
		for (i = 0; i < nz; i++) {
			row = Ri[i];
			Ri[i] = iPerm[row];
		}
	}

	/* R is packed */
	Rp[0] = 0;
	Rp[1] = nz;

	return(R);
}

static int
intcmp(const void *xa, const void *xb)
{
	const int *a, *b;
	a = xa; b = xb;

	if (*a < *b)
		return (-1);
	else if (*a > *b)
		return (1);
	else
		return (0);
}

cholmod_factor* 
modify_factor(cholmod_sparse *A, cholmod_factor *L,
    long *F, long *nF_, long *G, long *nG_, long *H1, long *nH1_,
    long *H2, long *nH2_, int verbose, cholmod_common *c)
{
	double flop_ratio;
	bool update;
	long nF, nG, nH1, nH2;
	int n_threads;

	nF = *nF_;
	nG = *nG_;
	nH1 = *nH1_;
	nH2 = *nH2_;

	update = false;
	flop_ratio = 0.0;
	/*
	 * If available, use the flop counts from the most recent call to
	 * cholmod_analyze or cholmod_rowadd/rowdel to estimate the expense
	 * of the job.
	 *
	 * XXX: assumptions here are:
	 * 1. we can get ~ 16x speedup from 16 threads in GotoBLAS. 
	 * 2. the size of F is roughly constant: the next factorization will
	 *    be about as expensive as the previous one.
	 */

	/* 
	 * The operations performed inside of GotoBLAS are ~ 9 times faster
	 * than the plain-C operations in rowadd/rowdel, and we're running
	 * with 16 threads, so we expect functions taking the same number of
	 * flops to be ~144 times faster when using cholmod_factorize than
	 * rowadd/rowdel. Twiddle as appropriate.
	 */
	#define GOTO_SPEEDUP 9.0
	n_threads = get_nthreads();

	if ((c->modfl <= 0) && (c->lnz > 0)) {
		/* 
		 * Rowadd/rowdel are inherently single-threaded and use scalar
		 * arithmetic; switching modes too early in the problem can be
		 * disastrously slow. Here, we compute a conservative upper 
		 * bound on the work required for rowadd/rowdel, ~ O(|L|). 
		 */ 
		c->modfl = c->lnz;
	}

	if (L == NULL || nF == 0) {
		/* First iteration, always do a full factorization */
		update = false;
	} else if ((c->fl > 0) && (c->modfl > 0)) {
		/*
		 * Estimate work based on previous iterations.
		 */
		flop_ratio = (c->fl) / ( GOTO_SPEEDUP * n_threads * 
		    ((double)(nH1 + nH2))*(c->modfl) );
		update = (flop_ratio > 1.0);
	}

	#undef GOTO_SPEEDUP

	if (verbose)
		printf("\tFactor work: %.0lf Mod work: %.0lf\n",
		    c->fl, c->modfl);

	/* short-circuit for rank-1 updates */
	if (nH1 + nH2 == 1) update = true;

	if (!update && verbose)
		printf("\tRecomputing factorization from scratch "
		    "(F[%ld], G[%ld], H1[%ld], H2[%ld], ratio=%.2lf)\n",
		    nF, nG, nH1, nH2, flop_ratio);

	return(modify_factor_p(A, L, F, nF_, G, nG_, H1, nH1_, H2, nH2_,
	    update, verbose, c));

}

cholmod_factor* 
modify_factor_p(cholmod_sparse *A, cholmod_factor *L,
    long *F, long *nF_, long *G, long *nG_, long *H1, long *nH1_,
    long *H2, long *nH2_, bool update, bool verbose, cholmod_common *c)
{
	cholmod_sparse *col;
	int i, j, k, update_ready;
	long *iPerm;
	long nF, nG, nH1, nH2;
	clock_t t0, t1;
	
	nF = *nF_;
	nG = *nG_;
	nH1 = *nH1_;
	nH2 = *nH2_;

	iPerm = NULL;

	update_ready = (L && (L->n == (nF + nG)));

	/* Compute the inverse of the fill-reducing permutation */
	if (L && L->Perm) {
		iPerm = (long*)malloc(sizeof(long)*L->n);
		for (i = 0; i < L->n; i++)
			iPerm[((long*)(L->Perm))[i]] = i;
	}

	t0 = clock();

	/*
	 * Remove elements in H1 from F, and add them to G,
	 * exploiting the fact that H1 elements are in order
	 * relative to their order in F.
	 */
	for (i = 0, j = 0; i < nH1; i++) {
		G[nG++] = H1[i];
		/*
		 * NB: the following is _not_ an infinite loop, as every element 
		 * of H1 is guaranteed to also be in F and both are monotonic.
		 */
		while (F[j] != H1[i])
			j++;
		for (k = j+i; k+1 < nF; k++)
			F[k-i] = F[k-i+1];

		if (update && update_ready) {
			/* remove the row from the factorization */
			/* XXX: pass non-zero pattern of row, instead of NULL */
			cholmod_l_rowdel((iPerm != NULL) ? iPerm[H1[i]] : H1[i],
			    NULL, L, c);
		}
	}
	if (verbose && (update && update_ready) && (nH1 > 0)) {
		t1 = clock();
		printf("\tDelete %ld rows: %.2f s\n",nH1,
		    (double)(t1-t0)/(CLOCKS_PER_SEC));
	}
	nF -= nH1;
	nH1 = 0;

	t0 = clock();

	/* And vice versa */
	for (i = 0, j = 0; i < nH2; i++) {
		F[nF++] = H2[i];
		while (G[j] != H2[i])
			j++;
		for (k = j+i; k+1 < nG; k++)
			G[k-i] = G[k-i+1];

		if (update && update_ready) {
			/* 
			 * Extract column H2[i] from A, zeroing any
			 * row not in F and permuting the rows to
			 * match the ordering of the columns in L.
			 */
			col = get_column(A, H2[i], iPerm, F, nF, c);
			/* insert the column at the permuted index */
			cholmod_l_rowadd((iPerm != NULL) ? iPerm[H2[i]] : H2[i], 
			    col, L, c);
			cholmod_l_free_sparse(&col, c);
		}
	}
	if (verbose && (update && update_ready) && (nH2 > 0)) {
		t1 = clock();
		printf("\tAdd %ld rows: %.2f s\n",nH2,
		    (double)(t1-t0)/(CLOCKS_PER_SEC));
	}
	nG -= nH2;
	nH2 = 0;

	qsort(G, nG, sizeof(G[0]), intcmp);
	qsort(F, nF, sizeof(F[0]), intcmp);

	/* 
	 * If an update was requested, but L was only computed for a submatrix,
	 * then we have to compute it from scratch once in simplical form. This
	 * is a huge waste for single-row updates, but will allow use of
	 * rowadd/rowdel in the next iteration.
	 */
	if (update && !update_ready) {
		t0 = clock();
		cholmod_l_free_factor(&L, c);
		L = cholmod_l_analyze(A, c);

		if (iPerm != NULL)
			free(iPerm);

		iPerm = (long *)malloc(sizeof(long)*L->n);
		for (i = 0; i < L->n; i++)
			iPerm[((long*)(L->Perm))[i]] = i;

		L = recompute_factor(A, L, iPerm, F, nF, c);
		if (verbose) {
			t1 = clock();
			printf("\tFactorize[%ld]: %.2f s\n",nF,
			    (double)(t1-t0)/(CLOCKS_PER_SEC));
		}
	}

	/*
	 * If no single-row updates were requested, recompute the
	 * factorization of A[:,F][F,:] from scratch by the fastest and least
	 * memory-hungy method available.
	 */
	if (!update) {
		cholmod_sparse *A_F;

		t0 = clock();
		cholmod_l_free_factor(&L, c);
		/* XXX: A must secretly be symmetric */
		A->stype = 0;
		A_F = cholmod_l_submatrix(A, F, nF, F, nF, 1, 1, c);
		A->stype = 1;
		A_F->stype = 1;
		L = cholmod_l_analyze(A_F, c);
		cholmod_l_factorize(A_F, L, c);
		cholmod_l_free_sparse(&A_F, c);

		if (verbose) {
			t1 = clock();
			printf("\tFactorize[%ld]: %.2f s\n",nF,
			    (double)(t1-t0)/(CLOCKS_PER_SEC));
		}
	}

	*nF_  = nF;
	*nG_  = nG;
	*nH1_ = nH1;
	*nH2_ = nH2;

	if (iPerm != NULL)
		free(iPerm);
	
	return(L);	
}

cholmod_factor * 
recompute_factor(cholmod_sparse *A, cholmod_factor *L, long *iPerm,
    long *F, unsigned long nF, cholmod_common *c)
{
	long *iF, *FPerm;
	unsigned long *Lrows;
	long *LPerm, *LColCount, *L_FColCount;
	long *Li, *Lp, *Lnz, *Lnext, *L_Fi, *L_Fp, *L_Fnz;
	double *Lx, *L_Fx;
	cholmod_sparse *A_F;
	cholmod_factor *L_F;
	unsigned long i, j, nz, nFPerm;
	int common_nmethods, Astype;
	bool common_postorder;

	/* 
	 * Scortched earth: free the unneeded numeric 
	 * bits of L, should they exist.
	 */
	cholmod_l_change_factor(CHOLMOD_PATTERN, 
	    false, /* to LL' */
	    false, /* to supernodal */
	    false, /* to packed */
	    false, /* to monotonic */
	    L, c);

	LPerm = (long*)(L->Perm);
	LColCount = (long*)(L->ColCount);

	/* build an inverse mapping for F */
	iF = (long*)malloc(sizeof(long)*L->n);
	for (i = 0; i < L->n; i++) iF[i] = -1;
	for (i = 0; i < nF; i++) iF[F[i]] = i;

	/* Permute F to match L->Perm */
	nFPerm = 0;
	assert(nF>0);
	FPerm = (long*)malloc(sizeof(long)*nF);
	Lrows = (unsigned long*)malloc(sizeof(unsigned long)*nF);
	if (iPerm != NULL) {
		/* calculate the permuation as it applies to subset F */
		for (i = 0; i < L->n; i++) {
			if (iF[LPerm[i]] >= 0) {
				FPerm[nFPerm] = iF[LPerm[i]];
				/* 
				 * This row of L_F (corrensponding to F[j]) 
				 * corresponds to the position of F[j] in the
				 * fill-reducing permuation.
				 */
				Lrows[nFPerm] = iPerm[F[FPerm[nFPerm]]];
				nFPerm++;
			}
		}
		
	} else {
		/* mostly a no-op in natural ordering */
		for (i = 0; i < nF; i++) {
			FPerm[nFPerm] = i;
			Lrows[nFPerm] = F[i];
			nFPerm++;
		}
	}
	assert(nFPerm==nF);

	Astype = A->stype;
	assert( Astype != 0);
	A->stype = 0; /* cholmod_submatrix doesn't like symmetric matrices */
	A_F = cholmod_l_submatrix(A, F, nF, F, nF, 1, 1, c);
	A->stype = Astype;
	A_F->stype = Astype;
	A_F->stype = 1;

	/* 
	 * Since we intend to copy the results of the submatrix factorization
	 * in to L, we need to force cholmod_analyze to use exactly the same
	 * permutation previously calculated for L, without postordering.
	 */
	common_nmethods = c->nmethods;
	common_postorder = c->postorder;
	c->nmethods = 1;
	c->postorder = false;

	/* calculate pattern of L in the given permuation */
	L_F = cholmod_l_analyze_p(A_F, FPerm, NULL, -1, c);
	L_FColCount = (long*)(L_F->ColCount);

	/* restore alternate orderings */
	c->nmethods = common_nmethods;
	c->postorder = common_postorder;

#ifndef NDEBUG
	/* sanity check: did we get the permutation we asked for? */
	long* L_FPerm = (long*)(L_F->Perm);
	for (i = 0; i < nF; i++) assert( FPerm[i] == L_FPerm[i] );
#endif

	cholmod_l_factorize(A_F, L_F, c);

	cholmod_l_free_sparse(&A_F, c);

	/* We really, really, really don't want L in supernodal form. */
	cholmod_l_change_factor(CHOLMOD_REAL, 
	    false, /* to LL' */
	    false, /* to supernodal */
	    false, /* to packed */
	    false, /* to monotonic */
	    L_F, c);

	/* 
	 * XXX: swizzle the nonzero pattern of L to match
	 * L'. Does this break anything down the line?
	 */
	for (i = 0; i < L->n; i++) LColCount[i] = 1;
	for (i = 0; i < nF; i++) LColCount[Lrows[i]] = L_FColCount[i];

	/* 
	 * Restore numeric bits of L, initialized to identity 
	 * and with the new nonzero pattern.
	 */
	cholmod_l_change_factor(CHOLMOD_REAL, 
	    false, /* to LL' */
	    false, /* to supernodal */
	    false, /* to packed */
	    false, /* to monotonic */
	    L, c);

	/* 
	 * Now that we have simplicial numeric factorizations
	 * of both L_F and L (currently identity), we can get 
	 * the appropriate pointers to the numerical bits of both.
	 */
	Lp =    (long*)(L->p);
	Lnz =   (long*)(L->nz);
	Li =    (long*)(L->i);
	Lnext = (long*)(L->next);
	Lx =    (double*)(L->x);
	L_Fp =  (long*)(L_F->p);
	L_Fnz = (long*)(L_F->nz);
	L_Fi =  (long*)(L_F->i);
	L_Fx =  (double*)(L_F->x);

	/* 
	 * Copy the numeric factorization from L_F to 
	 * the corresponding columns and rows of L.
	 */
	for (i = 0; i < nF; i++) {
		/* 
		 * NB: cholmod_factor is never explicitly packed
		 * (unlike cholmod_sparse, L->nz is always meaningful).
		 */
		if (L_Fnz) nz = L_Fnz[i];
		else nz = L_Fp[i+1] - L_Fp[i];

		if ( Lp[Lnext[Lrows[i]]] - Lp[Lrows[i]] < nz ) {
			cholmod_l_reallocate_column(Lrows[i], nz, L, c);
#if 0
			printf("L->nz[%ld] <= %ld, L_F->nz[%d] = %ld\n", 
		    	    Lrows[i], Lp[Lnext[Lrows[i]]] - Lp[Lrows[i]],
			    i, L_Fnz[i]);
#endif
		}

		assert( Lp[Lnext[Lrows[i]]] - Lp[Lrows[i]] >= nz );

		for (j = 0; j < nz; j++) {
			/*
			 * Adjust index of each row to match
			 * its position in L.
			 */
			Li[Lp[Lrows[i]]+j] = Lrows[L_Fi[L_Fp[i]+j]] ;	
			/*
			 * Numerical values can just be copied
			 * to the appropriate column.
			 */
			Lx[Lp[Lrows[i]]+j] = L_Fx[L_Fp[i]+j];	
		}
		Lnz[Lrows[i]] = nz;	
	}

	/* we're done with L_F */
	cholmod_l_free_factor(&L_F, c);
	free(iF);
	free(FPerm);
	free(Lrows);

	return(L);
}

cholmod_sparse *
submatrix_symm(cholmod_sparse *A, long *rows, unsigned long nrows,
    long *cols, unsigned long ncols, cholmod_common *c)
{
	cholmod_sparse *C;
	int i, j, pend, pstart, last;
	long *Ap, *Ai, *Anz, *Cp, *Ci;
	double *Ax, *Cx;
	long *cnz, *irows, *icols;
	int *ccol_start, *ccol_stop;
	bool csorted;
	long nz;

	assert( A->stype != 0 );
	assert( A->nrow == A->ncol );
	assert( A->sorted );
	assert( ncols <= A->ncol );

	Ap = (long*)(A->p);
	Ai = (long*)(A->i);
	Anz = (long*)(A->nz);
	Ax  = (double*)(A->x);

	irows = (long*)malloc(sizeof(long)*A->nrow);
	icols = (long*)malloc(sizeof(long)*A->ncol);
	cnz   = (long*)malloc(sizeof(long)*ncols);
	ccol_start = (int*)malloc(sizeof(int)*ncols);
	ccol_stop  = (int*)malloc(sizeof(int)*ncols);

	for (i = 0; i < A->nrow; i++) irows[i] = icols[i] = -1;
	for (i = 0; i < nrows; i++) irows[rows[i]] = i;
	for (i = 0; i < ncols; i++) {
		icols[cols[i]] = i;
		ccol_start[i] = ccol_stop[i] = -1;
		cnz[i] = 0;
	}

	nz = 0;
	csorted = true;
	for (i = 0; i < A->ncol; i++) {
		/* Redundant with the assert above, but keeps
		 * clang-analyzer happy */
		assert(i < A->nrow);
		if (icols[i] < 0) continue;

		pend = (A->packed) ? Ap[i+1] : Anz[i];
		pstart = Ap[i];
		j = 0;

		/* find the diagonal entry in this column */
		while ((Ai[pstart+j] < i) && (j < pend)) j++;

		/* there's nothing below the diagonal */
		if ((A->stype < 0) && (Ai[Ap[i]+j] < i)) continue;

		/* set up bounds for the given symmetry */
		pstart = (A->stype < 0) ? j : 0;
		pend = (A->stype < 0) ? pend-Ap[i] : j+1;

		j = pstart;
		pstart = -1;
		last = pstart;
		for ( ; j < pend; j++) {
			/* select rows in the set */
			if (irows[Ai[Ap[i]+j]] < 0) continue;
			if (pstart < 0) pstart = j;
			nz += 1;
			cnz[icols[i]] += 1;
			/* check for out-of-order rows */	
			if (j < last) csorted = false;
			last = j;
		}
		pend = last+1;

		/* set up the bounds for the next traversal */
		if (cnz[icols[i]] > 0) {
			ccol_start[icols[i]] = Ap[i]+pstart;
			ccol_stop[icols[i]] = Ap[i]+pend;
		}
		
	}

	/* nrows-by-ncols, rows not sorted, packed, stype of input */
	C = cholmod_l_allocate_sparse(nrows, ncols, nz, 
	    csorted, true, A->stype, CHOLMOD_REAL, c);

	Cp = (long*)(C->p);
	Ci = (long*)(C->i);
	Cx  = (double*)(C->x);

	Cp[0] = 0;
	int p;
	for (i = 0; i < ncols; i++) {
		if (cnz[i] == 0) {
			Cp[i+1] = Cp[i];
			continue;
		}	
		p = Cp[i];
		for (j = ccol_start[i]; j < ccol_stop[i]; j++) {
			/* is this row one we want? */
			if (irows[Ai[j]] < 0) continue;
			/* set row indices and copy values */
			Ci[p] = irows[Ai[j]];
			Cx[p] = Ax[j];
			p++;
		}
		Cp[i+1] = p;
	}

	free(irows);
	free(icols);
	free(cnz);
	free(ccol_start);
	free(ccol_stop);

	if (!csorted)
		cholmod_l_sort(C, c);

	return(C);
}

/*
 * Calculate the 2-norm ||Ax - b|| to within an additive constant (b'b)
 *
 * ||Ax - b|| - b'b = x'A'Ax - 2x'A'b = x'(A'Ax - 2A'b)
 */

double
calc_residual(cholmod_sparse *AtA, cholmod_dense *Atb, cholmod_dense *x,
    cholmod_common *c)
{
	int i;
	double result;
	cholmod_dense *AtAx;
	int nvar = x->nrow;
	double ones[2]  = { 1., 0};
	double mtwos[2] = {-2., 0};

	/* AtAx = 0.5*A'Ax - A'b */
	AtAx = cholmod_l_copy_dense(Atb, c);
	cholmod_l_sdmult(AtA, 0, ones, mtwos, x, AtAx, c);
	
	/* XXX: hopefully vector optimizations will pick this up */
	result = 0;
	for (i = 0; i < nvar; i++)
		result += ((double*)(x->x))[i] * ((double*)(AtAx->x))[i];

	cholmod_l_free_dense(&AtAx, c);
	
	return(result);
}

void
evaluate_descent(void *trial_)
{
	descent_trial *trial = (descent_trial*)trial_;
	int i, err;
	double *xptr;
	cholmod_dense *x, *x_F;
	const long *F;
	long nF;

	/*
	 * Bind to the CPU corresponding to the thread id, unless we're on
	 * OS X, where the thread affinity API is...special.
	 */
#ifndef __APPLE__
	cpu_set_t cpus;
	CPU_ZERO(&cpus);
	CPU_SET(trial->id, &cpus);
#ifdef __FreeBSD__
	err = cpuset_setaffinity(CPU_LEVEL_WHICH, CPU_WHICH_TID, -1,
	    sizeof(cpus), &cpus);
#else
	err = sched_setaffinity(0, sizeof(cpus), &cpus);
#endif
#else
	err = 0;
#endif
	
	if (err)
		printf("\tsched_setaffinity failed! cpu: %d status: %d\n",
		    trial->id, err);

	while (1) {

		/* Wait for instructions */
		pthread_mutex_lock(trial->mutex);
		while (trial->state == WAIT)
			pthread_cond_wait(trial->cv, trial->mutex);

		if (trial->state == TERMINATE) {
			pthread_mutex_unlock(trial->mutex);
			pthread_exit(NULL);
			break;
		}

		pthread_mutex_unlock(trial->mutex);

		F = trial->F;
		nF = trial->nF;

		x = trial->x;
		x_F = trial->x_F;

		/* Set up list of infeasibles if necessary */
		if (!trial->H1)
			trial->H1 = (long*)malloc((trial->nF)*sizeof(long));
		else
			trial->H1 = (long*)realloc(trial->H1, 
			    (trial->nF)*sizeof(long));
		trial->nH1 = 0;

		/* Allocate trial descent if necessary */
		if (!trial->x_c)
			trial->x_c = cholmod_l_allocate_dense(nF, 1, nF,
			    CHOLMOD_REAL, trial->c);
		else 
			assert(trial->x_c->nrow == trial->nF);

		/* Calculate trial descent */
		for (i = 0; i < trial->nF; i++) {
			xptr = &((double*)(trial->x_c->x))[i];
			*xptr = (1.0 - trial->alpha[0]) * 
			    ((double*)(x->x))[F[i]] + trial->alpha[0] *
			    ((double*)(x_F->x))[i];

			/*
			 * Project into feasible space, keeping track of which
			 * coefficients become constrained.
			 */

			if (*xptr < 0.0) {
				*xptr = 0.0;
				trial->H1[trial->nH1++] = F[i];
			}
		}

		/* Calculate residual for this descent */
		trial->residual = calc_residual(trial->AtA_F, trial->Atb_F,
		    trial->x_c, trial->c);

		/* Report back */
		pthread_mutex_lock(trial->mutex);
		trial->state = WAIT;
		pthread_cond_broadcast(trial->cv);
		pthread_mutex_unlock(trial->mutex);
	}
}

int
get_nthreads(void)
{
	char* str;
	int nthreads = -1;

	str = getenv("GOTO_NUM_THREADS");
	if (str == NULL)
		str = getenv("OMP_NUM_THREADS");
	if (str == NULL)
		goto fail;

	nthreads = atoi(str);

	if (nthreads >= 1)
		return (nthreads);

fail:
	/* Use the number of CPUs in the system */
#ifdef _SC_NPROCESSORS_ONLN
	nthreads = sysconf(_SC_NPROCESSORS_ONLN);
#else
	#warn Do not know how to determine CPU count on this platform
#endif
	if (nthreads < 1)
		nthreads = 1;

	return (nthreads);
}

/* Sort doubles in descending order */
static int
double_rcmp(const void *xa, const void *xb)
{
	const double *a, *b;
	a = xa; b = xb;

	if (*a < *b)
		return (1);
	else if (*a > *b)
		return (-1);
	else
		return (0);
}

int
walk_descents(cholmod_sparse *AtA_F, 
    cholmod_dense *Atb_F, cholmod_dense *x, 
    cholmod_dense *x_F, long *F, long *nF_, long *H1, long *nH1_,
    double *residual, int *residual_calcs, int verbose, cholmod_common *c)
{

	long nF, nH1;
	clock_t t0, t1;

	double *alpha, res=NAN;
	int n_alpha, n_threads, n_blocks, success, feasible;
	int i, j=-2, k;
	pthread_t *threads;
	pthread_attr_t thread_attr;
	descent_trial *descent_trials;
				
	nF = *nF_;
	nH1 = *nH1_;

	t0 = clock();

	feasible = false;

	alpha = (double*)malloc((nF+2) 
	    * sizeof(double));
	alpha[0] = 0;
	alpha[1] = 1;
	n_alpha = 2;

	/*
	 * For each infeasible coordinate in x_F[i], 
	 * alpha is the relative distance along the
	 * descent vector between the current solution
	 * and the origin in dimension F[i].
	 */
	for (i = 0; i < nF; i++) {
		/* alpha = x[F]/(x[F]-x_F) */
		if (((double*)(x_F->x))[i] < 0) {
			alpha[n_alpha] = 
			    ((double*)(x->x))[F[i]] /
			    (((double*)(x->x))[F[i]]
			    -((double*)(x_F->x))[i]);
			if ((alpha[n_alpha] < 1) &&
			    (alpha[n_alpha] > 
			    0)) ++n_alpha;
		}
	}

	/* Sort the scales (except 0 and 1) in descending order */
	qsort(&alpha[2], n_alpha-2, sizeof(alpha[0]),
	    double_rcmp);

	/* Begin threadening */
	n_threads = get_nthreads();
	threads = (pthread_t*)malloc(n_threads * 
	    sizeof(pthread_t));

	descent_trials = (descent_trial*)malloc(
	    n_threads*sizeof(descent_trial));

	/* Set up thread attributes */
	pthread_attr_init(&thread_attr);
	pthread_attr_setdetachstate(&thread_attr, PTHREAD_CREATE_JOINABLE);

	pthread_mutex_t mutex;
	pthread_cond_t cv;
	pthread_mutex_init(&mutex, NULL);
	pthread_cond_init(&cv, NULL);

	/* Initialize argument structs for each thread */
	descent_trials[0].x        = x;
	descent_trials[0].x_F      = x_F;
	descent_trials[0].AtA_F    = AtA_F;
	descent_trials[0].Atb_F    = Atb_F;
	descent_trials[0].c        = c;
	descent_trials[0].F        = F;
	descent_trials[0].nF       = nF;
	descent_trials[0].alpha    = NULL;
	descent_trials[0].x_c      = NULL;
	descent_trials[0].residual = 0;
	descent_trials[0].H1       = NULL;
	descent_trials[0].nH1      = 0;
	descent_trials[0].state    = WAIT;
	descent_trials[0].mutex    = &mutex;
	descent_trials[0].cv       = &cv;
	descent_trials[0].id       = 0;
	for (i = 1; i < n_threads; i++) {
		memcpy(&descent_trials[i], &descent_trials[0],
		    sizeof(descent_trial));
		descent_trials[i].id = i;
	}
							
	n_blocks = (int)ceil(n_alpha/((double)(n_threads)));

	/* Spawn worker threads */
	for (i = 0; i < n_threads; i++)
		pthread_create(&threads[i], &thread_attr, 
		    (void*)(&evaluate_descent), (void*)(&descent_trials[i]));

	/*
	 * Find the furthest distance we can move along
	 * the descent vector while still reducing the
	 * magnitude of the residual vector.
	 *
	 * Calculate the residual in blocks of n_threads alphas each.
	 */
	success = false;
	for (i = 0; i < n_blocks; i++) {
		if (success)
			break;
		/* Set alpha and start threads running */
		pthread_mutex_lock(&mutex);
		for (j = 0; j < n_threads; j++) {
			if (i*n_threads + j >= n_alpha)
				break;
			descent_trials[j].alpha =
			    &alpha[i*n_threads + j];
			descent_trials[j].state = RUN;
		} 
		pthread_cond_broadcast(&cv);
		pthread_mutex_unlock(&mutex);

		/* Wait for threads to finish calculations */
		int done = false;
		pthread_mutex_lock(&mutex);
		while (!done) {
			pthread_cond_wait(&cv, &mutex);
			done = true;
			for (j = 0; j < n_threads; j++) {
				if (i*n_threads + j >= n_alpha)
					break;
				if (descent_trials[j].state != WAIT)
					done = false;
			}
		}
		pthread_mutex_unlock(&mutex);

		(*residual_calcs) += n_threads;

		/*
		 * See if we reduced the residual
		 */
		for (j = 0; j < n_threads; j++) {
			if (i*n_threads + j >= n_alpha)
				break;
			/* 
			 * NB: the first alpha is zero, which gives us the
			 * residual at the current solution. Use this for
			 * comparisons.
			 */
			if ((i == 0) && (j == 0)) {
				assert( descent_trials[j].alpha[0] == 0 );
				res = descent_trials[j].residual;
			} else if ((descent_trials[j].residual < res) ||
			    (i*n_threads + j == n_alpha-1)) {
				success = true;
				/* 
			 	 * NB: x_c is feasible by construction.
			 	 */
				for (k = 0; k < nF; k++) {
					((double*)(x->x))[F[k]] =
					    ((double*)
					    (descent_trials[j].x_c->x))[k];
				}
				assert( nH1 == 0 );
				for (k = 0; k < descent_trials[j].nH1; k++) 
					H1[nH1++] = descent_trials[j].H1[k];
				
				/*
				 * NB: if we've arrived at the last (smallest)
				 * alpha and have not reduced the residual, we
				 * bind the negative coefficients in the trial
				 * solution and try again. This is equivalent 
				 * to the Lawson-Hanson single-pivot 
				 * interpolation step.
				 */
				if (descent_trials[j].residual < res) {
					feasible = true;
					*residual = descent_trials[j].residual;
				} else
					feasible = false;

				if (verbose)
					printf("\talpha[%d] = "
					    "%.3e, d_res = "
					    "%e\n", i*n_threads + j, 
					    alpha[i*n_threads + j],
					    descent_trials[j].residual - res);
				break;
			}

		} /* for (n_threads) */
	} /* for (n_blocks) */

	/* Shut down worker threads */
	pthread_mutex_lock(&mutex);
	for (k = 0; k < n_threads; k++) 
		descent_trials[k].state = TERMINATE;
	pthread_cond_broadcast(&cv);
	pthread_mutex_unlock(&mutex);

	/* Wait for threads to exit */
	for (k = 0; k < n_threads; k++) 
		pthread_join(threads[k], NULL);

	/* Clean up thread data */
	for (k = 0; k < n_threads; k++) {
		if (descent_trials[k].H1)
			free(descent_trials[k].H1);
		if (descent_trials[k].x_c)
			cholmod_l_free_dense(
			    &(descent_trials[k].x_c), c);
	}

	/* Clean up pthreads-related detritus */
	pthread_cond_destroy(&cv);
	pthread_mutex_destroy(&mutex);
	pthread_attr_destroy(&thread_attr);
	free(descent_trials);
	free(threads);
	free(alpha);

	if (verbose) {
		t1 = clock();
		printf("\tCompare %d descents: "
		    "%.2f s\n", i*n_threads + j-1, (double)(t1-t0) / 
		    (CLOCKS_PER_SEC));
	}

	*nH1_ = nH1;
	return (feasible);
}

