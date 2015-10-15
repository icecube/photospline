/*
 * nnls.c: Sparse non-negative least squares algorithms. Except for
 *  Lawson-Hanson, all algorithms requrie full-rank problems
 *  (exactly or over-determined). For well-conditioned full-rank problems,
 *  BLOCK3 is the fastest.
 */

#include <cholmod.h>
#include <SuiteSparseQR_C.h>

#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#include <time.h>
#include <float.h>
#include <assert.h>

#include "photospline/detail/splineutil.h"
#include "cholesky_solve.h"

#define MAX_TRIALS 5
#define N_RESOLVES 0
#define KKT_TOL 1e-6

static int intcmp(const void *xa, const void *xb);

/*
 * Lawson-Hanson NNLS
 *
 * Algorithm NNLS from "Solving Least Squares Problems", Charles Lawson and
 *  Richard Hanson. Prentice-Hall, 1974.
 */

cholmod_dense *
nnls_lawson_hanson(cholmod_sparse *A, cholmod_dense *y, double tolerance,
    int min_iterations, int max_iterations, unsigned int npos, int normaleq,
    int verbose, cholmod_common *c)
{
        cholmod_dense *x, *w, *p, *yp;
        cholmod_sparse *Ap;
        long P[A->ncol], Z[A->ncol];
        unsigned int nP, nZ;
        double wmax, wpmin, alpha, qtemp;
        int i, j, n, t, qmax;
        int last_freed = -1;

        /* By default, all coefficients are positive */
        if (npos == 0)
                npos = A->ncol;

        /* Step 1: All elements actively constrained */
        nP = A->ncol - npos; nZ = npos;
        for (i = 0; i < nZ; i++)
                Z[i] = i;
        for (i = 0; i < nP; i++)
                P[i] = npos + i;

        /* Initialize coefficient and negative gradient vectors */
        x = cholmod_l_zeros(A->ncol, 1, CHOLMOD_REAL, c);
        w = cholmod_l_zeros(A->ncol, 1, CHOLMOD_REAL, c);

        for (n = 0; n < max_iterations || max_iterations == 0; n++) { 
                /* Step 2: compute the negative gradient of the residuals, w */
                if (normaleq) {
                         /* If the normal equations are pre-formulated, 
                          * w = y - Ax. */
                        double alpha[2] = {1.0, 0.0}, beta[2] = {-1.0, 0.0};

                        memcpy(w->x, y->x, sizeof(double)*A->ncol);
                        cholmod_l_sdmult(A, 0 /* no transpose */, beta, alpha,
                            x, w, c);
                } else {
                        /* Otherwise, take the derivative: w = At(y - Ax) */
                        cholmod_dense *wtemp;
                        double alpha[2] = {1.0, 0.0}, beta[2] = {-1.0, 0.0};

                        wtemp = cholmod_l_copy_dense(y, c);
                        cholmod_l_sdmult(A, 0 /* no transpose */, beta, alpha,
                            x, wtemp, c);
                        beta[0] = 0;
                        cholmod_l_sdmult(A, 1 /* transpose */, alpha, beta,
                            wtemp, w, c);
                        cholmod_l_free_dense(&wtemp, c);
                }

                /* Step 3a: Check for completion */
                if (nZ == 0)
                        break;

                /* Steps 3b,4: Find maximum w, complete if <= 0 */
                wmax = ((double *)(w->x))[Z[0]];
                t = 0;
                for (i = 1; i < nZ; i++) {
                        if (((double *)(w->x))[Z[i]] > wmax &&
                            last_freed != Z[i]) {
                                t = i;
                                wmax = ((double *)(w->x))[Z[t]];
                        }
                }

                if (wmax <= 0)
                        break;

                /* See if we might be within fit tolerance */
                if (wmax < tolerance && n >= min_iterations) {
                        /* Check if any passive coefficients need to be reduced
                         * due to clipped ringing */
                        if (nP == 0)
                                break;
                        
						assert(nP>0);
                        wpmin = ((double *)(w->x))[P[0]];
                        for (i = 1; i < nP; i++) {
                                if (((double *)(w->x))[P[i]] < wpmin)
                                        wpmin = ((double *)(w->x))[P[i]];
                        }

                        if (-wpmin < tolerance)
                                break;
                }

                /* Step 5: Move coefficient Z[t] into the passive set P */
                if (verbose)
                        printf("Freeing coefficient %ld (active: %d, "
                            "passive: %d, wmax: %e)\n", Z[t], nZ, nP, wmax);
                last_freed = Z[t];
                alpha = -1;
                P[nP++] = Z[t];
                nZ--;
                for (i = t; i < nZ; i++)
                        Z[i] = Z[i+1];

                /*
                 * Steps 6-11: Move coefficients from the passive to the
                 * active set
                 */
                while (1) {        
                        /*
                         * Step 6: Solve the least-squares subproblem for the
                         * columns of A in the passive set P by QR
                         * decomposition.
                         */
                        if (normaleq) {
                                Ap = cholmod_l_submatrix(A, P, nP, P, nP, 1,
                                    1, c);
                                yp = cholmod_l_allocate_dense(nP, 1, nP,
                                    CHOLMOD_REAL, c);
                                for (i = 0; i < nP; i++)
                                        ((double *)(yp->x))[i] =
                                            ((double *)(y->x))[P[i]];
                                p = SuiteSparseQR_C_backslash_default(Ap, yp, c);
                                cholmod_l_free_sparse(&Ap, c);
                                cholmod_l_free_dense(&yp, c);
                        } else {
                                Ap = cholmod_l_submatrix(A, NULL, -1, P, nP,
                                    1, 1, c);
                                p = SuiteSparseQR_C_backslash_default(Ap, y, c);
                                cholmod_l_free_sparse(&Ap, c);
                        }

                        /*
                         * Step 7: Check if any coefficients need be constrained
                         */
                        for (i = 0; i < nP; i++)
                                if (P[i] < npos && ((double *)(p->x))[i] <= 0)
                                        break;
                        if (i == nP) {
                                /*
                                 * All were positive. Cycle back for the next.
                                 */
                                bzero(x->x, sizeof(double)*x->nrow);
                                for (i = 0; i < nP; i++)
                                        ((double *)(x->x))[P[i]] =
                                            ((double *)(p->x))[i];
                                cholmod_l_free_dense(&p, c);
                                break; /* Break loop to step 2 */
                        }

                        /*
                         * Step 8-9: Compute q, alpha
                         */
                        alpha = 2; /* All computed values <= 1 */
                        qmax = -1;
                        for (i = 0; i < nP; i++) {
                                if (P[i] >= npos || ((double *)(p->x))[i] > 0)
                                        continue;

                                qtemp = ((double *)(x->x))[P[i]]/
                                    (((double *)(x->x))[P[i]] -
                                    ((double *)(p->x))[i]);

                                if (qtemp < alpha && qtemp != 0) {
                                        qmax = P[i];
                                        alpha = qtemp;
                                } else if (last_freed == P[i]) {
                                        /* Anti-cycling advice from LH */
                                        alpha = 0;
                                        qmax = P[i];
                                        break;
                                }
                        }

                        if (qmax < 0) {
                                fprintf(stderr, "%s line %d: Math has failed\n",
                                    __FILE__, __LINE__);
                                exit(1);
                        }

                        /*
                         * Step 10: Recompute x
                         */
                        for (i = 0; i < nP; i++)
                                ((double *)(x->x))[P[i]] += alpha*
                                    (((double *)(p->x))[i] -
                                    ((double *)(x->x))[P[i]]);

                        /* Avoid rounding errors above */
                        ((double *)(x->x))[qmax] = 0;
                        cholmod_l_free_dense(&p, c);

                        /*
                         * Step 11: Move coefficients equal to zero to the
                         * active set.
                         */
                        for (i = 0; i < nP; i++) {
                                if (P[i] >= npos || ((double *)(x->x))[P[i]] > 0)
                                        continue;

                                if (verbose)
                                        printf("\tConstraining coefficient %ld "
                                            "(active: %d, passive: %d, "
                                            "value: %e)\n", P[i], nZ, nP,
                                            ((double *)(x->x))[P[i]]);

                                ((double *)(x->x))[P[i]] = 0;
                                Z[nZ++] = P[i];
                                nP--;
                                for (j = i; j < nP; j++)
                                        P[j] = P[j+1];
                                i--;
                        }

                        /* If alpha = 0, we've reached equilibrium */
                        if (alpha == 0)
                                break;
                }

                /* Exit to the caller in equilibrium */
                if (alpha == 0)
                        break;
        }

        /* Step 12: return */
        cholmod_l_free_dense(&w, c);
        return (x);
}

/*
 * An implementation of the Portugal/Judice/Vicente block-pivoting algorithm for
 * pre-formulated normal equations
 *
 *    See:
 *       A Comparison of Block Pivoting and Interior-Point Algorithms for
 *       Linear Least Squares Problems with Nonnegative Variables
 *
 *       Author(s): Luis F. Portugal, Joaquim J. Judice, Luis N. Vicente
 *       Source: Mathematics of Computation, Vol. 63, No. 208 (Oct., 1994),
 *         pp. 625-643
 *       Published by: American Mathematical Society
 *       Stable URL: http://www.jstor.org/stable/2153286
 */

cholmod_dense *
nnls_normal_block(cholmod_sparse *AtA, cholmod_dense *Atb, int verbose,
   cholmod_common *c)
{
        int nvar = AtA->nrow;
        long F[nvar], G[nvar];
        int H1[nvar], H2[nvar];
        cholmod_dense *x, *y, *x_F, *Atb_F;
        cholmod_sparse *AtA_F;
        int nF, nG, nH1, nH2, ninf;
        int i, j, k, trials, murty_steps, iter;

        /* XXX: make these settable? */
        iter = 3*nvar;                /* Maximum number of iterations */
        trials = MAX_TRIALS;        /* Runs without progress before reverting
                                 * to a deterministic algorithm */
        murty_steps = MAX_TRIALS;

        x = cholmod_l_zeros(nvar, 1, CHOLMOD_REAL, c);
        y = cholmod_l_allocate_dense(nvar, 1, nvar, CHOLMOD_REAL, c);

        /*
         * Initialize variables such that all coefficients are in the
         * active set.
         */

        nF = 0;                        /* No feasible coefficients */
        nG = nvar;
        for (i = 0; i < nvar; i++)
                G[i] = i;
        ninf = nvar + 1;        /* Number of infeasible coefficients */

        /* Drop small entries from the problem */
        //cholmod_l_drop(DBL_EPSILON, AtA, c);

        /* Set up the dual vector */
        for (i = 0; i < nvar; i++)
                ((double *)(y->x))[i] = -((double *)(Atb->x))[i];

        while (iter-- > 0) {
                /*
                 * Fill H1, H2 with the points that must be moved from the
                 * passive to active and active to passive sets, respectively.
                 */

                nH1 = nH2 = 0;
                for (i = 0; i < nF; i++)
                        if (((double *)(x->x))[F[i]] < -KKT_TOL)
                                H1[nH1++] = F[i];
                for (i = 0; i < nG; i++)
                        if (((double *)(y->x))[G[i]] < -KKT_TOL)
                                H2[nH2++] = G[i];

                /*
                 * If all coefficients were positive, then we are done.
                 */

                if (nH1 == 0 && nH2 == 0)
                        break;

                if (ninf <= murty_steps)
                        trials = -1;
                /*
                 * Check the status of the bulk set switching.
                 * After MAX_TRIALS iterations of Murty's finite method,
                 * revert to block switching.
                 */
                
                if (ninf > murty_steps &&
                    (nH2 + nH1 < ninf || trials < -murty_steps)) {
                        if (nH2 + nH1 <= ninf)
                            murty_steps++;
                        ninf = nH2 + nH1;
                        trials = MAX_TRIALS;
                } else {
                        /* Stuck, check if we need to try something else */
                        trials--;
                        if (verbose)
                                printf("Stuck! trials: %d nH1: %d nH2: %d\n",trials,nH1,nH2);
                        if (trials < 0) {
                                /*
                                 * Out of luck. Fall back to slow but
                                 * guaranteed method (picking the last
                                 * infeasible coordinate).
                                 */
                                
                                if (nH2 == 0) {
                                        goto maxh1;
                                } else if (nH1 == 0) {
                                        goto maxh2;
                                } else if (H1[nH1 - 1] > H2[nH2 - 1]) {
                                        maxh1:
                                        H1[0] = H1[nH1 - 1];
                                        nH1 = 1; nH2 = 0;
                                        if (verbose)
                                                printf("H1: %d (%e)\n",H1[0],((double *)(x->x))[H1[0]]);
                                } else {
                                        maxh2:
                                        H2[0] = H2[nH2 - 1];
                                        nH2 = 1; nH1 = 0;
                                        if (verbose)
                                                printf("H2: %d (%e)\n",H2[0],((double *)(y->x))[H2[0]]);
                                }
                        }
                }

                if (verbose)
                        printf("Iteration %d Infeasibles: %d\n", (3*nvar - iter), ninf);

                /*
                 * Next, remove elements in H1 from F, and add them to G,
                 * exploiting the fact that H1 elements are in order
                 * relative to their order in F.
                 */
                for (i = 0, j = 0; i < nH1; i++) {
                        G[nG++] = H1[i];
                        while (F[j] != H1[i]) j++;
                        for (k = j+i; k+1 < nF; k++)
                                F[k-i] = F[k-i+1];
                }
                nF -= nH1;

                /* And vice versa */
                for (i = 0, j = 0; i < nH2; i++) {
                        F[nF++] = H2[i];
                        while (G[j] != H2[i]) j++;
                        for (k = j+i; k+1 < nG; k++)
                                G[k-i] = G[k-i+1];
                }
                nG -= nH2;

                qsort(G, nG, sizeof(G[0]), intcmp);
                qsort(F, nF, sizeof(F[0]), intcmp);

                /* Solve the unconstrained part */

                if (verbose)
                    printf("Unconstrained solve for %d of %d coefficients\n",
                      nF, nvar);
                AtA_F = cholmod_l_submatrix(AtA, F, nF, F, nF, 1, 1, c);
                AtA_F->stype = 1;
                Atb_F = cholmod_l_allocate_dense(nF, 1, nF, CHOLMOD_REAL, c);
                for (i = 0; i < nF; i++)
                        ((double *)(Atb_F->x))[i] = ((double *)(Atb->x))[F[i]];

                /* Solve the system AtA_F*x = Atb_F, refining the solution iteratively */
                x_F = cholesky_solve(AtA_F, Atb_F, c, verbose, N_RESOLVES);

                for (i = 0; i < nF; i++)
                        ((double *)(x->x))[F[i]] = ((double *)(x_F->x))[i];
                cholmod_l_free_sparse(&AtA_F, c);
                cholmod_l_free_dense(&Atb_F, c);

                /* Update the constrained part */

                for (i = 0; i < nG; i++)
                        ((double *)(x->x))[G[i]] = 0;
                for (i = 0; i < nF; i++)
                        ((double *)(y->x))[F[i]] = 0;

                
                {
                        cholmod_sparse *AtA_FG;
                        cholmod_dense *Atb_G;
                        double ones[2] = {1., 0}, mones[2] = {-1., 0};

                        Atb_G = cholmod_l_allocate_dense(nG, 1, nG,
                            CHOLMOD_REAL, c);
                        for (i = 0; i < nG; i++)
                                ((double *)(Atb_G->x))[i] =
                                    ((double *)(Atb->x))[G[i]];

                        AtA_FG = cholmod_l_submatrix(AtA, G, nG, F, nF, 1,1,c);
                        cholmod_l_sdmult(AtA_FG, 0, ones, mones, x_F, Atb_G, c);
                        
                        for (i = 0; i < nG; i++)
                                ((double *)(y->x))[G[i]] =
                                    ((double *)(Atb_G->x))[i];

                        cholmod_l_free_dense(&Atb_G, c);
                        cholmod_l_free_sparse(&AtA_FG, c);
                }
                
#if 0
                double* x_full = (double*)(x->x);
                double* y_full = (double*)(y->x);
                DPRINT(x_full, nvar, " %- .1e");
                DPRINT(y_full, nvar, " %- .1e");
#endif
                cholmod_l_free_dense(&x_F, c);
        }

        cholmod_l_free_dense(&y, c);

        return (x);
}

/*
 * An implementation of the Portugal/Judice/Vicente block-pivoting algorithm for
 * pre-formulated normal equations, with single-row up/down-dates.
 *
 *    See:
 *       A Comparison of Block Pivoting and Interior-Point Algorithms for
 *       Linear Least Squares Problems with Nonnegative Variables
 *
 *       Author(s): Luis F. Portugal, Joaquim J. Judice, Luis N. Vicente
 *       Source: Mathematics of Computation, Vol. 63, No. 208 (Oct., 1994),
 *         pp. 625-643
 *       Published by: American Mathematical Society
 *       Stable URL: http://www.jstor.org/stable/2153286
 */


cholmod_dense *
nnls_normal_block_updown(cholmod_sparse *AtA, cholmod_dense *Atb, int verbose,
   cholmod_common *c)
{
        int nvar = AtA->nrow;
        long *F, *G, *H1, *H2;
        cholmod_dense *x, *y;
        cholmod_factor *L;
        long nF, nG, nH1, nH2, ninf;
        int i, trials, murty_steps, iter;
        clock_t t0, t1;

        /* XXX: make these settable? */
        iter = 3*nvar;                /* Maximum number of iterations */
        trials = MAX_TRIALS;        /* Runs without progress before reverting
                                 * to a deterministic algorithm */
        murty_steps = MAX_TRIALS;

        /*
         * The METIS graph-partitioning code tends to segfault. Don't use it.
         */
        c->nmethods = 2;

        F  = (long*)malloc(sizeof(long)*nvar);
        G  = (long*)malloc(sizeof(long)*nvar);
        H1 = (long*)malloc(sizeof(long)*nvar);
        H2 = (long*)malloc(sizeof(long)*nvar);

        x = cholmod_l_zeros(nvar, 1, CHOLMOD_REAL, c);
        y = cholmod_l_allocate_dense(nvar, 1, nvar, CHOLMOD_REAL, c);

        /*
         * Initialize variables such that all coefficients are in the
         * active set.
         */

        nF = 0;                        /* No feasible coefficients */
        nG = nvar;
        for (i = 0; i < nvar; i++)
                G[i] = i;
        ninf = nvar + 1;        /* Number of infeasible coefficients */

        /* Drop small entries from the problem, replacing them with eels. 
         * Also drops the lower half of AtA if AtA->stype is 1. */
        cholmod_l_drop(DBL_EPSILON, AtA, c);
        AtA->stype = 1;

        /* Set up the dual vector */
        for (i = 0; i < nvar; i++)
                ((double *)(y->x))[i] = -((double *)(Atb->x))[i];

        L = NULL;
        t0 = clock();

        while (iter-- > 0) {
                /*
                 * Fill H1, H2 with the points that must be moved from the
                 * passive to active and active to passive sets, respectively.
                 */

                nH1 = nH2 = 0;
                for (i = 0; i < nF; i++)
                        if (((double *)(x->x))[F[i]] < -KKT_TOL)
                                H1[nH1++] = F[i];
                for (i = 0; i < nG; i++)
                        if (((double *)(y->x))[G[i]] < -KKT_TOL)
                                H2[nH2++] = G[i];

                /*
                 * If all coefficients were positive, then we are done.
                 */

                if (nH1 == 0 && nH2 == 0)
                        break;

                if ((ninf <= murty_steps) )
                        trials = -1;

                /*
                 * Check the status of the bulk set switching.
                 * After MAX_TRIALS iterations of Murty's finite method,
                 * revert to block switching.
                 */
                
                if (ninf > murty_steps && nH2 + nH1 < ninf) {
                        if (nH2 + nH1 <= ninf)
                            murty_steps++;
                        ninf = nH2 + nH1;
                        trials = MAX_TRIALS;
                } else {
                        /* Stuck, check if we need to try something else */
                        trials--;
                        if (verbose)
                                printf("\tStuck! trials: %d nH1: %ld nH2: %ld"
                                    "\n", trials, nH1, nH2);
                        if (trials < 0) {
                                /*
                                 * Out of luck. Fall back to slow but
                                 * guaranteed method (picking the last
                                 * infeasible coordinate).
                                 */
                                
                                if (nH2 == 0) {
                                        goto maxh1;
                                } else if (nH1 == 0) {
                                        goto maxh2;
                                } else if (H1[nH1 - 1] > H2[nH2 - 1]) {
                                        maxh1:
                                        H1[0] = H1[nH1 - 1];
                                        nH1 = 1; nH2 = 0;
                                        if (verbose)
                                                printf("\tH1: %ld (%e)\n",
                                                    H1[0],
                                                    ((double *)(x->x))[H1[0]]);
                                } else {
                                        maxh2:
                                        H2[0] = H2[nH2 - 1];
                                        nH2 = 1; nH1 = 0;
                                        if (verbose)
                                                printf("\tH2: %ld (%e)\n",
                                                    H2[0],
                                                    ((double *)(y->x))[H2[0]]);
                                }
                        }
                }

                if (verbose)
                        printf("Iteration %d Infeasibles: %ld\n",
                            (3*nvar - iter), ninf);

                /* We're about to recompute x, so free it */
                cholmod_l_free_dense(&x, c);

                if (verbose)
                    printf("\tUnconstrained solve for %ld of %d coefficients\n",
                      nF - nH1 + nH2, nvar);

                /*
                 * Next, update the factorization of AtA[:,F][F,:], removing
                 * rows in H1 and adding those in H2. For large updates, the
                 * factorization will be recomputed from scratch; for small
                 * updates, the rows will be added and deleted from the
                 * factorization one-by-one, potentially saving a lot of
                 * computing time.
                 *
                 * F, G, H1, H2, and the associated counts are updated on
                 * exit.
                 */
                L = modify_factor(AtA, L, F, &nF, G, &nG, H1, &nH1, H2, &nH2,
                    verbose, c);
        
                if (verbose)
                        t0 = clock();
                
                /*
                 * Solve the full system, but with the rows of L corresponding
                 * to G set to identity. 
                 */
                if (L->n == nvar) {
                        x = cholmod_l_solve(CHOLMOD_A, L, Atb, c);
                } else {
                        /*
                         * If L is only a subset of the full matrix, solve
                         * only the passive set.
                         */
                        cholmod_dense *Atb_F, *x_F;

                        Atb_F = cholmod_l_allocate_dense(nF, 1, nF,
                            CHOLMOD_REAL, c);
                        for (i = 0; i < nF; i++)
                                ((double*)(Atb_F->x))[i] =
                                    ((double*)(Atb->x))[F[i]];

                        x_F = cholmod_l_solve(CHOLMOD_A, L, Atb_F, c);

                        cholmod_l_free_dense(&Atb_F, c);
                        x = cholmod_l_allocate_dense(nvar, 1, nvar,
                            CHOLMOD_REAL, c);
                        for (i = 0; i < nF; i++)
                                ((double*)(x->x))[F[i]] =
                                    ((double*)(x_F->x))[i];
                        cholmod_l_free_dense(&x_F, c);
                }

                if (verbose) {
                        t1 = clock();
                        printf("\tSolve[%d] (%ld free): %.2f s\n",nvar,nF,
                            (double)(t1-t0)/(CLOCKS_PER_SEC));
                }

                /* Explicitly zero x[G] and y[F] */
                for (i = 0; i < nG; i++)
                        ((double *)(x->x))[G[i]] = 0;
                for (i = 0; i < nF; i++)
                        ((double *)(y->x))[F[i]] = 0;

                /* Update the constrained part */
                {
                        cholmod_sparse *AtA_FG;
                        cholmod_dense *Atb_G, *x_F;
                
                        /*
                         * NB: We happen to know that AtA is actually
                         * symmetric, but cholmod_submatrix will only
                         * handle unsymmetric ones. The submatrix constructed
                         * this way will also be secretly symmetric.
                         */
                        AtA->stype = 0;
                        double ones[2] = {1., 0}, mones[2] = {-1., 0};

                        if (verbose)
                                t0 = clock();

                        Atb_G = cholmod_l_allocate_dense(nG, 1, nG,
                            CHOLMOD_REAL, c);
                        for (i = 0; i < nG; i++)
                                ((double *)(Atb_G->x))[i] =
                                    ((double *)(Atb->x))[G[i]];

                        x_F = cholmod_l_allocate_dense(nF, 1, nF,
                            CHOLMOD_REAL, c);
                        for (i = 0; i < nF; i++)
                                ((double *)(x_F->x))[i] =
                                    ((double *)(x->x))[F[i]];

                        AtA_FG = cholmod_l_submatrix(AtA, G, nG, F, nF, 1,1,c);
                        cholmod_l_sdmult(AtA_FG, 0, ones, mones, x_F, Atb_G, c);
                        
                        for (i = 0; i < nG; i++)
                                ((double *)(y->x))[G[i]] =
                                    ((double *)(Atb_G->x))[i];

                        cholmod_l_free_dense(&Atb_G, c);
                        cholmod_l_free_dense(&x_F, c);
                        cholmod_l_free_sparse(&AtA_FG, c);

                        if (verbose) {
                                t1 = clock();
                                printf("\tUpdate y[%ld]: %.2f s\n",
                                    nG,(double)(t1-t0)/(CLOCKS_PER_SEC));
                        }
                        AtA->stype = 1;
                        
                }
#if 0
                double* x_full = (double*)(x->x);
                double* y_full = (double*)(y->x);
                DPRINT(x_full, nvar, " %- .1e");
                DPRINT(y_full, nvar, " %- .1e");
#endif
        }

        cholmod_l_free_dense(&y, c);
        cholmod_l_free_factor(&L, c);

        free(F);
        free(G);
        free(H1);
        free(H2);

        return (x);
}

/*
 * An implementation of Adlers' BLOCK3 algorithm for NNLS with
 * pre-formulated normal equations and single row additions/deletions.
 *
 * This algorithm always reduces the magnitude of the residual vector
 * and is thus finite.
 *
 *     See:
 *         http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.32.8173
 */


cholmod_dense *
nnls_normal_block3(cholmod_sparse *AtA, cholmod_dense *Atb, int verbose,
   cholmod_common *c)
{
        int nvar = AtA->nrow;
        long *F, *G, *H1, *H2;
        long *F_, *G_, *Fprime, *Gprime;
        cholmod_dense *x, *y, *x_F, *Atb_F;
        cholmod_sparse *AtA_F;
        cholmod_factor *L;
        long nF, nG, nH1, nH2, ninf;
        long nFprime, nGprime, nF_, nG_;
        int i, j, k;
        int iter, max_iter, solves, residual_calcs;
        int feasible;
        clock_t t0, t1;
        double kkt_tolerance, y_min, residual;

        /* XXX: make these settable? */
        max_iter = 120;                /* Maximum number of iterations */
        solves = 0;
        residual_calcs = 0;
        residual = DBL_MAX;

        /* Heuristic stopping tolerance inspired from Adlers' thesis */
        kkt_tolerance = ((double)(nvar)) * DBL_EPSILON * 1e5;
        if (verbose)
                printf("Stopping tolerance: %e\n",kkt_tolerance);

        /* F and G track the state of the factorization */
        F  = (long*)malloc(sizeof(long)*nvar);
        G  = (long*)malloc(sizeof(long)*nvar);
        /* 
         * When the factorization is out of sync with the column
         * partition, Fprime and Gprime can be used instead.
         */
        Fprime = (long*)malloc(sizeof(long)*nvar);
        Gprime = (long*)malloc(sizeof(long)*nvar);
        /* 
         * H1 and H2 hold pending changes to the column partition
         * of the factorization.
         */
        H1 = (long*)malloc(sizeof(long)*nvar);
        H2 = (long*)malloc(sizeof(long)*nvar);

        nF = nG = nH1 = nH2 = 0;
        nGprime = -1;

        t0 = clock();

        /* Drop small entries from the problem, replacing them with eels. 
         * Also drops the lower half of AtA if AtA->stype is 1. */
        cholmod_l_drop(DBL_EPSILON, AtA, c);
        AtA->stype = 1;

        /*
         * This algorithm is under a tremendous amount of memory pressure, 
         * especially in the initial solve, which operates on the full matrix.
         * NESDIS can sometimes find a much better ordering than AMD, but is 
         * ~ 5x slower. The extra time spend trying different orderings is 
         * well worth ~ 2x savings in memory.
         * 
         * XXX: METIS graph-partitioning is known to segfault on large
         * matrices, so we remove it from the list of orderings.
         */
        c->method[2] = c->method[3]; /* NESDIS, default parameters */
        c->nmethods = 3;

	/*
	 * Start with x equal to the null vector. This implies
	 * that y = -Atb.
	 */
	AtA_F = NULL; x_F = NULL; Atb_F = NULL; L = NULL; 
	x = cholmod_l_zeros(nvar, 1, CHOLMOD_REAL, c); 	
	y = cholmod_l_copy_dense(Atb, c);
	for (i = 0; i < nvar; i++) {
		G[nG++] = i;
		((double*)(y->x))[i] *= -1.0;
	}

        for (iter = 0; iter < max_iter; iter++) {

                nH2 = 0;
                /*
                 * Fill H2 with the points that have large negative Lagrange
                 * multipliers and thus have to be released into the passive
                 * set.
                 *
                 * If the state of the factorization is out of sync
                 * with the current column partition, use Gprime to 
                 * fill H2. Otherwise, use G.
                 */

                if (nGprime < 0) {
                        G_ = G;
                        nG_ = nG;
                } else {
                        G_ = Gprime;
                        nG_ = nGprime;
                }

                for (i = 0; i < nG_; i++)
                        if (((double *)(y->x))[G_[i]] < -kkt_tolerance) {
                                H2[nH2++] = G_[i];
                        }

#ifndef NDEBUG
                {
                /* Paranoia: check for duplicate entries. */
                        int count;
                        for (i = 0; i < nH1; i++) {
                                count = 0;
                                for (j = 0; j < nH1; j++) {
                                        if (H1[j] == H1[i]) count++;
                                }
                                assert(count == 1);
                        }
                        for (i = 0; i < nH2; i++) {
                                count = 0;
                                for (j = 0; j < nH2; j++) {
                                        if (H2[j] == H2[i]) count++;
                                }
                                assert(count == 1);
                        }
                }
#endif

                /*
                 * If a coefficient was marked for constraint at the end of
                 * the inner loop and then freed again here, H1 and H2 can
                 * contain common elements. This is silly, since the change
                 * is a no-op. Make the sets disjoint again.
                 */
                if (nH1 > 0)
                        for (i = 0, j = 0; i < nH2; ) {
                                while((H1[j] < H2[i]) && (j < nH1)) j++;
                                if ((j < nH1) && (H2[i] == H1[j])) {
                                        /* Remove the element from both */
                                        for (k = j; k+1 < nH1; k++)
                                                H1[k] = H1[k+1];
                                        for (k = i; k+1 < nH2; k++)
                                                H2[k] = H2[k+1];
                                        --nH1;
                                        --nH2;
                                } else {
                                        ++i;
                                }
                        }

                y_min = 0;
                for (i = 0; i < nH2; i++)
                        if ((((double *)(y->x))[H2[i]] < y_min))
                                y_min = ((double *)(y->x))[H2[i]];

                /*
                 * If we've satisfied the KKT conditions, we're done. 
                 */

                if (nH2 == 0) break;

                ninf = nH1 + nH2;

                if (verbose)
                        printf("\tFreeing %ld coefficients (y_min = %e)\n",
                            nH2, y_min);

                if (verbose)
                        printf("Iteration %d Infeasibles: %ld Residual: %.20e\n",
                            iter, ninf, residual);

                feasible = false;
                while (!feasible) {

                /*
                 * Next, update the factorization of AtA[:,F][F,:], removing
                 * rows in H1 and adding those in H2. For large updates, the
                 * factorization will be recomputed from scratch; for small
                 * updates, the rows will be added and deleted from the
                 * factorization one-by-one, potentially saving a lot of
                 * computing time.
                 *
                 * F, G, H1, H2, and the associated counts are updated on
                 * exit.
                 */
                        L = modify_factor(AtA, L, F, &nF, G, &nG, H1, &nH1, 
                            H2, &nH2, verbose, c);

                        if (verbose) t0 = clock();

                        AtA_F = NULL;
                        Atb_F = x_F = NULL;

                        /*
                         * Solve the full system, but with the rows of L 
                         * corresponding to G set to identity. 
                         */
                        if (L->n == nvar) {
                                cholmod_dense *x_temp;

                                x_temp = cholmod_l_solve(CHOLMOD_A, L, Atb, c);
                                /*
                                 * NB: x_F must be the same length as F.
                                 * Extract the free set from the full solution.
                                 */
                                x_F = cholmod_l_allocate_dense(nF, 1, nF,
                                    CHOLMOD_REAL, c);
                                for (i = 0; i < nF; i++)
                                        ((double*)(x_F->x))[i] =
                                            ((double*)(x_temp->x))[F[i]];
                                cholmod_l_free_dense(&x_temp, c);
                        } else {
                                /*
                                 * If L is only a subset of the full matrix,
                                 * solve only the passive set.
                                 */
                                cholmod_dense *Atb_F;
        
                                Atb_F = cholmod_l_allocate_dense(nF, 1, nF,
                                    CHOLMOD_REAL, c);
                                for (i = 0; i < nF; i++)
                                        ((double*)(Atb_F->x))[i] =
                                            ((double*)(Atb->x))[F[i]];
        
                                x_F = cholmod_l_solve(CHOLMOD_A, L, Atb_F, c);
                        }

                        ++solves;

                        if (verbose) {
                                t1 = clock();
                                printf("\tSolve[%d] (%ld free): %.2f s\n",
                                    nvar,nF,(double)(t1-t0)/(CLOCKS_PER_SEC));
                        }

                        /*
                         * Characterize the new solution: how many of the new
                         * coefficients are negative, and of those, how many
                         * correspond to elements of the current solution close
                         * enough to zero that any movement along the descent
                         * vector will make those coefficients negative?
                         */
                        int nF_inf, nF_inf_boundary;
                        nF_inf = nF_inf_boundary = 0;
                        for (i = 0; i < nF; i++)
                                if (((double*)(x_F->x))[i] < 0) {
                                        nF_inf++;
                                        if (((double*)(x->x))[F[i]] < 
                                            kkt_tolerance)
                                                nF_inf_boundary++;
                                        
                                }

                        if (nF_inf == 0) {
                                /*
                                 * The new solution doesn't violate any
                                 * constraints. Accept it.
                                 */
                                for (i = 0; i < nF; i++)
                                        ((double*)(x->x))[F[i]] =
                                            ((double*)(x_F->x))[i];
                                cholmod_l_free_dense(&x_F, c);
                                feasible = true;

                                if (verbose)
                                        printf("\tSolution entirely "
                                            "feasible\n");

                        } else if (nF_inf == nF_inf_boundary) {
                                /*
                                 * Part of the new solution is negative, but
                                 * the corresponding coefficients in the
                                 * current solution are already nearly 0. Any
                                 * movement along the descent vector will
                                 * encounter a constraint, so we can't reduce
                                 * the residual with the current constraints.
                                 * Constrain the corresponding coefficients
                                 * and try again.
                                 */
                                for (i = 0; i < nF; i++) {
                                        if (((double*)(x_F->x))[i] < 0) {
                                                H1[nH1++] = F[i];
                                                ((double*)(x->x))[F[i]] = 0;
                                        }
                                }
                                cholmod_l_free_dense(&x_F, c);
                                feasible = false;

                                if (verbose)
                                        printf("\tConstraining %ld coefficients"
                                            " (descent at boundary)\n", nH1);
                                
                        } else {
                                /*
                                 * Attempt to reduce the magnitude of the
                                 * residual vector by moving along the descent
                                 * vector. Ideally, we'd like to minimize the
                                 * the residual along the descent vector
                                 * orthogonally projected into the feasible
                                 * region, but that way lies madness. Instead,
                                 * we choose from amongst the set of distances
                                 * along the descent vector that violate
                                 * exactly one constraint, trying each distance
                                 * in descending order until we reduce the 
                                 * residual with respect to the current
                                 * solution.
                                 */

                                assert( nH1 == 0 );
                                assert( nH2 == 0 );

                                /* AtA is symmetric, but stored in full form */
                                AtA->stype = 0;
                                /* FIXME: re-use AtA_F from modify call? */
                                AtA_F = cholmod_l_submatrix(AtA, F, nF, F, nF,
                                    1, 1, c);
                                AtA->stype = 1;

                                if (!Atb_F) {
                                        Atb_F = cholmod_l_allocate_dense(nF, 1,
                                             nF, CHOLMOD_REAL, c);
                                        for (i = 0; i < nF; i++) 
                                                ((double*)(Atb_F->x))[i] =
                                                    ((double*)(Atb->x))[F[i]];
                                }

                                feasible = walk_descents(AtA_F, Atb_F, x, x_F,
                                    F, &nF, H1, &nH1, &residual,
                                    &residual_calcs, verbose, c);

                        } /* if (nF_inf == 0) */

                        /*
                         * At this point, we've either rejected x_F completely
                         * or copied the feasible parts into x.
                         */
                        cholmod_l_free_dense(&x_F, c);
                        cholmod_l_free_dense(&Atb_F, c);
                        cholmod_l_free_sparse(&AtA_F, c);

                } /* while (!feasible) */

                if (nH1 == 0) {
                        /* 
                         * The passive set hasn't changed since the last
                         * factorization. 
                         */
                        nGprime = nFprime = -1;
                        F_ = F; nF_ = nF;
                        G_ = G; nG_ = nG;
                } else {
                        /*
                         * The passive set has changed since the last
                         * factorization. Copy F and G to Fprime and Gprime,
                         * and move elements in H1 from Fprime to Gprime.
                         */
                        memcpy(Fprime, F, nF*sizeof(long));
                        memcpy(Gprime, G, nG*sizeof(long));
                        F_ = Fprime; nF_ = nF;
                        G_ = Gprime; nG_ = nG;
                        /* Move elements from F_ into G_ */
                        for (i = 0, j = 0; i < nH1; i++) {
                                G_[nG_++] = H1[i];
                                while (F_[j] != H1[i]) j++;
                                for (k = j+i; k+1 < nF_; k++)
                                        F_[k-i] = F_[k-i+1];
                        }
                        /*
                         * NB: we do not reset nH1 here, since we still have
                         * to apply the changes to F before the next solve.
                         */
                        nF_ -= nH1;
                        nFprime = nF_;
                        nGprime = nG_;

                        assert( nFprime + nGprime == nvar );
                
                        qsort(G_, nG_, sizeof(G_[0]), intcmp);
                }

                /* Update the constrained part */
                {
                        cholmod_sparse *AtA_FG;
                        cholmod_dense *Atb_G, *x_F;
                
                        /*
                         * NB: We happen to know that AtA is actually
                         * symmetric, but cholmod_submatrix will only
                         * handle unsymmetric ones. The submatrix constructed
                         * this way will also be secretly symmetric.
                         */
                        AtA->stype = 0;
                        double ones[2] = {1., 0}, mones[2] = {-1., 0};

                        if (verbose) t0 = clock();

                        Atb_G = cholmod_l_allocate_dense(nG_, 1, nG_,
                            CHOLMOD_REAL, c);
                        for (i = 0; i < nG_; i++)
                                ((double *)(Atb_G->x))[i] =
                                    ((double *)(Atb->x))[G_[i]];

                        x_F = cholmod_l_allocate_dense(nF_, 1, nF_,
                            CHOLMOD_REAL, c);
                        for (i = 0; i < nF_; i++)
                                ((double *)(x_F->x))[i] =
                                    ((double *)(x->x))[F_[i]];

                        AtA_FG = cholmod_l_submatrix(AtA, G_, nG_, F_, nF_, 
                            1, 1, c);
                        cholmod_l_sdmult(AtA_FG, 0, ones, mones, x_F, Atb_G, c);
                        
                        for (i = 0; i < nG_; i++)
                                ((double *)(y->x))[G_[i]] =
                                    ((double *)(Atb_G->x))[i];

                        cholmod_l_free_dense(&Atb_G, c);
                        cholmod_l_free_dense(&x_F, c);
                        cholmod_l_free_sparse(&AtA_FG, c);

                        if (verbose) {
                                t1 = clock();
                                printf("\tUpdate y[%ld]: %.2f s\n",
                                    nG_,(double)(t1-t0)/(CLOCKS_PER_SEC));
                        }
                        AtA->stype = 1;
                        
                }

                /* Explicitly zero x[G] and y[F] */
                for (i = 0; i < nG_; i++)
                        ((double *)(x->x))[G_[i]] = 0;
                for (i = 0; i < nF_; i++)
                        ((double *)(y->x))[F_[i]] = 0;

#if 0
                double* x_full = (double*)(x->x);
                double* y_full = (double*)(y->x);
                DPRINT(x_full, nvar, " %- .1e");
                DPRINT(y_full, nvar, " %- .1e");
#endif
        }

        cholmod_l_free_dense(&y, c);
        cholmod_l_free_factor(&L, c);

        free(F);
        free(G);
        free(Fprime);
        free(Gprime);
        free(H1);
        free(H2);

        if (verbose) {
                if (iter == max_iter)
                        printf("VARNING! Failed to converge after ");
                else
                        printf("Finished in ");
                printf("%d iterations (%d factorizations, %d "
                    "residual calculations).\n", iter+1,solves,residual_calcs);
        }

        return (x);
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

