
#ifndef PHOTOSPLINE_CFITTER_CHOLESKY_SOLVE_H
#define PHOTOSPLINE_CFITTER_CHOLESKY_SOLVE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <cholmod.h>
#include <stdbool.h>
#include <pthread.h>

cholmod_dense *
cholesky_solve(cholmod_sparse *AtA, cholmod_dense *Atb, cholmod_common *c, int verbose, int n_resolves);

cholmod_sparse *
get_column(cholmod_sparse *A, long k, long *iPerm, 
    long *Fset, long nF, cholmod_common *c);

cholmod_sparse *
submatrix_symm(cholmod_sparse *A, long *rows, unsigned long nrows,
    long *cols, unsigned long ncols, cholmod_common *c);

cholmod_factor * 
modify_factor(cholmod_sparse* A, cholmod_factor *L,
    long *F, long *nF, long *G, long *nG, long *H1, long *nH1,
    long *H2, long *nH2, int verbose, cholmod_common *c);

cholmod_factor * 
modify_factor_p(cholmod_sparse *A, cholmod_factor *L,
    long *F, long *nF_, long *G, long *nG_, long *H1, long *nH1_,
    long *H2, long *nH2_, bool update, bool verbose, cholmod_common *c);

cholmod_factor * 
recompute_factor(cholmod_sparse *A, cholmod_factor *L, long *iPerm,
    long *F, unsigned long nF, cholmod_common *c);

double
calc_residual(cholmod_sparse *AtA, cholmod_dense *Atb, cholmod_dense *x,
    cholmod_common *c);

#define FACTOR_INFO(L) printf(#L " is_ll: %d, is_super: %d, is_monotonic: %d, xtype: %d, ordering: %d\n",L->is_ll,L->is_super,L->is_monotonic,L->xtype,L->ordering)

#define DPRINT(vec,size,format)\
	{\
	int i;\
	printf(#vec":");\
	for (i = 0; i < size; i++) printf(format,vec[i]);\
	printf("\n");\
	}

/*
 * Get the number of threads to use from the environment variables
 * GOTO_NUM_THREADS or OMP_NUM_THREADS. Defaults to 1.
 */
int
get_nthreads(void);

typedef struct {
	
	cholmod_dense *x;   /* current solution */
	cholmod_dense *x_F; /* unconstrained local minimum */

	cholmod_sparse *AtA_F;
	cholmod_dense *Atb_F;
	cholmod_common *c; /* Caution to the wind. */

	const long *F;
	long nF;

	const double *alpha;   /* Descent scale to use */

	cholmod_dense *x_c; /* scaled descent, projected into feasible space */
	double residual; /* Resulting residual after orthogonal projection */

	long *H1; /* Infeasible coefficients for this scaled descent */
	long nH1;

	pthread_cond_t *cv;
	pthread_mutex_t *mutex;
	int state;
	int id;

} descent_trial;

enum { WAIT, RUN, TERMINATE } worker_thread_state;

void
evaluate_descent(void *trial_);

int
walk_descents(cholmod_sparse *AtA_F, 
    cholmod_dense *Atb_F, cholmod_dense *x, 
    cholmod_dense *x_F, long *F, long *nF_, long *H1, long *nH1_,
    double *residual, int *residual_calcs, int verbose, cholmod_common *c);

#ifdef __cplusplus
}
#endif

#endif /* PHOTOSPLINE_CFITTER_CHOLESKY_SOLVE_H */

