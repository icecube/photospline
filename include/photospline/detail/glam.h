/*
 * GLAM (Generalized Linear Array Model) is an algorithm for efficiently
 *   computing N-D least squares fits on a grid. More information can be
 *   found at http://www.ma.hw.ac.uk/~iain/research/GLAM.html
 */

#ifndef PHOTOSPLINE_CFITTER_GLAM_H
#define PHOTOSPLINE_CFITTER_GLAM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <photospline/detail/splineutil.h>

#include <cholmod.h>
	
#define PHOTOSPLINE_GLAM_NO_MONODIM ((uint32_t)-1)
   
///\return 0 on success, otherwise nonzero
int glamfit_complex(
	/*input data*/
	const struct ndsparse* data, const double* weights, const double* const* coords,
	/*spline table components*/
	uint32_t ndim, const uint64_t* nknots, const double* const* knots, const uint64_t* naxes,
	float* out_coefficients,
	/*fit properties*/
    const uint32_t* order, cholmod_sparse* penalty,
    uint32_t monodim, int verbose, cholmod_common *c);

cholmod_sparse* add_penalty_term(uint64_t* nsplines, double* knots, uint32_t ndim,
    uint32_t dim, uint32_t order, uint32_t porder, double scale, int mono,
    cholmod_sparse* penalty, cholmod_common* c);

#ifdef __cplusplus
}
#endif

#endif /* PHOTOSPLINE_CFITTER_GLAM_H */

