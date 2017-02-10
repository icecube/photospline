#ifndef PHOTOSPLINE_C_SPLINETABLE_H
#define PHOTOSPLINE_C_SPLINETABLE_H

#include <stdint.h>
#include <stddef.h>

#ifdef PHOTOSPLINE_INCLUDES_SPGLAM
#include "photospline/detail/splineutil.h"
#endif //PHOTOSPLINE_INCLUDES_SPGLAM

#ifdef __cplusplus
extern "C" {
#endif

struct splinetable{
	void* data;
};
	
// Initialization and I/O
	
typedef enum {
	SPLINETABLE_INT,
	SPLINETABLE_DOUBLE
} splinetable_dtype;

int splinetable_init(struct splinetable* table);

void splinetable_free(struct splinetable* table);
	
int readsplinefitstable(const char* path, struct splinetable* table);

int writesplinefitstable(const char* path, const struct splinetable* table);

struct splinetable_buffer {
	void* data;
	size_t size;
};
	
int readsplinefitstable_mem(const struct splinetable_buffer* buffer,
                            struct splinetable* table);

int writesplinefitstable_mem(struct splinetable_buffer* buffer,
                             const struct splinetable* table);

const char* splinetable_get_key(const struct splinetable* table, const char* key);

int splinetable_read_key(const struct splinetable* table, splinetable_dtype type,
                         const char* key, void* result);
	
int splinetable_write_key(struct splinetable* table, splinetable_dtype type,
                          const char* key, const void* value);

// Access to spline properties

uint32_t splinetable_ndim(const struct splinetable* table);
uint32_t splinetable_order(const struct splinetable* table, uint32_t dim);
uint64_t splinetable_nknots(const struct splinetable* table, uint32_t dim);
const double* splinetable_knots(const struct splinetable* table, uint32_t dim);
double splinetable_knot(const struct splinetable* table, uint32_t dim,
                        uint64_t knot);
double splinetable_lower_extent(const struct splinetable* table, uint32_t dim);
double splinetable_upper_extent(const struct splinetable* table, uint32_t dim);
double splinetable_period(const struct splinetable* table, uint32_t dim);
uint64_t splinetable_ncoeffs(const struct splinetable* table, uint32_t dim);
uint64_t splinetable_total_ncoeffs(const struct splinetable* table);
uint64_t splinetable_stride(const struct splinetable* table, uint32_t dim);
const float* splinetable_coefficients(const struct splinetable* table);

// Spline evaluation
	
int tablesearchcenters(const struct splinetable* table, const double* x,
                       int* centers);
	
double ndsplineeval(const struct splinetable* table, const double* x,
                    const int* centers, int derivatives);

void ndsplineeval_gradient(const struct splinetable* table, const double* x,
                           const int* centers, double* evaluates);
	
double ndsplineeval_deriv2(const struct splinetable* table, const double* x,
                           const int* centers, int derivatives);

// Convolution with another spline
	
int splinetable_convolve(struct splinetable* table, const int dim,
                         const double* knots, size_t n_knots);

#ifdef PHOTOSPLINE_INCLUDES_SPGLAM
// Fitting to data

//unlike C++ interface smoothing and penaltyOrder must point to
//arrays of exactly the same sizes as data->ndim.
int splinetable_glamfit(struct splinetable* table, const struct ndsparse* data,
                        const double* weights, const double* const* coords,
                        const uint32_t* splineOrder, const double* const* knots,
                        const uint64_t* nknots,
                        const double* smoothing, const uint32_t* penaltyOrder,
                        uint32_t monodim=-1, bool verbose=true);
	
/// Evaluate a spline on a grid
///\param table the spline to evaluate
///\param coords an array of arrays of coordinates at which to evaluate in each
///       of the spline's dimensions. The length of this array must match the
///       spline's dimension.
///\param ncoords an array whose entries specify the number of entries in each
///       of the per-dimension arrays of coords. The length of this array must
///       also match the spline's dimension.
///\param A pointer to an ndsparse* which will be updated by this function to
///       point to a new ndsaprse, containing the results, which the caller is
///       responsible for destroying by passing it to ndsparse_destroy, unless
///       a failure occurs, in which case it will be set to NULL.
///\return 0 on success, non-zero on failure.
int splinetable_grideval(struct splinetable* table, const double* const* coords,
                         const uint32_t* ncoords, struct ndsparse** result);

///Deallocate an ndsparse allocated by splinetable_grideval
void ndsparse_destroy(struct ndsparse* nd);
#endif //PHOTOSPLINE_INCLUDES_SPGLAM

//TODO:
//	sampling?
	
int splinetable_permute(struct splinetable* table, size_t* permutation);
	
#ifdef __cplusplus
} //extern "C"
#endif

#endif