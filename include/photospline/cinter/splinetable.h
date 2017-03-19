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

/// Construct a new spline table
int splinetable_init(struct splinetable* table);

/// Destroy spline table
void splinetable_free(struct splinetable* table);

/// Read a table from a FITS file on disk
int readsplinefitstable(const char* path, struct splinetable* table);

/// Write the table to a FITS file on disk
int writesplinefitstable(const char* path, const struct splinetable* table);

struct splinetable_buffer {
	void* data;
	size_t size;
};

/// Read a table from a FITS 'file' in memory
int readsplinefitstable_mem(const struct splinetable_buffer* buffer,
                            struct splinetable* table);

/// Write the table to a FITS 'file' in memory
/// The destination buffer should be freed with free().
int writesplinefitstable_mem(struct splinetable_buffer* buffer,
                             const struct splinetable* table);

/// Get the specified FITS header key as a C string
const char* splinetable_get_key(const struct splinetable* table, const char* key);

/// Read a FITS header key
/// \param[in] type the type of the value
/// \param[in] key the key
/// \param[out] result where to store the result
int splinetable_read_key(const struct splinetable* table, splinetable_dtype type,
                         const char* key, void* result);
	
int splinetable_write_key(struct splinetable* table, splinetable_dtype type,
                          const char* key, const void* value);

// Access to spline properties

/// Number of dimensions
uint32_t splinetable_ndim(const struct splinetable* table);
/// Spline order in dimension \p dim
uint32_t splinetable_order(const struct splinetable* table, uint32_t dim);
/// Number of knots in dimension \p dim
uint64_t splinetable_nknots(const struct splinetable* table, uint32_t dim);
/// Knot vector in dimension \p dim
const double* splinetable_knots(const struct splinetable* table, uint32_t dim);
/// Knot \p knot in dimenson \p dim
double splinetable_knot(const struct splinetable* table, uint32_t dim,
                        uint64_t knot);
/// Lower extent of the spline's support in dimension \p dim
double splinetable_lower_extent(const struct splinetable* table, uint32_t dim);
/// Upper extent of the spline's support in dimension \p dim
double splinetable_upper_extent(const struct splinetable* table, uint32_t dim);
/// Period of the spline in dimension \p dim (currently unused)
double splinetable_period(const struct splinetable* table, uint32_t dim);
/// Number of splines along dimension \p dim
uint64_t splinetable_ncoeffs(const struct splinetable* table, uint32_t dim);
/// Total size of coefficient array
uint64_t splinetable_total_ncoeffs(const struct splinetable* table);
/// Stride of coefficient array in dimension \p dim
uint64_t splinetable_stride(const struct splinetable* table, uint32_t dim);
/// Coefficient array
const float* splinetable_coefficients(const struct splinetable* table);

// Spline evaluation

/// Find the indices of the central splines at coordinates \p x
///
/// \param[in] x coordinates
/// \param[out] centers indices of central splines
/// \returns 0 if \p is within the partial support of the spline
int tablesearchcenters(const struct splinetable* table, const double* x,
                       int* centers);

/// Evaluate the spline surface at coordinates \p x, with optional single differentiation
/// 
/// \param[in] x           coordinates at which to evaluate
/// \param[in] centers     indices of the central splines
/// \param[in] derivatives bitmask indicating which dimensions to differentiate
/// 
/// \pre \p centers has been filled via a call to tablesearchcenters()
/// \returns the value of the spline surface (or one of its derivatives) at \p x
double ndsplineeval(const struct splinetable* table, const double* x,
                    const int* centers, int derivatives);

/// Evaluate the spline surface and its gradient at coordinates \p x
/// 
/// \param[in] x         coordinates at which to evaluate
/// \param[in] centers   indices of the central splines
/// \param[in] evaluates storage for the evaluate and gradient (must have length at least ndim+1) 
/// 
/// \pre \p centers has been filled via a call to tablesearchcenters()
void ndsplineeval_gradient(const struct splinetable* table, const double* x,
                           const int* centers, double* evaluates);

/// Evaluate the spline surface at coordinates \p x, with optional double differentiation
/// 
/// \param[in] x           coordinates at which to evaluate
/// \param[in] centers     indices of the central splines
/// \param[in] derivatives order of derivative to calculate in each dimension
/// 
/// \pre \p centers has been filled via a call to tablesearchcenters()
/// \returns the value of the spline surface (or one of its second derivatives) at \p x
double ndsplineeval_deriv(const struct splinetable* table, const double* x,
                           const int* centers, const unsigned int *derivatives);

/// Convolve with another spline
///
/// \param[in] dim   dimension to convolve
/// \param[in] knots knot vector of convolution kernel
/// \param[in] n_knots size of kernel knot vector
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

/// Reorder the dimensions of the spline
///
/// \param[in] permutation an array of indices giving the new order of the dimensions
int splinetable_permute(struct splinetable* table, size_t* permutation);
	
#ifdef __cplusplus
} //extern "C"
#endif

#endif