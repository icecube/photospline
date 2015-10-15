#ifndef PHOTOSPLINE_C_SPLINETABLE_H
#define PHOTOSPLINE_C_SPLINETABLE_H

#include <stdint.h>
#include <stddef.h>

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

//TODO:
//	sampling?
//	fitting
//	FITS serialization/deserialization

struct splinetable_buffer {
	void* data;
	size_t size;
};
	
int readsplinefitstable_mem(const struct splinetable_buffer* buffer,
                            struct splinetable* table);

int writesplinefitstable_mem(struct splinetable_buffer* buffer,
                             const struct splinetable* table);
	
#ifdef __cplusplus
} //extern "C"
#endif

#endif