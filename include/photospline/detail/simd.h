#ifndef PHOTOSPLINE_SIMD_H
#define PHOTOSPLINE_SIMD_H

#if defined(__i386__) || defined (__x86_64__)
#ifdef __GLIBC__
#include <alloca.h>
#endif
#include <xmmintrin.h>
#elif defined(__powerpc__)
#include <altivec.h>
#endif

#include "photospline/bspline.h"

#define PHOTOSPLINE_MAXDIM	    8
#define PHOTOSPLINE_VECTOR_SIZE 4

#if __GNUC__ == 3
#if PHOTOSPLINE_VECTOR_SIZE != 4
	#error On GCC 3, PHOTOSPLINE_VECTOR_SIZE must be 4!
#endif
typedef float v4sf __attribute__(( mode(V4SF) ));
#else
typedef float v4sf __attribute__((vector_size(PHOTOSPLINE_VECTOR_SIZE*sizeof(float))));
#endif

#define PHOTOSPLINE_NVECS PHOTOSPLINE_MAXDIM/PHOTOSPLINE_VECTOR_SIZE

#if defined(__i386__) || defined (__x86_64__)
#define v4sf_init(a, b) a = _mm_set1_ps(b)
#elif defined(__powerpc__)
#ifdef vec_splats
#define v4sf_init(a, b) a = vec_splats(b)
#else
#define v4sf_init(a, b) { float b_tmp __aligned(16) = b; \
	a = vec_splat(*((v4sf *)(&b_tmp)), 0); }
#endif
#else
#define v4sf_init(a, b) { \
	((float *)(&a))[0] = b; \
	((float *)(&a))[1] = b; \
	((float *)(&a))[2] = b; \
	((float *)(&a))[3] = b; \
}
#endif

#endif //PHOTOSPLINE_SIMD_H
