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
#define PHOTOSPLINE_NVECS PHOTOSPLINE_MAXDIM/PHOTOSPLINE_VECTOR_SIZE

namespace photospline { namespace detail {

template <typename Float>
struct simd_vector;

template <>
struct simd_vector<float> {
#if __GNUC__ == 3
#if PHOTOSPLINE_VECTOR_SIZE != 4
	#error On GCC 3, PHOTOSPLINE_VECTOR_SIZE must be 4!
#endif
	typedef float type __attribute__(( mode(V4SF) ));
#else
	typedef float type __attribute__((vector_size(PHOTOSPLINE_VECTOR_SIZE*sizeof(float))));
#endif
	static void init(type &a, float b)
	{
#if defined(__i386__) || defined (__x86_64__)
		a = _mm_set1_ps(b);
#elif defined(__powerpc__)
#ifdef vec_splats
		a = vec_splats(b);
#else
		float b_tmp __aligned(16) = b;
		a = vec_splat(*((v4sf *)(&b_tmp)), 0);
#endif
#else
		for (int i=0; i<PHOTOSPLINE_VECTOR_SIZE; i++)
			((float *)(&a))[i] = b;
#endif
	}
};

}}

#endif //PHOTOSPLINE_SIMD_H
