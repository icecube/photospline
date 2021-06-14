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
struct simd_vector {
#ifdef __clang__
	typedef Float type __attribute__((ext_vector_type(PHOTOSPLINE_VECTOR_SIZE)));
#else
	typedef Float type __attribute__((vector_size(PHOTOSPLINE_VECTOR_SIZE*sizeof(Float))));
#endif

	static void init(type &a, Float b)
	{
		a = b - type{};
	}
};

}}

#endif //PHOTOSPLINE_SIMD_H
