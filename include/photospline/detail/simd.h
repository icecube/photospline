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

#if defined(__i386__) || defined(__amd64__)
	#define PHOTOSPLINE_VECTOR_ISN_VARIANTS "avx512f","avx2","avx","sse4.2","default"
	#ifdef __clang__ //clang, obviously
		// this feature exists and works nicely in clang 14, except that
		// "multiversioned functions do not yet support function templates"
		// which is most of the places we want to use this
		#define PHOTOSPLINE_USE_TARGET_CLONING 0
	#elif defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) //gcc
		//Activating this feature causes gcc 8.3 to crash, but 8.5 works
		#if __GNUC__ >= 9 || (__GNUC__ == 8 && __GNUC_MINOR__ > 3)
			#define PHOTOSPLINE_USE_TARGET_CLONING 1
		#else
			#define PHOTOSPLINE_USE_TARGET_CLONING 0
		#endif
	#else
		//for other compilers, assume we don't have this
		#define PHOTOSPLINE_USE_TARGET_CLONING 0
	#endif
#else
	//For other architectures, leave this alone for now
	#define PHOTOSPLINE_USE_TARGET_CLONING 0
#endif

#if PHOTOSPLINE_USE_TARGET_CLONING
	#define PHOTOSPLINE_TARGET_CLONE __attribute__ ((target_clones( PHOTOSPLINE_VECTOR_ISN_VARIANTS )))
#else
	//make this a no-op
	#define PHOTOSPLINE_TARGET_CLONE
#endif

#endif //PHOTOSPLINE_SIMD_H
