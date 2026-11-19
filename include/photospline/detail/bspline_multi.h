#ifndef PHOTOSPLINE_BSPLINE_MULTI_H
#define PHOTOSPLINE_BSPLINE_MULTI_H

#include "photospline/detail/simd.h"

namespace photospline{

template <typename Alloc>
template <typename Float, unsigned int Row, unsigned int NV>
void splinetable<Alloc>::ndsplineeval_multibasis_core_row(const int *centers, const typename detail::simd_vector<Float>::type*** localbasis, typename detail::simd_vector<Float>::type* result) const{
#if defined(__i386__) && defined(__ELF__)
	/*
	 * Work around GCC ABI-compliance issue with SSE on x86 by
	 * forcibly realigning the stack to a 16-byte boundary.
	 */
	volatile register unsigned long sp __asm("esp");
	__asm("" : "=r"(sp));
	if (__builtin_expect(sp & 15UL, 0))
		(void)alloca(16 - (sp & 15UL));
#endif
	detail::contract<Row, NV>(&coefficients[0], &strides[0], &order[0], ndim, centers,
	    [&](unsigned d, unsigned i){ return localbasis[d][i]; }, result);
}

template <typename Alloc>
template <typename Float>
typename splinetable<Alloc>::template multibasis_core_ptr<Float> splinetable<Alloc>::select_multibasis_core() const{
	// value plus ndim derivatives, PHOTOSPLINE_VECTOR_SIZE lanes per vector
	const unsigned nv = (ndim + PHOTOSPLINE_VECTOR_SIZE) / PHOTOSPLINE_VECTOR_SIZE;
#ifndef PHOTOSPLINE_NO_EVAL_TEMPLATES
	static_assert(PHOTOSPLINE_NVECS == 2, "variants below assume at most two vectors");
	switch ((order[ndim-1] + 1)*10 + nv) {
		case 21: return &splinetable::ndsplineeval_multibasis_core_row<Float, 2, 1>;
		case 22: return &splinetable::ndsplineeval_multibasis_core_row<Float, 2, 2>;
		case 31: return &splinetable::ndsplineeval_multibasis_core_row<Float, 3, 1>;
		case 32: return &splinetable::ndsplineeval_multibasis_core_row<Float, 3, 2>;
		case 41: return &splinetable::ndsplineeval_multibasis_core_row<Float, 4, 1>;
		case 42: return &splinetable::ndsplineeval_multibasis_core_row<Float, 4, 2>;
		case 51: return &splinetable::ndsplineeval_multibasis_core_row<Float, 5, 1>;
		case 52: return &splinetable::ndsplineeval_multibasis_core_row<Float, 5, 2>;
		case 61: return &splinetable::ndsplineeval_multibasis_core_row<Float, 6, 1>;
		case 62: return &splinetable::ndsplineeval_multibasis_core_row<Float, 6, 2>;
	}
#endif
	(void)nv;
	return &splinetable::ndsplineeval_multibasis_core_row<Float, 0, PHOTOSPLINE_NVECS>;
}

template <typename Alloc>
template <typename Float>
void splinetable<Alloc>::ndsplineeval_multibasis_core(const int *centers, const typename detail::simd_vector<Float>::type*** localbasis, typename detail::simd_vector<Float>::type* result) const{
	(this->*select_multibasis_core<Float>())(centers, localbasis, result);
}

/* Evaluate the spline surface and all its derivatives at x */

template<typename Alloc>
template<typename Float>
void
splinetable<Alloc>::ndsplineeval_gradient(const double* x, const int* centers, double* evaluates) const
{
	uint32_t maxdegree = *std::max_element(order,order+ndim) + 1;
	uint32_t nbases = ndim + 1;
	typename detail::simd_vector<Float>::type acc[PHOTOSPLINE_NVECS];
	Float valbasis[maxdegree];
	Float gradbasis[maxdegree];
	typename detail::simd_vector<Float>::type localbasis[ndim][maxdegree][PHOTOSPLINE_NVECS];
	const typename detail::simd_vector<Float>::type* localbasis_rowptr[ndim][maxdegree];
	const typename detail::simd_vector<Float>::type** localbasis_ptr[ndim];

	assert(ndim > 0);
	if (ndim+1 > PHOTOSPLINE_MAXDIM)
		throw std::runtime_error("Error: ndsplineeval_gradient() can only "
		    "process up to "+std::to_string(PHOTOSPLINE_MAXDIM-1)+"-dimensional tables. "
		    "Adjust PHOTOSPLINE_MAXDIM in detail/simd.h to change this.");
		
	for (uint32_t n = 0; n < ndim; n++) {

		/* 
		 * Compute the values and derivatives of the table->order[n]+1 non-zero
		 * splines at x[n], filling them into valbasis and gradbasis.
		 */
		bspline_nonzero(&*knots[n], nknots[n],
		    x[n], centers[n], order[n], valbasis, gradbasis);
	
		for (uint32_t i = 0; i <= order[n]; i++) {
			
			((Float*)(localbasis[n][i]))[0] = valbasis[i];
			
			for (uint32_t j = 1; j < ndim+1; j++) {
				if (j == 1+n)
					((Float*)(localbasis[n][i]))[j] = gradbasis[i];
				else
					((Float*)(localbasis[n][i]))[j] = valbasis[i];
			}
			
			localbasis_rowptr[n][i] = localbasis[n][i];
		}
		
		localbasis_ptr[n] = localbasis_rowptr[n];
	}

	Float* acc_ptr = (Float*)acc;

	for (uint32_t i = 0; i < nbases; i++)
		acc_ptr[i] = 0;

	ndsplineeval_multibasis_core<Float>(centers, localbasis_ptr, acc);

	for (uint32_t i = 0; i < nbases; i++)
		evaluates[i] = acc_ptr[i];
}

template<typename Alloc>
template<typename Float>
void splinetable<Alloc>::evaluator_type<Float>::ndsplineeval_gradient(const double* x, const int* centers, double* evaluates) const
{
	uint32_t maxdegree = *std::max_element(table.order,table.order+table.ndim) + 1;
	uint32_t nbases = table.ndim + 1;
	typename detail::simd_vector<Float>::type acc[PHOTOSPLINE_NVECS];
	Float valbasis[maxdegree];
	Float gradbasis[maxdegree];
	typename detail::simd_vector<Float>::type localbasis[table.ndim][maxdegree][PHOTOSPLINE_NVECS];
	const typename detail::simd_vector<Float>::type* localbasis_rowptr[table.ndim][maxdegree];
	const typename detail::simd_vector<Float>::type** localbasis_ptr[table.ndim];

	assert(table.ndim > 0);
	if (table.ndim+1 > PHOTOSPLINE_MAXDIM)
		throw std::runtime_error("Error: ndsplineeval_gradient() can only "
		    "process up to "+std::to_string(PHOTOSPLINE_MAXDIM-1)+"-dimensional tables. "
		    "Adjust PHOTOSPLINE_MAXDIM in detail/simd.h to change this.");
		
	for (uint32_t n = 0; n < table.ndim; n++) {

		/* 
		 * Compute the values and derivatives of the table->order[n]+1 non-zero
		 * splines at x[n], filling them into valbasis and gradbasis.
		 */
		bspline_nonzero(&*table.knots[n], table.nknots[n],
		    x[n], centers[n], table.order[n], valbasis, gradbasis);
	
		for (uint32_t i = 0; i <= table.order[n]; i++) {
			
			((Float*)(localbasis[n][i]))[0] = valbasis[i];
			
			for (uint32_t j = 1; j < table.ndim+1; j++) {
				if (j == 1+n)
					((Float*)(localbasis[n][i]))[j] = gradbasis[i];
				else
					((Float*)(localbasis[n][i]))[j] = valbasis[i];
			}
			
			localbasis_rowptr[n][i] = localbasis[n][i];
		}
		
		localbasis_ptr[n] = localbasis_rowptr[n];
	}

	Float* acc_ptr = (Float*)acc;

	for (uint32_t i = 0; i < nbases; i++)
		acc_ptr[i] = 0;

	(table.*(v_eval_ptr))(centers, localbasis_ptr, acc);

	for (uint32_t i = 0; i < nbases; i++)
		evaluates[i] = acc_ptr[i];
}

} //namespace photospline

#endif
