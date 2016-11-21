#ifndef PHOTOSPLINE_BSPLINE_MULTI_H
#define PHOTOSPLINE_BSPLINE_MULTI_H

#include "photospline/detail/simd.h"

namespace photospline{
	
void
bspline_nonzero(const double* knots, const unsigned nknots,
                const double x, int left, const int n,
                float* values, float* derivs);
	
template <typename Alloc>
void splinetable<Alloc>::ndsplineeval_multibasis_core(const int *centers, const v4sf*** localbasis, v4sf* result) const{
#if (defined(__i386__) || defined (__x86_64__)) && defined(__ELF__)
	/*
	 * Work around GCC ABI-compliance issue with SSE on x86 by
	 * forcibly realigning the stack to a 16-byte boundary.
	 */
	volatile register unsigned long sp __asm("esp");
	__asm("" : "=r"(sp));
	if (__builtin_expect(sp & 15UL, 0))
		(void)alloca(16 - (sp & 15UL));
#endif
	v4sf basis_tree[ndim+1][PHOTOSPLINE_NVECS];
	int decomposedposition[ndim];
	
	int64_t tablepos = 0;
	for (uint32_t n = 0; n < ndim; n++) {
		decomposedposition[n] = 0;
		tablepos += (centers[n] - order[n])*strides[n];
	}
	
	for (uint32_t k = 0; k < PHOTOSPLINE_NVECS; k++) {
		v4sf_init(basis_tree[0][k], 1);
		for (uint32_t n = 0; n < ndim; n++)
			basis_tree[n+1][k] = basis_tree[n][k]*localbasis[n][0][k];
	}
	
	uint32_t nchunks = 1;
	for (uint32_t n = 0; n < ndim - 1; n++)
		nchunks *= (order[n] + 1);
	
	uint32_t n = 0;
	while (1) {
		for (uint32_t i = 0; __builtin_expect(i < order[ndim-1] + 1, 1); i++) {
			v4sf weights;
			v4sf_init(weights, coefficients[tablepos + i]);
			for (uint32_t k = 0; k < PHOTOSPLINE_NVECS; k++)
				result[k] += basis_tree[ndim-1][k]*
				localbasis[ndim-1][i][k]*weights;
		}
		
		if (__builtin_expect(++n == nchunks, 0))
			break;
		
		tablepos += strides[ndim-2];
		decomposedposition[ndim-2]++;
		
		/* Carry to higher dimensions */
		uint32_t i;
		for (i = ndim-2;
			 decomposedposition[i] > order[i]; i--) {
			decomposedposition[i-1]++;
			tablepos += (strides[i-1]
						 - decomposedposition[i]*strides[i]);
			decomposedposition[i] = 0;
		}
		for (uint32_t j = i; __builtin_expect(j < ndim-1, 1); j++)
			for (uint32_t k = 0; k < PHOTOSPLINE_NVECS; k++)
				basis_tree[j+1][k] = basis_tree[j][k]*
				localbasis[j][decomposedposition[j]][k];
	}
}
	
namespace{
	template <unsigned int D>
	struct vectorCountHelper{
		static constexpr unsigned int VC =
		((D+1) / PHOTOSPLINE_VECTOR_SIZE)
		+ ((D+1) % PHOTOSPLINE_VECTOR_SIZE ? 1 : 0);
	};
}
	
template <typename Alloc>
template <unsigned int D>
void splinetable<Alloc>::ndsplineeval_multibasis_coreD(const int *centers, const v4sf*** localbasis, v4sf* result) const{
#if (defined(__i386__) || defined (__x86_64__)) && defined(__ELF__)
	/*
	 * Work around GCC ABI-compliance issue with SSE on x86 by
	 * forcibly realigning the stack to a 16-byte boundary.
	 */
	volatile register unsigned long sp __asm("esp");
	if (__builtin_expect(sp & 15UL, 0))
		(void)alloca(16 - (sp & 15UL));
#endif
	const unsigned int VC=vectorCountHelper<D>::VC;
	v4sf basis_tree[D+1][VC];
	int decomposedposition[D];
	
	int64_t tablepos = 0;
	for (uint32_t n = 0; n < D; n++) {
		decomposedposition[n] = 0;
		tablepos += (centers[n] - order[n])*strides[n];
	}
	
	for (uint32_t k = 0; k < VC; k++) {
		v4sf_init(basis_tree[0][k], 1);
		for (uint32_t n = 0; n < D; n++)
			basis_tree[n+1][k] = basis_tree[n][k]*localbasis[n][0][k];
	}
	
	uint32_t nchunks = 1;
	for (uint32_t n = 0; n < D - 1; n++)
		nchunks *= (order[n] + 1);
	
	uint32_t n = 0;
	while (1) {
		for (uint32_t i = 0; __builtin_expect(i < order[D-1] + 1, 1); i++) {
			v4sf weights;
			v4sf_init(weights, coefficients[tablepos + i]);
			for (uint32_t k = 0; k < VC; k++)
				result[k] += basis_tree[D-1][k]*localbasis[D-1][i][k]*weights;
		}
		
		if (__builtin_expect(++n == nchunks, 0))
			break;
		
		tablepos += strides[D-2];
#ifdef __clang__ //this code is unreachable if D<2
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Warray-bounds"
#endif
		decomposedposition[D-2]++;
#ifdef __clang__
	#pragma clang diagnostic pop
#endif
		
		/* Carry to higher dimensions */
		uint32_t i;
		for (i = D-2; decomposedposition[i] > order[i]; i--) {
			decomposedposition[i-1]++;
			tablepos += (strides[i-1] - decomposedposition[i]*strides[i]);
			decomposedposition[i] = 0;
		}
		for (uint32_t j = i; __builtin_expect(j < D-1, 1); j++)
			for (uint32_t k = 0; k < VC; k++)
				basis_tree[j+1][k] = basis_tree[j][k]*
				localbasis[j][decomposedposition[j]][k];
	}
}

template <typename Alloc>
template <unsigned int D, unsigned int Order>
void splinetable<Alloc>::ndsplineeval_multibasis_coreD_FixedOrder(const int *centers, const v4sf*** localbasis, v4sf* result) const{
#if (defined(__i386__) || defined (__x86_64__)) && defined(__ELF__)
	/*
	 * Work around GCC ABI-compliance issue with SSE on x86 by
	 * forcibly realigning the stack to a 16-byte boundary.
	 */
	volatile register unsigned long sp __asm("esp");
	if (__builtin_expect(sp & 15UL, 0))
		(void)alloca(16 - (sp & 15UL));
#endif
	const unsigned int VC=vectorCountHelper<D>::VC;
	v4sf basis_tree[D+1][VC];
	int decomposedposition[D];
	
	int64_t tablepos = 0;
	for (uint32_t n = 0; n < D; n++) {
		decomposedposition[n] = 0;
		tablepos += (centers[n] - Order)*strides[n];
	}
	
	for (uint32_t k = 0; k < VC; k++) {
		v4sf_init(basis_tree[0][k], 1);
		for (uint32_t n = 0; n < D; n++)
			basis_tree[n+1][k] = basis_tree[n][k]*localbasis[n][0][k];
	}
	
	uint32_t nchunks = 1;
	for (uint32_t n = 0; n < D - 1; n++)
		nchunks *= (Order + 1);
	
	uint32_t n = 0;
	while (1) {
		for (uint32_t i = 0; __builtin_expect(i < Order + 1, 1); i++) {
			v4sf weights;
			v4sf_init(weights, coefficients[tablepos + i]);
			for (uint32_t k = 0; k < VC; k++)
				result[k] += basis_tree[D-1][k]*
				localbasis[D-1][i][k]*weights;
		}
		
		if (__builtin_expect(++n == nchunks, 0))
			break;
		
		tablepos += strides[D-2];
#ifdef __clang__ //this code is unreachable if D<2
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Warray-bounds"
#endif
		decomposedposition[D-2]++;
#ifdef __clang__
	#pragma clang diagnostic pop
#endif
		
		/* Carry to higher dimensions */
		uint32_t i;
		for (i = D-2;
			 decomposedposition[i] > Order; i--) {
			decomposedposition[i-1]++;
			tablepos += (strides[i-1]
						 - decomposedposition[i]*strides[i]);
			decomposedposition[i] = 0;
		}
		for (uint32_t j = i; __builtin_expect(j < D-1, 1); j++)
			for (uint32_t k = 0; k < VC; k++)
				basis_tree[j+1][k] = basis_tree[j][k]*
				localbasis[j][decomposedposition[j]][k];
	}
}
	
/* Evaluate the spline surface and all its derivatives at x */

template<typename Alloc>
void
splinetable<Alloc>::ndsplineeval_gradient(const double* x, const int* centers, double* evaluates) const
{
	uint32_t maxdegree = *std::max_element(order,order+ndim) + 1;
	uint32_t nbases = ndim + 1;
	v4sf acc[PHOTOSPLINE_NVECS];
	float valbasis[maxdegree];
	float gradbasis[maxdegree];
	v4sf localbasis[ndim][maxdegree][PHOTOSPLINE_NVECS];
	const v4sf* localbasis_rowptr[ndim][maxdegree];
	const v4sf** localbasis_ptr[ndim];

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
		bspline_nonzero(knots[n], nknots[n],
		    x[n], centers[n], order[n], valbasis, gradbasis);
	
		for (uint32_t i = 0; i <= order[n]; i++) {
			
			((float*)(localbasis[n][i]))[0] = valbasis[i];
			
			for (uint32_t j = 1; j < ndim+1; j++) {
				if (j == 1+n)
					((float*)(localbasis[n][i]))[j] = gradbasis[i];
				else
					((float*)(localbasis[n][i]))[j] = valbasis[i];
			}
			
			localbasis_rowptr[n][i] = localbasis[n][i];
		}
		
		localbasis_ptr[n] = localbasis_rowptr[n];
	}

	float* acc_ptr = (float*)acc;

	for (uint32_t i = 0; i < nbases; i++)
		acc_ptr[i] = 0;

	ndsplineeval_multibasis_core(centers, localbasis_ptr, acc);

	for (uint32_t i = 0; i < nbases; i++)
		evaluates[i] = acc_ptr[i];
}

template<typename Alloc>
void splinetable<Alloc>::evaluator::ndsplineeval_gradient(const double* x, const int* centers, double* evaluates) const{
	uint32_t maxdegree = *std::max_element(table.order,table.order+table.ndim) + 1;
	uint32_t nbases = table.ndim + 1;
	v4sf acc[PHOTOSPLINE_NVECS];
	float valbasis[maxdegree];
	float gradbasis[maxdegree];
	v4sf localbasis[table.ndim][maxdegree][PHOTOSPLINE_NVECS];
	const v4sf* localbasis_rowptr[table.ndim][maxdegree];
	const v4sf** localbasis_ptr[table.ndim];
	
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
		bspline_nonzero(table.knots[n], table.nknots[n],
						x[n], centers[n], table.order[n], valbasis, gradbasis);
		
		for (uint32_t i = 0; i <= table.order[n]; i++) {
			
			((float*)(localbasis[n][i]))[0] = valbasis[i];
			
			for (uint32_t j = 1; j < table.ndim+1; j++) {
				if (j == 1+n)
					((float*)(localbasis[n][i]))[j] = gradbasis[i];
				else
					((float*)(localbasis[n][i]))[j] = valbasis[i];
			}
			
			localbasis_rowptr[n][i] = localbasis[n][i];
		}
		
		localbasis_ptr[n] = localbasis_rowptr[n];
	}
	
	float* acc_ptr = (float*)acc;
	
	for (uint32_t i = 0; i < nbases; i++)
		acc_ptr[i] = 0;
	
	(table.*(v_eval_ptr))(centers, localbasis_ptr, acc);
	
	for (uint32_t i = 0; i < nbases; i++)
		evaluates[i] = acc_ptr[i];
}

} //namespace photospline

#endif
