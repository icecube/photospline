#ifndef PHOTOSPLINE_CONVOLVE_H
#define PHOTOSPLINE_CONVOLVE_H

#include "photospline/splinetable.h"

namespace photospline{

double divdiff(const double *x, const double *y, size_t n);
unsigned int factorial(unsigned int n);
double convoluted_blossom(const double *x, size_t nx, const double *y, size_t ny,
                          double z, const double *bags, size_t nbags);

template <typename Alloc>
void splinetable<Alloc>::convolve(const uint32_t dim, const double* conv_knots, size_t n_conv_knots)
{
	/* Construct the new knot field. */
	size_t n_rho = 0;
	const uint32_t convorder = order[dim] + n_conv_knots - 1;
	std::unique_ptr<double[]> rho_scratch(new double[nknots[dim]*n_conv_knots + 2*convorder]);
	double* rho = rho_scratch.get() + convorder;
	for (uint32_t i = 0; i < nknots[dim]; i++)
		for (uint32_t j = 0; j < n_conv_knots; j++)
			rho[n_rho++] = this->knots[dim][i] + conv_knots[j];
	
	/* Order the new knot field and remove any duplicates. */
	std::sort(rho,rho+n_rho);
	//TODO: would removing duplicates, should they exist, be useful and safe?
	
	/* Set up space for the convolved coefficients */
	std::unique_ptr<uint64_t[]> naxes(new uint64_t[ndim]);
	std::unique_ptr<uint64_t[]> strides(new uint64_t[ndim]);
	
	std::copy(this->naxes,this->naxes+ndim,naxes.get());
	naxes[dim] = n_rho - convorder - 1;
	
	size_t arraysize = 1;
	strides[ndim - 1] = 1;
	for (int32_t i = ndim-1; i >= 0; i--) {
		arraysize *= naxes[i];
		if (i > 0)
			strides[i-1] = arraysize;
	}
	
	/*
	 * Now, calculate a transformation from coefficients on the raw knot
	 * field to coefficients on the convoluted knot field. Since the knots
	 * are on a
	 * grid, this transformation can be applied to each slice of the array.
	 * The coefficient of a spline in the convolved basis is a linear
	 * combination of the blossoms of each spline in the un-convolved basis
	 * convolved with the kernel spline. Here we just store the raw blossoms
	 * and multiply by the coefficients of each un-convolved spline later.
	 *
	 * This is analogous Stroem, Proposition 10, but with the prefactor
	 * adapted to account for the fact that one of the splines is de-Boor
	 * normalized (all supported splines add up to one -- "partition of
	 * unity") and the other unit normalized (each basis function integrates
	 * to one).
	 *
	 * NB: we're convolving a de-Boor spline with a unit-norm spline,
	 * hence q!(k-1)! rather than (q-1)!(k-1)! (as for two de-Boor splines).
	 */
	const uint32_t k = order[dim] + 1;
	const uint32_t q = n_conv_knots - 1;
	
	/* Norm is Stroem, Equation 13 */
	double norm = ((double)(factorial(q)*factorial(k-1)))/((double)factorial(k+q-1));
	
	std::unique_ptr<float[]> coefficients(new float[arraysize]);
	std::fill_n(coefficients.get(),arraysize,0.f);
	
	uint64_t stride1 = 1, stride2 = 1;
	for (uint32_t i = 0; i < ndim; i++) {
		if (i < dim)
			stride1 *= naxes[i];
		else if (i > dim)
			stride2 *= naxes[i];
	}
	
	std::unique_ptr<double[]> trafo(new double[naxes[dim]*this->naxes[dim]]);
	/*
	 * Fill the transformation matrix ahead of time to avoid recomputing
	 * it *stride1* times.
	 */
	for (uint32_t i = 0; i < naxes[dim]; i++) {
		for (uint32_t j = 0; j < this->naxes[dim]; j++) {
			trafo[i*this->naxes[dim] + j] =
			norm*convoluted_blossom(&this->knots[dim][j],
			                        k+1, conv_knots, n_conv_knots, rho[i], &rho[i+1], k+q-1);
		}
	}
	
	/*
	 * Multiply each vector of coefficients along dimension *dim*
	 * by the transformation matrix.
	 */
	for (uint32_t i = 0; i < stride1; i++)
		for (uint32_t j = 0; j < naxes[dim]; j++)
			for (uint32_t l = 0; l < this->naxes[dim]; l++)
				for (uint32_t k = 0; k < stride2; k++)
					coefficients[i*stride2*naxes[dim] + j*stride2 + k] +=
					trafo[j*this->naxes[dim] + l] *
					this->coefficients[i*stride2*this->naxes[dim] +
					l*stride2 + k];
	
	/*
	 * If the extent already had partial support at the lower end,
	 * let the new table extend to the limit of support. Otherwise,
	 * retain only full support.
	 */
	if (extents[dim][0] < this->knots[dim][order[dim]])
		extents[dim][0] = rho[0];
	else
		extents[dim][0] = rho[convorder];
	
	//In case we are using an allocator with limited total memory available
	//we need to avoid fragmentation. To do this, we need to deallocate all
	//memory currently used by the coefficiencts and knots before allocating
	//space for the new ones. Most of the old knot data we still need, so we
	//have to make temporary buffers for it.
	
	deallocate(this->coefficients,this->naxes[0]*this->strides[0]);
	
	std::unique_ptr<std::unique_ptr<double[]>[]> knots_store(new std::unique_ptr<double[]>[ndim]);
	for (uint32_t i = 0; i < ndim; i++) {
		//copy the old knots, except in the convolution dimension
		if (i!=dim) {
			knots_store[i].reset(new double[nknots[i]]);
			std::copy(knots[i],knots[i]+nknots[i],knots_store[i].get());
		}
		deallocate(knots[i]-order[i],nknots[i]+2*order[i]);
	}
	
	this->nknots[dim] = n_rho;
	this->order[dim] = convorder;
	this->naxes[dim] = naxes[dim];
	std::copy(strides.get(),strides.get()+ndim,this->strides);
	
	this->coefficients = allocate<float>(arraysize);
	std::copy(coefficients.get(),coefficients.get()+arraysize,this->coefficients);
	
	for (uint32_t i = 0; i < ndim; i++) {
		knots[i] = allocate<double>(nknots[i]+2*order[i]) + order[i];
		double* src = (i!=dim ? knots_store[i].get() : rho);
		std::copy(src,src+nknots[i],&knots[i][0]);
	}
	
	/*
	 * NB: A monotonic function remains monotonic after convolution
	 * with a strictly positive kernel. However, a spline cannot increase
	 * monotonically beyond its last fully-supported knot. Here, we reduce
	 * the extent of the spline by half the support of the spline kernel so
	 * that the surface will remain monotonic over its full extent.
	 */
	this->extents[dim][1] += conv_knots[0];
}

} //namespace photospline

#endif //PHOTOSPLINE_CONVOLVE_H
