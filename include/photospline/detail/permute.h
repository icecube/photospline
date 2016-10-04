#ifndef PHOTOSPLINE_PERMUTE_H
#define PHOTOSPLINE_PERMUTE_H

#include "photospline/splinetable.h"

namespace photospline{
	
template<typename Alloc>
void splinetable<Alloc>::permuteDimensions(const std::vector<size_t>& permutation){
	{
		if(permutation.size()!=ndim)
			throw std::runtime_error("Wrong number of indices passed to permuteDimensions");
		std::vector<bool> permutation_test(permutation.size(),false);
		for(size_t i=0; i<permutation.size(); i++){
			size_t j=permutation[i];
			if(j>=ndim)
				throw std::runtime_error("Too large index passed to permuteDimensions");
			if(permutation_test[j])
				throw std::runtime_error("Duplicate index passed to permuteDimensions");
			permutation_test[j]=true;
		}
		for(size_t i=0; i<permutation_test.size(); i++){
			if(!permutation_test[i])
				throw std::runtime_error("Missing index in permutation passed to permuteDimensions");
		}
	}
	
	//Note that we use regular pointers because these allocations will be 'local'
	//to this function.
	std::unique_ptr<uint32_t[]> t_order(new uint32_t[ndim]);
	std::vector<uint64_t> t_naxes(ndim);
	std::unique_ptr<uint64_t[]> t_strides(new uint64_t[ndim]);
	std::unique_ptr<uint64_t[]> t_nknots(new uint64_t[ndim]);
	std::unique_ptr<double_ptr[]> t_knots(new double_ptr[ndim]);
	std::unique_ptr<double*[],void(*)(double**)> t_extents(new double*[ndim],
		[](double** p){
			if(p && p[0])
				delete[] p[0];
			delete[] p;
		});
	t_extents[0]=new double[2*ndim];
	for(uint32_t i=1; i<ndim; i++)
		t_extents[i] = &t_extents[0][2*i];
	
	// Permute various per-axis properties
	uint32_t iperm[ndim];
	for(uint32_t i=0; i<ndim; i++){
		uint32_t j = permutation[i];
		iperm[j] = i;
		t_order[i] = order[j];
		t_naxes[i] = naxes[j];
		t_nknots[i] = nknots[j];
		t_knots[i] = knots[j];
		t_extents[i][0] = extents[j][0];
		t_extents[i][1] = extents[j][1];
	}
	
	// Compute new strides
	t_strides[0]=1;
	std::partial_sum(t_naxes.rbegin(),t_naxes.rend()-1,t_strides.get()+1,std::multiplies<uint64_t>());
	std::reverse(t_strides.get(),t_strides.get()+ndim);
	uint64_t ncoeffs=t_strides[0]*t_naxes[0];
	
	// Re-order coefficient array
	std::unique_ptr<float[]> t_coefficients(new float[ncoeffs]);
	for(uint64_t pos=0; pos<ncoeffs; pos++){
		uint64_t npos = 0;
		// Multiply index of point in old shape by new stride
		for(uint32_t i=0; i<ndim; i++)
			npos += (pos / strides[i] % naxes[i])*t_strides[iperm[i]];
		t_coefficients[npos] = coefficients[pos];
	}
	
	// Copy all data back to main storage
	std::copy(t_order.get(),t_order.get()+ndim,order);
	std::copy(t_naxes.begin(),t_naxes.end(),naxes);
	std::copy(t_strides.get(),t_strides.get()+ndim,strides);
	std::copy(t_nknots.get(),t_nknots.get()+ndim,nknots);
	std::copy(t_knots.get(),t_knots.get()+ndim,knots);
	for(uint32_t i=0; i<ndim; i++){
		extents[i][0]=t_extents[i][0];
		extents[i][1]=t_extents[i][1];
	}
	std::copy(t_coefficients.get(),t_coefficients.get()+ncoeffs,coefficients);
}
	
} //namespace photospline

#endif