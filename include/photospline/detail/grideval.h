#ifndef PHOTOSPLINE_DETAIL_GRIDEVAL_H
#define PHOTOSPLINE_DETAIL_GRIDEVAL_H

#include "photospline/splinetable.h"

namespace photospline{

template<typename Alloc>
template<typename DoubleContCont>
std::unique_ptr<ndsparse> splinetable<Alloc>::grideval(const DoubleContCont& coords) const{
	
	typedef typename DoubleContCont::value_type DoubleCont;
	static_assert(std::is_same<DoubleCont,typename std::remove_const<typename DoubleContCont::value_type>::type>::value,
	              "DoubleContCont must be a container of DoubleCont values");
	static_assert(std::is_same<double,typename std::remove_const<typename DoubleCont::value_type>::type>::value,
	              "DoubleCont must be a container of double values");

	if (coords.size() != ndim)
		throw(std::logic_error("Number of coordinate vectors ("
			+std::to_string(coords.size())+
			") must match dimensions ("+std::to_string(ndim)+")"));
	
	size_t size = naxes[0]*strides[0];
	size_t nnz = 0;
	for (size_t i=0; i<size; i++)
		if (coefficients[i] != 0)
			nnz++;
	std::unique_ptr<ndsparse> nd(new ndsparse(nnz, ndim));
	{
		std::vector<unsigned int> indices(ndim);
		for (size_t i=0; i<size; i++) {
			if (coefficients[i] != 0) {
				size_t coord = i;
				for (unsigned int dim = 0; dim < ndim; dim++) {
					indices[dim] = coord / strides[dim];
					coord = coord % strides[dim];
				}
				nd->insertEntry(coefficients[i], indices.data());
			}
		}
		// Set dimensions by hand to account for possible zeroes at the edges
		for (unsigned int dim = 0; dim < ndim; dim++) {
			nd->ranges[dim] = naxes[dim];
		}
		
	}
	
	cholmod_common cholmod_state;
	cholmod_l_start(&cholmod_state);
	
	for (unsigned i=0; i<ndim; i++) {
		const DoubleCont &coord_vec = coords[i];
		cholmod_sparse *basis, *basist;
		
		basis = bsplinebasis(knots[i], nknots[i],
		    coord_vec.data(), coord_vec.size(), order[i], &cholmod_state);
		basist = cholmod_l_transpose(basis, 1, &cholmod_state);
		cholmod_l_free_sparse(&basis, &cholmod_state);

		slicemultiply(nd.get(), basist, i, &cholmod_state);

		cholmod_l_free_sparse(&basist, &cholmod_state);
	}
	
	cholmod_l_finish(&cholmod_state);
	
	return std::move(nd);
}
	
}

#endif // PHOTOSPLINE_DETAIL_GRIDEVAL_H
