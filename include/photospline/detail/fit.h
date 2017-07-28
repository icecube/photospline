#ifndef PHOTOSPLINE_DETAIL_FIT_H
#define PHOTOSPLINE_DETAIL_FIT_H

#include "photospline/splinetable.h"

namespace photospline{

template<typename Alloc>
template<typename DoubleCont, typename UInt32Cont, typename DoubleContCont>
void splinetable<Alloc>::fit(const ::ndsparse& data,
                             const DoubleCont& weights,
                             const DoubleContCont& coords,
                             const UInt32Cont& splineOrder,
                             const DoubleContCont& knots,
                             const DoubleCont& smoothing,
                             const UInt32Cont& penaltyOrder,
                             uint32_t monodim, bool verbose){
	static_assert(std::is_same<double,typename std::remove_const<typename DoubleCont::value_type>::type>::value,
	              "DoubleCont must be a container of double values");
	static_assert(std::is_same<uint32_t,typename std::remove_const<typename UInt32Cont::value_type>::type>::value,
	              "UInt32Cont must be a container of uint32_t values");
	static_assert(std::is_same<DoubleCont,typename std::remove_const<typename DoubleContCont::value_type>::type>::value,
	              "DoubleContCont must be a container of DoubleCont values");
	
	//Sanity checking
	if(data.rows!=weights.size())
		throw std::logic_error("Number of weights ("
		                       +std::to_string(weights.size())
		                       +") does not equal number of data points ("
		                       +std::to_string(data.rows)+")");
	for(uint32_t i=0; i<data.ndim; i++){
		unsigned int maxIdx=*std::max_element(data.i[i],data.i[i]+data.rows);
		if(maxIdx>=data.ranges[i])
			throw std::logic_error("Range of coordinate indices ("
			                       +std::to_string(data.ranges[i])
			                       +") in dimension "+std::to_string(i)
			                       +" does not include the maximum index used ("
			                       +std::to_string(maxIdx)+")");
	}
	if(coords.size()!=data.ndim)
		throw std::logic_error("Number of coordinate vectors ("
		                       +std::to_string(coords.size())
		                       +") does not equal dimension of input data ("
		                       +std::to_string(data.ndim)+")");
	if(splineOrder.size()!=data.ndim)
		throw std::logic_error("Number of spline orders ("
		                       +std::to_string(splineOrder.size())
		                       +") does not equal dimension of input data ("
		                       +std::to_string(data.ndim)+")");
	if(knots.size()!=data.ndim)
		throw std::logic_error("Number of knot vectors ("
		                       +std::to_string(knots.size())
		                       +") does not equal dimension of input data ("
		                       +std::to_string(data.ndim)+")");
	for(uint32_t i=0; i<data.ndim; i++){
		if(!std::is_sorted(knots[i].begin(),knots[i].end()))
			throw std::logic_error("Knot vector for dimension "
			                       +std::to_string(i)+
			                       " is not in sorted order");
	}
	if(smoothing.size()!=data.ndim && smoothing.size()!=1)
		throw std::logic_error("Number of smoothing strengths specified ("
		                       +std::to_string(smoothing.size())
		                       +") should be 1 or the number of spline dimensions ("
		                       +std::to_string(data.ndim)+")");
	if(penaltyOrder.size()!=data.ndim && penaltyOrder.size()!=1)
		throw std::logic_error("Number of penalty orders specified ("
		                       +std::to_string(penaltyOrder.size())
		                       +") should be 1 or the number of spline dimensions ("
		                       +std::to_string(data.ndim)+")");
	if(monodim!=no_monodim && monodim>=data.ndim)
		throw std::logic_error("Requested monotonic dimension ("
		                       +std::to_string(monodim)
		                       +") shoulb be less than the number of spline dimensions ("
		                       +std::to_string(data.ndim)+")");
	
	//Initialize variables
	ndim=data.ndim;
	order = allocate<uint32_t>(ndim);
	std::copy(splineOrder.begin(),splineOrder.end(),order);
	this->knots = allocate<double_ptr>(ndim);
	nknots = allocate<uint64_t>(ndim);
	for(uint32_t i=0; i<ndim; i++)
		nknots[i]=knots[i].size();
	extents = allocate<double_ptr>(ndim);
	extents[0] = allocate<double>(2*ndim);
	naxes = allocate<uint64_t>(ndim);
	for(uint32_t i=0; i<ndim; i++)
		naxes[i]=nknots[i]-order[i]-1;
	strides = allocate<uint64_t>(ndim);
	strides[ndim-1]=1;
	for(uint32_t i = ndim-1; i>0; i--)
		strides[i-1] = strides[i]*naxes[i];
	uint64_t ncoeffs=strides[0]*naxes[0];
	coefficients = allocate<float>(ncoeffs);
	
	//glamfit_complex really wants a double** for the knots
	//so we set up this shim
	std::unique_ptr<double*[]> dummy_knots(new double*[ndim]);
	for(uint32_t i=0; i<ndim; i++){
		this->knots[i]=allocate<double>(nknots[i]+2*order[i]) + order[i];
		std::copy(knots[i].begin(),knots[i].end(),this->knots[i]);
		dummy_knots[i]=&this->knots[i][0];
	}
	//same deal for the coordinates
	std::unique_ptr<const double*[]> dummy_coords(new const double*[ndim]);
	for(uint32_t i=0; i<ndim; i++)
		dummy_coords[i]=coords[i].data();
	
	//make up extents, hopefully reasonably
	for(uint32_t i=0; i<ndim; i++){
		extents[i] = &extents[0][2*i];
		extents[i][0] = this->knots[i][order[i]];
		extents[i][1] = this->knots[i][nknots[i] - order[i] - 1];
	}
	
	cholmod_common cholmod_state;
	cholmod_l_start(&cholmod_state);
	//build penalty matrix
	cholmod_sparse* penalty;
	{
		std::unique_ptr<uint64_t[]> nsplines(new uint64_t[ndim]);
		uint64_t sidelen = 1;
		for(uint32_t i = 0; i < ndim; i++) {
			nsplines[i] = nknots[i] - order[i] - 1;
			sidelen *= nsplines[i];
		}
		penalty=cholmod_l_spzeros(sidelen, sidelen, 1, CHOLMOD_REAL, &cholmod_state);
		for(uint32_t i = 0; i < ndim; i++){
			penalty = add_penalty_term(nsplines.get(), &this->knots[i][0], ndim,
			                           i, order[i],
			                           (penaltyOrder.size()>1?penaltyOrder[i]:penaltyOrder[0]),
			                           (smoothing.size()>1?smoothing[i]:smoothing[0]),
			                           i==monodim, penalty,
			                           &cholmod_state);
		}
	}
	
	//do the fit
	int result=
	glamfit_complex(&data,weights.data(),dummy_coords.get(),
	                ndim,&nknots[0],dummy_knots.get(),&naxes[0],
	                &coefficients[0],
	                &order[0],penalty,(monodim==no_monodim?-1:(int)monodim),
	                verbose,&cholmod_state);
	//clean up
	cholmod_l_free_sparse(&penalty, &cholmod_state);
	cholmod_l_finish(&cholmod_state);
	if(result!=0)
		throw std::runtime_error("GLAM fit failed");
}
	
} //namespace photospline

#endif