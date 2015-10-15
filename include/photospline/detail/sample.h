#ifndef PHOTOSPLINE_SAMPLE_H
#define PHOTOSPLINE_SAMPLE_H

#include "photospline/splinetable-mod.h"

#include <array>
#include <cassert>
#include <random>

namespace photospline{

///\tparam N the number of dimensions in which to sample
///\tparam Distribution a type which can both evaluate a proposal pdf
///        via double operator()(const std::vector<double>& coordinates) const
///        and can draw samples from that same pdf
///        via std::vector<double> sample(RNG rng)
///\tparam RNG a random number generator, conforming to the standard interface
///\tparam Transform 
///\todo document parameters
///\pre The entries of coordinates for dimensions which do not appear in 
///        samplingDimensions must be with in the extents of spline
template<typename Alloc>
template<size_t N, typename Distribution, typename RNG, typename Transform>
std::vector<std::array<double,N>> splinetable<Alloc>::sample(
  size_t nresults, size_t burnin, std::array<size_t,N> samplingDimensions,
  std::vector<double> coordinates, Distribution distribution, RNG& rng,
  Transform transform) const{
	bool accept;
	const uint32_t tN=get_ndim();
	std::vector<double> x(tN), xp(tN);
	std::array<double,N> s, sp;
	std::vector<int> cx(tN), cxp(tN);
	double px, pxp, propx, propxp;

	//check preconditions
	assert(coordinates.size()==tN);
	for(uint32_t i=0; i<tN; i++){
		//if this dimension is not being sampled, the user supplied coordinate needs to be valid
		if(std::find(samplingDimensions.begin(),samplingDimensions.end(),i)==samplingDimensions.end()){
			assert(coordinates[i]>=lower_extent(i));
			assert(coordinates[i]<=upper_extent(i));
		}
	}

	x = xp = coordinates;
	//sample an intial point
	do{
		s=distribution.sample(rng);
		//assert(s.size()==N);
		accept=true;
		//demand that the sampled point be within the table extents
		for(uint32_t i=0; i<N; i++){
			if(s[i]<lower_extent(samplingDimensions[i])
			   || s[i]>upper_extent(samplingDimensions[i])){
				accept=false;
				break;
			}
		}

		if(accept){
			//copy the sampled coordinates into their places in the full set
			for(uint32_t i=0; i<N; i++)
				x[samplingDimensions[i]]=s[i];
		
			accept=searchcenters(x.data(),cx.data());
		}
	} while(!accept);

	propx=distribution(s);
	px=transform(ndsplineeval(x.data(),cx.data(),0));

	std::uniform_real_distribution<> acceptor(0.,1.);

	std::vector<std::array<double,N>> results(nresults);
	//for the number of results requested by the user
	for(size_t k=0; k<nresults; k++){
		for(size_t j=0; j<=burnin; j++){
			sp=distribution.sample(rng);
			accept=true;
			//demand that the sampled point be within the table extents
			for(uint32_t i=0; i<N; i++){
				if(sp[i]<lower_extent(samplingDimensions[i])
				   || sp[i]>upper_extent(samplingDimensions[i])){
					accept=false;
					break;
				}
			}
			if(!accept)
				continue;
	
			//copy the sampled coordinates into their places in the full set
			for(uint32_t i=0; i<N; i++)
				xp[samplingDimensions[i]]=sp[i];
			
			accept=searchcenters(xp.data(),cxp.data());
			if(!accept)
				continue;
	
			propxp=distribution(sp);
			pxp=transform(ndsplineeval(xp.data(),cxp.data(),0));
	
			double odds=(pxp/px)*(propx/propxp);
			accept=((odds>1.) || acceptor(rng)<odds);

			if(accept){
				s=sp;
				x=xp;
				propx=propxp;
				px=pxp;
			}
		}
		//record whatever we sampled
		for(size_t i=0; i<N; i++)
			results[k][i]=s[i];
	}

	return(results);
}

} //namesapce photospline

#endif            
