//Create a variety of splines for testing purposes

#include <iostream>
#include "photospline/splinetable-mod.h"

int main(){
	
	for(size_t dim=1; dim<7; dim++){
		std::cout << "Dimension " << dim  << std::endl;
		photospline::ndsparse data;
		std::vector<double> weights;
		
		std::vector<uint32_t> orders(dim,2);
		
		//pick a different domain to cover in each dimension
		std::vector<double> width(dim);
		//std::cout << " widths:" << std::endl;
		for(size_t i=0; i<dim; i++){
			width[i]=2*log10(i+1)+1;
			//std::cout << "  " << width[i] << std::endl;
		}
		//choose some knot locations
		std::vector<std::vector<double>> knots;
		//std::cout << " knots:" << std::endl;
		for(size_t i=0; i<dim; i++){
			//std::cout << "  dimension " << i  << std::endl;
			size_t nknots=(20-i*2);
			double offset=!(nknots%2);
			double scale=(4./((nknots-1)*(nknots-1)))*width[i];
			for(size_t j=0; j<nknots; j++){
				double knot=(j<(nknots/2)?-1:1)*scale*(j+.5-nknots/2.)*(j+.5-nknots/2.);
				//std::cout << "   " << knot << std::endl;
				knots[i].push_back(knot);
			}
		}
		//define the grid on which we will sample
		const size_t nsamples=40;
		std::vector<std::vector<double>> coordinates(dim);
		for(size_t i=0; i<dim; i++){
			double step=2*width[i]/(nsamples-1);
			for(size_t j=0; j<nsamples; j++)
				coordinates[i].push_back(-width[i]+j*step);
		}
		
		//invent some data to fit
		//const size_t total_samples=pow(nsamples,dim);
		//for(size_t i=0; i<total_samples)
		
		//-----
		orders.back()=3;
		
	}
}