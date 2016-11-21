//Create a variety of splines for testing purposes

#include <iostream>
#include "photospline/splinetable-mod.h"

int main(){
	
	for(size_t dim=1; dim<6; dim++){
		std::cout << "Dimension " << dim  << std::endl;
		
		std::cout << "Generating data to fit..." << std::endl;
		
		std::vector<uint32_t> orders(dim,2);
		
		//pick a different domain to cover in each dimension
		std::vector<double> width(dim);
		//std::cout << " widths:" << std::endl;
		for(size_t i=0; i<dim; i++){
			width[i]=2*log10(i+1)+1;
			//std::cout << "  " << width[i] << std::endl;
		}
		//choose some knot locations
		std::vector<std::vector<double>> knots(dim);
		//std::cout << " knots:" << std::endl;
		for(size_t i=0; i<dim; i++){
			//std::cout << "  dimension " << i  << std::endl;
			size_t nknots=16;//(16-i*2);
			double offset=!(nknots%2);
			double scale=2*width[i]/nknots;
			for(size_t j=0; j<nknots; j++){
				double knot=scale*(j+.5-nknots/2.);
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
		size_t total_samples=1;
		for(size_t i=0; i<dim; i++)
			total_samples*=nsamples;
		photospline::ndsparse data(total_samples,dim);
		std::vector<double> weights(total_samples,1.);
		
		std::vector<unsigned int> indices(dim,0);
		for(size_t i=0; i<total_samples; i++){
			double value=1;
			for(size_t j=0; j<dim; j++){
				double x=coordinates[j][indices[j]];
				double s=(j+1)/2.;
				value*=exp(-(x*x)/(s*s)/2);
			}
			data.insertEntry(value,&indices.front());
			
			for(size_t j=0; j<dim; j++){
				if((++indices[dim-j-1])>=nsamples)
					indices[dim-j-1]=0;
				else
					break;
			}
		}
		
		size_t dataSize=0;
		dataSize+=total_samples*sizeof(double);
		dataSize+=total_samples*sizeof(double);
		dataSize+=total_samples*dim*sizeof(unsigned int);
		dataSize+=40*dim*sizeof(double);
		std::cout << "Input data is " << dataSize << " bytes" << std::endl;
		size_t ncoeffs=1;
		for(size_t i=0; i<dim; i++)
			ncoeffs*=knots[i].size()-orders[i]-1;
		
		std::cout << "Output will have " << ncoeffs << " coefficients" << std::endl;
		
		photospline::splinetable<> spline;
		spline.fit(data,weights,coordinates,orders,knots,{dim*1e-8},{2});
		spline.write_fits("test_spline_"+std::to_string(dim)+"d.fits");
		std::cout << "Output has " << spline.get_ncoeffs() << " coefficients" << std::endl;
		spline.benchmark_evaluation(1e4,true);
		
		//-----
		orders.back()=3;
		
		photospline::splinetable<> spline2;
		spline2.fit(data,weights,coordinates,orders,knots,{dim*1e-10},{2});
		spline2.write_fits("test_spline_"+std::to_string(dim)+"d_nco.fits");
		std::cout << "Output has " << spline2.get_ncoeffs() << " coefficients" << std::endl;
		spline2.benchmark_evaluation(1e4,true);
	}
}