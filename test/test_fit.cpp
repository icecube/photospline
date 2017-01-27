#include "test.h"

#include "photospline/splinetable-mod.h"

TEST(glamfit_2d){
	const uint32_t dim=2;
	std::vector<uint32_t> orders(dim,2);
	std::vector<std::vector<double> > knots(dim);
	for(size_t i=0; i<dim; i++){
		size_t nknots=69;
		for(size_t j=0; j<nknots; j++){
			double knot=-5.6666666666+(11.3333333333/(nknots-1))*j;
			knots[i].push_back(knot);
		}
	}
	
	const size_t nsamples=100;
	std::vector<std::vector<double> > coordinates(dim);
	for(size_t i=0; i<dim; i++){
		double step=12./(nsamples-1);
		for(size_t j=0; j<nsamples; j++){
			coordinates[i].push_back(-6+j*step);
		}
	}
	
	//make some data to fit
	size_t total_samples=1;
	for(size_t i=0; i<dim; i++)
		total_samples*=coordinates[i].size();
	photospline::ndsparse data(total_samples,dim);
	std::vector<double> weights(total_samples,1.);
	
	std::vector<unsigned int> indices(dim,0);
	for(size_t i=0; i<total_samples; i++){
		for(size_t j=0; j<dim; j++){
			if((++indices[dim-j-1])>=coordinates[j].size())
				indices[dim-j-1]=0;
			else
				break;
		}
		
		double value=cos(coordinates[0][indices[0]])
		             *sin(coordinates[1][indices[1]]);
		data.insertEntry(value,&indices.front());
	}
	
	std::cout << "Fitting..." << std::endl;
	photospline::splinetable<> spline;
	spline.fit(data,weights,coordinates,orders,knots,{0},{2});
	
	std::cout << "Testing fit..." << std::endl;
	std::vector<double> coords(dim);
	std::vector<int> centers(dim);
	
	double maxErr=0;
	
	indices=decltype(indices)(dim,0);
	for(size_t i=0; i<total_samples; i++){
		for(size_t j=0; j<dim; j++){
			if((++indices[dim-j-1])>=coordinates[j].size())
				indices[dim-j-1]=0;
			else
				break;
		}
		
		coords[0]=coordinates[0][indices[0]];
		coords[1]=coordinates[1][indices[1]];
		
		//The edges of this fit don't work very well, probably because the
		//target function is oscillitory. Only test in the interior.
		if(coords[0]<-4. || coords[0]>4.
		   || coords[1]<-4. || coords[1]>4.){
			continue;
		}
		
		double reference=cos(coordinates[0][indices[0]])
		*sin(coordinates[1][indices[1]]);
		
		ENSURE(spline.searchcenters(coords.data(), centers.data()), "Center lookup should succeed");
		double evaluate=spline.ndsplineeval(coords.data(), centers.data(), 0);
		
		if(std::abs(evaluate-reference)>std::abs(maxErr))
			maxErr=(evaluate-reference);
		ENSURE_DISTANCE(evaluate,reference,1e-4,"spline value should be close to true function value");
	}
	std::cout << "Maximimum error: " << maxErr << std::endl;
}

//Ensure that we can do a fit on inputs which we don't copy (important for
//making the bindings to other languges efficient)
TEST(glam_fit_unowned_data){
	const uint32_t dim=1;
	const double knot_data[10]={1.,2.,3.,4.,5.,6.,7.,8.,9.,10.};
	const uint32_t order_data[1]={2};
	const uint32_t porder_data[1]={2};
	double penalty_data[1]={0};
	double coord_data[100];
	double weight_data[100];
	photospline::ndsparse data(100,dim);
	for(unsigned int i=0; i<100; i++){
		double x=i/10.;
		coord_data[i]=x;
		data.insertEntry(x,&i);
		weight_data[i]=1;
	}
	
	using DoubleCont=photospline::detail::array_view<double>;
	using UInt32Cont=photospline::detail::array_view<uint32_t>;
	using DoubleContCont=std::vector<DoubleCont>;
	
	DoubleContCont knots={DoubleCont(knot_data,10)};
	UInt32Cont orders(order_data,1);
	DoubleContCont coordinates={DoubleCont(coord_data,100)};
	DoubleCont weights(weight_data,100);
	UInt32Cont porders(porder_data,1);
	DoubleCont penalties(penalty_data,1);
	
	std::cout << "Fitting..." << std::endl;
	photospline::splinetable<> spline;
	spline.fit(data,weights,coordinates,orders,knots,penalties,porders);
	
	std::cout << "Testing fit..." << std::endl;
	std::vector<double> coords(dim);
	std::vector<int> centers(dim);
}
