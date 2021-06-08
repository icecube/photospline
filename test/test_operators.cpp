#include "test.h"

#include "photospline/splinetable.h"

namespace{
	struct makeOpts{
		unsigned int nKnots;
		double knotScale;
		double dataTwiddle;
		
		makeOpts():nKnots(69),knotScale(11.33333333333),dataTwiddle(1){}
	};
	
	photospline::splinetable<> makeSpline(const makeOpts& opts){
		ENSURE(opts.nKnots%2==1,"Number of knots must be odd");
		const uint32_t dim=2;
		std::vector<uint32_t> orders(dim,2);
		std::vector<std::vector<double> > knots(dim);
		const double knotOffset=-(opts.knotScale/(opts.nKnots-1))*(opts.nKnots/2);
		for(size_t i=0; i<dim; i++){
			for(size_t j=0; j<opts.nKnots; j++){
				double knot=knotOffset+(opts.knotScale/(opts.nKnots-1))*j;
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
			
			double value=opts.dataTwiddle*cos(coordinates[0][indices[0]])
			*sin(coordinates[1][indices[1]]);
			data.insertEntry(value,&indices.front());
		}
		
		std::cout << "Fitting..." << std::endl;
		photospline::splinetable<> spline;
		spline.fit(data,weights,coordinates,orders,knots,{0},{2});
		return spline;
	}
}

TEST(Equality){
	{
		photospline::splinetable<> spline1("test_data/test_spline_4d.fits");
		photospline::splinetable<> spline2("test_data/test_spline_4d.fits");
		ENSURE(spline1==spline2, "Splines loaded from same data should compare equal");
	}
	{
		photospline::splinetable<> spline1("test_data/test_spline_1d.fits");
		photospline::splinetable<> spline2("test_data/test_spline_2d.fits");
		ENSURE(!(spline1==spline2), 
		       "Splines with different numbers of dimensions should not compare equal");
	}
	{
		photospline::splinetable<> spline1("test_data/test_spline_1d.fits");
		photospline::splinetable<> spline2("test_data/test_spline_2d.fits");
		ENSURE(!(spline1==spline2), 
		       "Splines with different numbers of dimensions should not compare equal");
	}
	{
		photospline::splinetable<> spline1("test_data/test_spline_2d.fits");
		photospline::splinetable<> spline2("test_data/test_spline_2d_nco.fits");
		ENSURE(!(spline1==spline2), 
		       "Splines with different orders should not compare equal");
	}
	{
		photospline::splinetable<> spline1=makeSpline(makeOpts());
		photospline::splinetable<> spline2=makeSpline(makeOpts());
		ENSURE(spline1==spline2, "Splines fit with same settings should compare equal");
	}
	{
		makeOpts opts;
		photospline::splinetable<> spline1=makeSpline(opts);
		opts.nKnots=61;
		photospline::splinetable<> spline2=makeSpline(opts);
		ENSURE(!(spline1==spline2), "Splines with different numbers of knots should not compare equal");
	}
	{
		makeOpts opts;
		photospline::splinetable<> spline1=makeSpline(opts);
		opts.knotScale=10;
		photospline::splinetable<> spline2=makeSpline(opts);
		ENSURE(!(spline1==spline2), "Splines with different knots should not compare equal");
	}
	{
		makeOpts opts;
		photospline::splinetable<> spline1=makeSpline(opts);
		opts.dataTwiddle=2;
		photospline::splinetable<> spline2=makeSpline(opts);
		ENSURE(!(spline1==spline2), "Splines with different coefficients should not compare equal");
	}
}

TEST(Inequality){
	{
		photospline::splinetable<> spline1("test_data/test_spline_4d.fits");
		photospline::splinetable<> spline2("test_data/test_spline_4d.fits");
		ENSURE(!(spline1!=spline2), "Splines loaded from same data should not compare unequal");
	}
	{
		photospline::splinetable<> spline1("test_data/test_spline_1d.fits");
		photospline::splinetable<> spline2("test_data/test_spline_2d.fits");
		ENSURE(spline1!=spline2, 
		       "Splines with different numbers of dimensions should compare unequal");
	}
	{
		photospline::splinetable<> spline1("test_data/test_spline_1d.fits");
		photospline::splinetable<> spline2("test_data/test_spline_2d.fits");
		ENSURE(spline1!=spline2, 
		       "Splines with different numbers of dimensions should compare unequal");
	}
	{
		photospline::splinetable<> spline1("test_data/test_spline_2d.fits");
		photospline::splinetable<> spline2("test_data/test_spline_2d_nco.fits");
		ENSURE(spline1!=spline2, 
		       "Splines with different orders should compare unequal");
	}
	{
		photospline::splinetable<> spline1=makeSpline(makeOpts());
		photospline::splinetable<> spline2=makeSpline(makeOpts());
		ENSURE(!(spline1!=spline2), "Splines fit with same settings should not compare unequal");
	}
	{
		makeOpts opts;
		photospline::splinetable<> spline1=makeSpline(opts);
		opts.nKnots=61;
		photospline::splinetable<> spline2=makeSpline(opts);
		ENSURE(spline1!=spline2, "Splines with different numbers of knots should compare unequal");
	}
	{
		makeOpts opts;
		photospline::splinetable<> spline1=makeSpline(opts);
		opts.knotScale=10;
		photospline::splinetable<> spline2=makeSpline(opts);
		ENSURE(spline1!=spline2, "Splines with different knots should compare unequal");
	}
	{
		makeOpts opts;
		photospline::splinetable<> spline1=makeSpline(opts);
		opts.dataTwiddle=2;
		photospline::splinetable<> spline2=makeSpline(opts);
		ENSURE(spline1!=spline2, "Splines with different coefficients should compare unequal");
	}
}
