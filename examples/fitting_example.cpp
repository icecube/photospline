
#include <photospline/splinetable.h>
#include <random>

int main (int argc, char const *argv[])
{
	std::default_random_engine rng(0);
	
	size_t numpts = 500;
	// Create a random coordinate vector
	std::vector<double> x1(numpts);
	{
		std::uniform_real_distribution<> uniform(-4,25);
		for (unsigned int i=0; i < numpts; i++)
			x1[i] = uniform(rng);
		std::sort(x1.begin(), x1.end());
	}
	// Create a sparse array to hold the data points and a vector to hold the
	// weights. The length of `weights` should be the number of entries in
	// `zsparse` with non-zero weights
	photospline::ndsparse zsparse(500,1);
	std::vector<double> weights(numpts);
	for (unsigned int i=0; i < numpts; i++) {
		double z = std::poisson_distribution<>(
		    std::pow(std::cos(x1[i]), 2) + std::pow(x1[i]-3.0, 2) + 10)(rng);
		zsparse.insertEntry(z, &i);
		// Use a minimum weight of zero, and weight the high-occupancy points up
		weights[i] = 1. + z;
	}
	
	std::vector<double> knots(30);
	for (unsigned int i=0; i < knots.size(); i++) {
		double lo(-8), hi(35);
		knots[i] = lo + i*(hi-lo)/(knots.size()-1);
	}
	
	std::array<uint32_t,1> order = {2};
	auto penalty_order = order;
	double smooth = 3.14159e3;
	
	photospline::splinetable<> spline;
	spline.fit(zsparse, weights, std::array<decltype(x1),1>({x1}), order, {knots}, {smooth}, penalty_order);
	
	spline.write_fits("splinefit-1d.fits");
	
	/* code */
	return 0;
}