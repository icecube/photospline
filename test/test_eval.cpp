#include "test.h"

#include <random>
#include <vector>

#include "photospline/splinetable-mod.h"

TEST(ndssplineeval_vs_ndssplineeval_gradient){
	photospline::splinetable<> spline("test_data/test_spline_4d.fits");
	const int ndim = spline.get_ndim();
	ENSURE(ndim < 6);
	
	std::mt19937 rng;
	rng.seed(42);
	
	//Create uniform distributions over the support of the spline in all dimensions
	std::vector<std::uniform_real_distribution<>> dists;
	for(size_t i=0; i<ndim; i++)
		dists.push_back(std::uniform_real_distribution<>(spline.lower_extent(i),spline.upper_extent(i)));
	
	std::vector<double> coords(ndim);
	std::vector<int> centers(ndim);
	std::vector<double> gradient(ndim);
	std::vector<double> evaluate_with_gradient(ndim+1);
	
	//Evaluate the spline both ways at random points
	for(size_t i=0; i<10000; i++) {
		for(size_t j=0; j<ndim; j++)
			coords[j]=dists[j](rng);
		
		ENSURE(spline.searchcenters(coords.data(), centers.data()), "Center lookup should succeed");
		double evaluate=spline.ndsplineeval(coords.data(), centers.data(), 0);
		for(uint32_t i=0; i<ndim; i++)
			gradient[i]=spline.ndsplineeval(coords.data(), centers.data(), 1u<<i);
		
		spline.ndsplineeval_gradient(coords.data(), centers.data(), evaluate_with_gradient.data());
		
		ENSURE_EQUAL(evaluate, evaluate_with_gradient[0],
					 "ndsplineeval() and ndssplineeval_gradient() yield identical evaluates");
		ENSURE_DISTANCE(evaluate, evaluate_with_gradient[0],std::abs(1e-4*evaluate),
						"ndsplineeval() and ndssplineeval_gradient() yield identical evaluates");
		for (int j=0; j < ndim; j++) {
			ENSURE_EQUAL(gradient[j], evaluate_with_gradient[j+1],
						 "ndsplineeval() and ndssplineeval_gradient() yield identical derivatives");
			ENSURE_DISTANCE(gradient[j], evaluate_with_gradient[j+1],std::abs(1e-4*gradient[j]),
							"ndsplineeval() and ndssplineeval_gradient() yield identical derivatives");
		}
	}
}

void test_evaluator_interface(const std::string& splinePath){
	std::cout << "Testing evaluation of " << splinePath << std::endl;
	photospline::splinetable<> spline(splinePath);
	photospline::splinetable<>::evaluator evaluator=spline.get_evaluator();
	const int ndim = spline.get_ndim();
	ENSURE(ndim < 6);
	
	std::mt19937 rng;
	rng.seed(52);
	
	//Create uniform distributions over the support of the spline in all dimensions
	std::vector<std::uniform_real_distribution<>> dists;
	for(size_t i=0; i<ndim; i++)
		dists.push_back(std::uniform_real_distribution<>(spline.lower_extent(i),spline.upper_extent(i)));
	
	std::vector<double> coords(ndim);
	std::vector<int> centers1(ndim), centers2(ndim);
	std::vector<double> gradient1(ndim), gradient2(ndim);
	std::vector<double> evaluate_with_gradient1(ndim+1), evaluate_with_gradient2(ndim+1);
	
	//Evaluate the spline with and without the evaluator
	for(size_t i=0; i<10000; i++) {
		for(size_t j=0; j<ndim; j++)
			coords[j]=dists[j](rng);
		
		ENSURE(spline.searchcenters(coords.data(), centers1.data()), "Center lookup should succeed");
		double evaluate1=spline.ndsplineeval(coords.data(), centers1.data(), 0);
		for(uint32_t i=0; i<ndim; i++)
			gradient1[i]=spline.ndsplineeval(coords.data(), centers1.data(), 1u<<i);
		
		spline.ndsplineeval_gradient(coords.data(), centers1.data(), evaluate_with_gradient1.data());
		
		ENSURE(evaluator.searchcenters(coords.data(), centers2.data()), "Center lookup should succeed");
		double evaluate2=spline.ndsplineeval(coords.data(), centers2.data(), 0);
		for(uint32_t i=0; i<ndim; i++)
			gradient2[i]=evaluator.ndsplineeval(coords.data(), centers2.data(), 1u<<i);
		
		evaluator.ndsplineeval_gradient(coords.data(), centers2.data(), evaluate_with_gradient2.data());
		
		for (int j=0; j < ndim; j++) {
			ENSURE_EQUAL(centers1[j],centers2[j],"Center lookups should yield same results");
		}
		ENSURE_EQUAL(evaluate1, evaluate2,
					 "splinetable::ndsplineeval() and evaluator::ndsplineeval() yield identical evaluates");
		ENSURE_EQUAL(evaluate_with_gradient1[0], evaluate_with_gradient2[0],
					 "splinetable::ndssplineeval_gradient() and evaluator::ndssplineeval_gradient() yield identical evaluates");
		for (int j=0; j < ndim; j++) {
			ENSURE_EQUAL(gradient1[j], gradient2[j],
						 "splinetable::ndsplineeval() and evaluator::ndsplineeval() yield identical derivatives");
			ENSURE_EQUAL(evaluate_with_gradient1[j+1], evaluate_with_gradient2[j+1],
						 "splinetable::ndssplineeval_gradient() and evaluator::ndssplineeval_gradient() yield identical derivatives");
		}
	}
}

TEST(evaluator_interface){
	test_evaluator_interface("test_data/test_spline_2d.fits");
	test_evaluator_interface("test_data/test_spline_2d_nco.fits");
	test_evaluator_interface("test_data/test_spline_3d.fits");
	test_evaluator_interface("test_data/test_spline_3d_nco.fits");
	test_evaluator_interface("test_data/test_spline_4d.fits");
	test_evaluator_interface("test_data/test_spline_4d_nco.fits");
}

TEST(bsplvb_simple_vs_bspline){
	const size_t n_knots = 10;
	const int order = 2;
	double x;
	const double* knots;
	int center, offset;
	std::vector<float> localbasis_bsplvb(order+1);
	std::vector<float> localbasis_bspline(order+1);
	std::vector<double> knotvec;
	
	std::mt19937 rng;
	rng.seed(142);
	
	// Generate a random knot field.
	{
		std::uniform_real_distribution<> uniform(0,1);
		for (size_t i = 0; i < n_knots; i++)
			knotvec.push_back(uniform(rng));
	}
	std::sort(knotvec.begin(), knotvec.end());
	knots = knotvec.data();
	
	// Before the first fully-supported knot.
	for (size_t i = 0; i < order+1; i++) {
		x = (knots[i]+knots[i+1])/2.0; // Less than fully-supported
		center = order; // First fully-supported spline.
		
		ENSURE(int(i) <= center);
		
		photospline::bsplvb_simple(knots, n_knots, x, center,
					               order+1, localbasis_bsplvb.data());
		
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			photospline::bspline(knots, x, center + offset, order);
		}
		
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset],
							std::numeric_limits<float>::epsilon());
		}
	}
	
	// Within the support.
	for (size_t i = order+1; i < n_knots-order-1; i++) {
		x = (knots[i]+knots[i+1])/2.0;
		center = i;
		
		photospline::bsplvb_simple(knots, n_knots, x, center,
					               order+1, localbasis_bsplvb.data());
		
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			photospline::bspline(knots, x, center + offset, order);
		}
		
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset],
							std::numeric_limits<float>::epsilon());
		}
	}
	
	// After the last first fully-supported knot.
	for (size_t i = n_knots-order-1; i < n_knots-2; i++) {
		x = (knots[i]+knots[i+1])/2.0; // Less than fully-supported
		center = n_knots-order-2; // Last fully-supported spline.
		
		ENSURE(int(i) >= center);
		
		photospline::bsplvb_simple(knots, n_knots, x, center,
					               order+1, localbasis_bsplvb.data());
		
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			photospline::bspline(knots, x, center + offset, order);
		}
		
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset],
							std::numeric_limits<float>::epsilon());
		}
	}
}

TEST(bspline_deriv_nonzero_vs_bspline_deriv){
	const size_t n_knots = 10;
	const uint32_t order = 2;
	double x;
	double* knots;
	int center, offset;
	std::vector<float> localbasis_bsplvb(order+1);
	std::vector<float> localbasis_bspline(order+1);
	std::vector<double> knotvec;
	// This calculation is less stable.
	float tol = 10*std::numeric_limits<float>::epsilon();
	
	std::mt19937 rng;
	rng.seed(242);
	
	// Generate a random knot field.
	{
		std::uniform_real_distribution<> uniform(0,1);
		for (size_t i = 0; i < n_knots; i++)
			knotvec.push_back(uniform(rng));
	}
	std::sort(knotvec.begin(), knotvec.end());
	knots = knotvec.data();
	
	// Before the first fully-supported knot.
	for (size_t i = 0; i < order+1; i++) {
		x = (knots[i]+knots[i+1])/2.0; // Less than fully-supported
		center = order; // First fully-supported spline.
		
		ENSURE(int(i) <= center);
		
		photospline::bspline_deriv_nonzero(knots, n_knots, x, center,
							               order, localbasis_bsplvb.data());
		
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			photospline::bspline_deriv(knots, x, center + offset, order);
		}
		
		for (offset = 0; offset < order+1; offset++)
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
	}
	
	// Within the support.
	for (size_t i = order+1; i < n_knots-order-1; i++) {
		x = (knots[i]+knots[i+1])/2.0;
		center = i;
		
		photospline::bspline_deriv_nonzero(knots, n_knots, x, center,
							               order, localbasis_bsplvb.data());
		
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			photospline::bspline_deriv(knots, x, center + offset, order);
		}
		
		for (offset = 0; offset < order+1; offset++)
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
	}
	
	// After the last first fully-supported knot.
	for (size_t i = n_knots-order-1; i < n_knots-2; i++) {
		x = (knots[i]+knots[i+1])/2.0; // Less than fully-supported
		center = n_knots-order-2; // Last fully-supported spline.
		
		ENSURE(int(i) >= center);
		
		photospline::bspline_deriv_nonzero(knots, n_knots, x, center,
							               order, localbasis_bsplvb.data());
		
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			photospline::bspline_deriv(knots, x, center + offset, order);
		}
		
		for (offset = 0; offset < order+1; offset++)
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
	}
}

TEST(bspline_nonzero_vs_bspline){
	const size_t n_knots = 10;
	const int order = 2;
	double x;
	double* knots;
	int center, offset;
	std::vector<float> localbasis_bsplvb(order+1);
	std::vector<float> localbasis_bspline(order+1);
	std::vector<float> localbasis_bsplvb_deriv(order+1);
	std::vector<float> localbasis_bspline_deriv(order+1);
	std::vector<double> knotvec;
	// This calculation is less stable.
	float tol = 10*std::numeric_limits<float>::epsilon();
	
	std::mt19937 rng;
	rng.seed(342);
	
	// Generate a random knot field.
	{
		std::uniform_real_distribution<> uniform(0,1);
		for (size_t i = 0; i < n_knots; i++)
			knotvec.push_back(uniform(rng));
	}
	std::sort(knotvec.begin(), knotvec.end());
	knots = knotvec.data();
	
	// Before the first fully-supported knot.
	for (size_t i = 0; i < order+1; i++) {
		x = (knots[i]+knots[i+1])/2.0; // Less than fully-supported
		center = order; // First fully-supported spline.
		
		ENSURE(int(i) <= center);
		
		photospline::bspline_nonzero(knots, n_knots, x, center,
						order, localbasis_bsplvb.data(), localbasis_bsplvb_deriv.data());
		
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			photospline::bspline(knots, x, center + offset, order);
			localbasis_bspline_deriv[offset+order] =
			photospline::bspline_deriv(knots, x, center + offset, order);
		}
		
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
			ENSURE_DISTANCE(localbasis_bspline_deriv[offset], localbasis_bsplvb_deriv[offset], tol);
		}
	}
	
	// Within the support.
	for (size_t i = order+1; i < n_knots-order-1; i++) {
		x = (knots[i]+knots[i+1])/2.0;
		center = i;
		
		photospline::bspline_nonzero(knots, n_knots, x, center,
						order, localbasis_bsplvb.data(), localbasis_bsplvb_deriv.data());
		
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			photospline::bspline(knots, x, center + offset, order);
			localbasis_bspline_deriv[offset+order] =
			photospline::bspline_deriv(knots, x, center + offset, order);
		}
		
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
			ENSURE_DISTANCE(localbasis_bspline_deriv[offset], localbasis_bsplvb_deriv[offset], tol);
		}
	}
	
	// After the last first fully-supported knot.
	for (size_t i = n_knots-order-1; i < n_knots-2; i++) {
		x = (knots[i]+knots[i+1])/2.0; // Less than fully-supported
		center = n_knots-order-2; // Last fully-supported spline.
		
		ENSURE(int(i) >= center);
		
		photospline::bspline_nonzero(knots, n_knots, x, center,
						order, localbasis_bsplvb.data(), localbasis_bsplvb_deriv.data());
		
		for (offset = -order; offset <= 0; offset++) {
			ENSURE(offset+order >= 0);
			localbasis_bspline[offset+order] =
			photospline::bspline(knots, x, center + offset, order);
			localbasis_bspline_deriv[offset+order] =
			photospline::bspline_deriv(knots, x, center + offset, order);
		}
		
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_DISTANCE(localbasis_bspline[offset], localbasis_bsplvb[offset], tol);
			ENSURE_DISTANCE(localbasis_bspline_deriv[offset], localbasis_bsplvb_deriv[offset], tol);
		}
	}
}

// bspline_nonzero() gives the same result as bsplvb_simple(), ergo the
// value bases used in ndsplineeval() and ndsplineeval_gradient() are
// identical. Roundoff error in the derivative basis calculation
// scales with the spline order, so we use a large number here.

TEST(single_basis_vs_multi){
	const size_t n_knots = 100;
	const int order = 5;
	double x, *knots;
	int center, offset;
	std::vector<float> localbasis_bsplvb_simple(order+1);
	std::vector<float> localbasis_bspline_nonzero(order+1);
	std::vector<float> localbasis_bsplvb_deriv(order+1);
	std::vector<float> localbasis_bspline_nonzero_deriv(order+1);
	std::vector<float> localbasis_bspline_deriv_nonzero(order+1);
	// bsplvb() may access up to *order* elements off either end
	// Pad accordingly.
	std::vector<double> knotvec(n_knots + 2*order, 0.);
	
	std::mt19937 rng;
	rng.seed(442);
	
	// Generate a random knot field.
	{
		std::uniform_real_distribution<> uniform(0,1);
		for (size_t i = 0; i < n_knots; i++)
			knotvec.push_back(uniform(rng));
	}
	std::sort(knotvec.begin(), knotvec.end());
	knots = knotvec.data();
	
	// Before the first fully-supported knot.
	for (size_t i = 0; i < order+1; i++) {
		x = (knots[i]+knots[i+1])/2.0; // Less than fully-supported
		center = order; // First fully-supported spline.
		
		ENSURE(int(i) <= center);
		
		// As used in ndssplineeval_gradient()
		photospline::bspline_nonzero(knots, n_knots, x, center,
						order, localbasis_bspline_nonzero.data(), localbasis_bspline_nonzero_deriv.data());
		// As used in ndsplineeval()
		photospline::bsplvb_simple(knots, n_knots, x, center,
					  order + 1, localbasis_bsplvb_simple.data());
		photospline::bspline_deriv_nonzero(knots, n_knots, x, center,
							  order, localbasis_bspline_deriv_nonzero.data());
		
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_EQUAL(localbasis_bspline_nonzero[offset], localbasis_bsplvb_simple[offset]);
			ENSURE_EQUAL(localbasis_bspline_nonzero_deriv[offset], localbasis_bspline_deriv_nonzero[offset]);
		}
	}
	
	// Within the support.
	for (size_t i = order+1; i < n_knots-order-1; i++) {
		x = (knots[i]+knots[i+1])/2.0;
		center = i;
		
		photospline::bspline_nonzero(knots, n_knots, x, center,
						order, localbasis_bspline_nonzero.data(), localbasis_bspline_nonzero_deriv.data());
		photospline::bsplvb_simple(knots, n_knots, x, center,
					  order + 1, localbasis_bsplvb_simple.data());
		photospline::bspline_deriv_nonzero(knots, n_knots, x, center,
							  order, localbasis_bspline_deriv_nonzero.data());
		
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_EQUAL(localbasis_bspline_nonzero[offset], localbasis_bsplvb_simple[offset]);
			ENSURE_EQUAL(localbasis_bspline_nonzero_deriv[offset], localbasis_bspline_deriv_nonzero[offset]);
		}
	}
	
	// After the last first fully-supported knot.
	for (size_t i = n_knots-order-1; i < n_knots-2; i++) {
		x = (knots[i]+knots[i+1])/2.0; // Less than fully-supported
		center = n_knots-order-2; // Last fully-supported spline.
		
		ENSURE(int(i) >= center);
		
		photospline::bspline_nonzero(knots, n_knots, x, center,
						order, localbasis_bspline_nonzero.data(), localbasis_bspline_nonzero_deriv.data());
		photospline::bsplvb_simple(knots, n_knots, x, center,
					  order + 1, localbasis_bsplvb_simple.data());
		photospline::bspline_deriv_nonzero(knots, n_knots, x, center,
							  order, localbasis_bspline_deriv_nonzero.data());
		
		for (offset = 0; offset < order+1; offset++) {
			ENSURE_EQUAL(localbasis_bspline_nonzero[offset], localbasis_bsplvb_simple[offset]);
			ENSURE_EQUAL(localbasis_bspline_nonzero_deriv[offset], localbasis_bspline_deriv_nonzero[offset]);
		}
	}
}

//hammer on evaluation a bit
TEST(evaluation_benchmark){
	//TODO: include dimension 1 once it is implemented
	for(size_t dim=2; dim<6; dim++){
		photospline::splinetable<> spline("test_data/test_spline_"+std::to_string(dim)+"d.fits");
		std::cout << "Dimension " << spline.get_ndim() << " spline:" << std::endl;
		spline.benchmark_evaluation((unsigned int)(8.e6*exp(-(double)dim/1.4427)),true);
	}
	for(size_t dim=2; dim<6; dim++){
		photospline::splinetable<> spline("test_data/test_spline_"+std::to_string(dim)+"d_nco.fits");
		std::cout << "Dimension " << spline.get_ndim() << " (variable order) spline:" << std::endl;
		spline.benchmark_evaluation((unsigned int)(8.e6*exp(-(double)dim/1.4427)),true);
	}
}

TEST(permutation){
	const std::string splinePath="test_data/test_spline_4d_nco.fits";
	
	photospline::splinetable<> spline(splinePath);
	photospline::splinetable<>::evaluator evaluator=spline.get_evaluator();
	
	photospline::splinetable<> splineP(splinePath);
	splineP.permuteDimensions(std::vector<size_t>{2,3,0,1});
	photospline::splinetable<>::evaluator evaluatorP=splineP.get_evaluator();
	
	std::mt19937 rng;
	rng.seed(93);
	
	//Create uniform distributions over the support of the spline in all dimensions
	std::vector<std::uniform_real_distribution<>> dists;
	for(size_t i=0; i<spline.get_ndim(); i++)
		dists.push_back(std::uniform_real_distribution<>(spline.lower_extent(i),spline.upper_extent(i)));
	
	std::vector<double> coords(spline.get_ndim()), coordsP(spline.get_ndim());
	std::vector<int> centers(spline.get_ndim()), centersP(spline.get_ndim());
	for(size_t i=0; i<10000; i++) {
		for(size_t j=0; j<spline.get_ndim(); j++)
			coords[j]=dists[j](rng);
		coordsP[0]=coords[2];
		coordsP[1]=coords[3];
		coordsP[2]=coords[0];
		coordsP[3]=coords[1];
		
		ENSURE(evaluator.searchcenters(coords.data(), centers.data()), "Center lookup should succeed");
		ENSURE(evaluatorP.searchcenters(coordsP.data(), centersP.data()), "Center lookup should succeed");
		
		double evaluate=evaluator.ndsplineeval(coords.data(), centers.data(), 0);
		double evaluateP=evaluatorP.ndsplineeval(coordsP.data(), centersP.data(), 0);
		
		//Reordering the dimensions changes how the calculation rounds off during
		//intermediate steps, so fairly generour error tolerances are needed here.
		ENSURE_DISTANCE(evaluate,evaluateP,std::max(std::abs(evaluate*1e-4),1e-4),"Permuted spline should give same result for permuted coordinates");
	}
}