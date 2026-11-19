#ifndef PHOTOSPLINE_BSPLINE_EVAL_H
#define PHOTOSPLINE_BSPLINE_EVAL_H

#include <random>
#include <chrono>
#include <limits>

namespace photospline{
	
template<typename Alloc>
bool splinetable<Alloc>::searchcenters(const double* x, int* centers) const
{
	for (uint32_t i = 0; i < ndim; i++) {
		
		/* Ensure we are actually inside the table. */
		if (x[i] <= knots[i][0] ||
			x[i] > knots[i][nknots[i]-1])
			return (false);
		
		/*
		 * If we're only a few knots in, take the center to be
		 * the nearest fully-supported knot.
		 */
		if (x[i] < knots[i][order[i]]) {
			centers[i] = order[i];
			continue;
		} else if (x[i] >= knots[i][naxes[i]]) {
			centers[i] = naxes[i]-1;
			continue;
		}
		
		uint32_t min = order[i];
		uint32_t max = naxes[i];
		double diff = x[i] - knots[i][min];
		uint32_t hi = (std::numeric_limits<uint32_t>::max() < diff*rmin_sep[i] ?
									 std::numeric_limits<uint32_t>::max() : static_cast<uint32_t>(diff*rmin_sep[i]));
		uint32_t lo = static_cast<uint32_t>(diff*rmax_sep[i]);
		if (hi < hi+min+1 && hi+min+1 < max) max = hi+min+1;
		min += lo;
		do {
			centers[i] = (max+min)/2;
			
			if (x[i] < knots[i][centers[i]])
				max = centers[i]-1;
			else
				min = centers[i]+1;
		} while (x[i] < knots[i][centers[i]] ||
				 x[i] >= knots[i][centers[i]+1]);
		
		/*
		 * B-splines are defined on a half-open interval. For the
		 * last point of the interval, move center one point to the
		 * left to get the limit of the sum without evaluating
		 * absent basis functions.
		 */
		if (unsigned(centers[i]) == naxes[i])
			centers[i]--;
	}
	
	return (true);
}

template<typename Alloc>
template<typename Float, unsigned int Row>
double splinetable<Alloc>::ndsplineeval_core_row(const int* centers, int maxdegree __attribute__((unused)), detail::buffer2d<Float> localbasis) const
{
	Float result = 0;
	detail::contract<Row, 1>(&coefficients[0], &strides[0], &order[0], ndim, centers,
	    [&](unsigned d, unsigned i){ return &localbasis[d][i]; }, &result);
	return result;
}

template<typename Alloc>
template<typename Float>
typename splinetable<Alloc>::template core_ptr<Float> splinetable<Alloc>::select_core() const
{
	switch (order[ndim-1] + 1) {
#ifndef PHOTOSPLINE_NO_EVAL_TEMPLATES
		case 2: return &splinetable::ndsplineeval_core_row<Float, 2>;
		case 3: return &splinetable::ndsplineeval_core_row<Float, 3>;
		case 4: return &splinetable::ndsplineeval_core_row<Float, 4>;
		case 5: return &splinetable::ndsplineeval_core_row<Float, 5>;
		case 6: return &splinetable::ndsplineeval_core_row<Float, 6>;
#endif
		default: return &splinetable::ndsplineeval_core_row<Float, 0>;
	}
}

template<typename Alloc>
template<typename Float>
double splinetable<Alloc>::ndsplineeval_core(const int* centers, int maxdegree, detail::buffer2d<Float> localbasis) const
{
	return (this->*select_core<Float>())(centers, maxdegree, localbasis);
}

template<typename Alloc>
template<typename Float>
double splinetable<Alloc>::ndsplineeval(const double* x, const int* centers, int derivatives) const
{
	uint32_t maxdegree = *std::max_element(order,order+ndim) + 1;
	Float localbasis_store[ndim*maxdegree];
	detail::buffer2d<Float> localbasis{localbasis_store,maxdegree};
	
	for (uint32_t n = 0; n < ndim; n++) {
		if (derivatives & (1 << n)) {
			bspline_deriv_nonzero(&knots[n][0],
								  nknots[n], x[n], centers[n],
								  order[n], localbasis[n]);
		} else {
			bsplvb_simple(&knots[n][0], nknots[n],
						  x[n], centers[n], order[n] + 1,
						  localbasis[n]);
		}
	}
	
	return(ndsplineeval_core(centers, maxdegree, localbasis));
}
	
template<typename Alloc>
double splinetable<Alloc>::operator()(const double* x) const{
	int centers[ndim];
	if(!searchcenters(x,centers))
		return(0);
	return(ndsplineeval(x,centers,0));
}

template<typename Alloc>
double splinetable<Alloc>::ndsplineeval_deriv(const double* x, const int* centers, const unsigned int *derivatives) const
{
	uint32_t maxdegree = *std::max_element(order,order+ndim) + 1;
	float localbasis_store[ndim*maxdegree];
	detail::buffer2d<float> localbasis{localbasis_store,maxdegree};
	
	for (uint32_t n = 0; n < ndim; n++) {
		if (derivatives == nullptr || derivatives[n] == 0) {
			bsplvb_simple(&knots[n][0], nknots[n],
						  x[n], centers[n], order[n] + 1,
						  localbasis[n]);
		} else if (derivatives[n] == 1) {
			bspline_deriv_nonzero(&knots[n][0], nknots[n],
						  x[n], centers[n], order[n],
						  localbasis[n]);
		} else {
			for (uint32_t i = 0; i <= order[n]; i++)
				localbasis[n][i] = bspline_deriv(
												   &knots[n][0], x[n],
												   centers[n] - order[n] + i, 
												   order[n], derivatives[n]);
		}
	}
	
	return ndsplineeval_core(centers, maxdegree, localbasis);
}
	
template<typename Alloc>
template<typename Float>
typename splinetable<Alloc>::template evaluator_type<Float>
splinetable<Alloc>::get_evaluator() const{
	evaluator_type<Float> eval(*this);
	eval.eval_ptr=select_core<Float>();
	eval.v_eval_ptr=select_multibasis_core<Float>();
	return(eval);
}

template<typename Alloc>
template<typename Float>
bool splinetable<Alloc>::evaluator_type<Float>::searchcenters(const double* x, int* centers) const{
	return(table.searchcenters(x,centers));
}
	
template<typename Alloc>
template<typename Float>
double splinetable<Alloc>::evaluator_type<Float>::ndsplineeval(const double* x, const int* centers, int derivatives) const{
	uint32_t maxdegree = *std::max_element(table.order,table.order+table.ndim) + 1;
	Float localbasis_store[table.ndim*maxdegree];
	detail::buffer2d<Float> localbasis{localbasis_store,maxdegree};
	
	for (uint32_t n = 0; n < table.ndim; n++) {
		if (derivatives & (1 << n)) {
			bspline_deriv_nonzero(&table.knots[n][0],
								  table.nknots[n], x[n], centers[n],
								  table.order[n], localbasis[n]);
		} else {
			bsplvb_simple(&table.knots[n][0], table.nknots[n],
						  x[n], centers[n], table.order[n] + 1,
						  localbasis[n]);
		}
	}
	
	return((table.*(eval_ptr))(centers, maxdegree, localbasis));
}
	
template<typename Alloc>
template<typename Float>
double splinetable<Alloc>::evaluator_type<Float>::operator()(const double* x, int derivatives) const{
	int centers[table.ndim];
	if(!table.searchcenters(x,centers))
		return(0);
	return(ndsplineeval(x,centers,derivatives));
}
	
template<typename Alloc>
template<typename Float>
double splinetable<Alloc>::evaluator_type<Float>::ndsplineeval_deriv(const double* x, const int* centers, const unsigned int *derivatives) const
{
	uint32_t maxdegree = *std::max_element(table.order,table.order+table.ndim) + 1;
	Float localbasis_store[table.ndim*maxdegree];
	detail::buffer2d<Float> localbasis{localbasis_store,maxdegree};
	
	for (uint32_t n = 0; n < table.ndim; n++) {
		if (derivatives == nullptr || derivatives[n] == 0) {
			bsplvb_simple(&table.knots[n][0], table.nknots[n],
						  x[n], centers[n], table.order[n] + 1,
						  localbasis[n]);
		} else if (derivatives[n] == 1) {
			bspline_deriv_nonzero(&table.knots[n][0], table.nknots[n],
						  x[n], centers[n], table.order[n],
						  localbasis[n]);
		} else {
			for (uint32_t i = 0; i <= table.order[n]; i++)
				localbasis[n][i] = bspline_deriv(
												   &table.knots[n][0], x[n],
												   centers[n] - table.order[n] + i, 
												   table.order[n], derivatives[n]);
		}
	}
	
	return((table.*(eval_ptr))(centers, maxdegree, localbasis));
}
	
template<typename Alloc>
template<typename Float>
typename splinetable<Alloc>::benchmark_results
splinetable<Alloc>::benchmark_evaluation(size_t trialCount, bool verbose){
	std::default_random_engine rng;
	
	volatile double dummy;
	benchmark_results result;
	evaluator_type<Float> eval=get_evaluator<Float>();
	
	std::vector<std::uniform_real_distribution<>> dists;
	for(size_t i=0; i<ndim; i++)
		dists.push_back(std::uniform_real_distribution<>(lower_extent(i),upper_extent(i)));
	
	std::vector<double> coords(ndim);
	std::vector<int> centers(ndim);
	std::vector<double> gradeval(ndim+1);
	std::chrono::high_resolution_clock::time_point t1, t2;
	
	//estimate time wasted running the RNG
	rng.seed(52);
	t1 = std::chrono::high_resolution_clock::now();
	for(size_t i=0; i<trialCount; i++){
		for(size_t j=0; j<ndim; j++)
			dummy=dists[j](rng);
	}
	t2 = std::chrono::high_resolution_clock::now();
	double rng_overhead=std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();
	//printf("RNG overhead: %lf s\n",rng_overhead);
	
	rng.seed(52);
	t1 = std::chrono::high_resolution_clock::now();
	for(size_t i=0; i<trialCount; i++){
		for(size_t j=0; j<ndim; j++)
			coords[j]=dists[j](rng);
		
		if(!searchcenters(coords.data(), centers.data()))
			throw std::logic_error("center lookup failure for point which should be in bounds");
		
		dummy=eval.ndsplineeval(coords.data(), centers.data(), 0);
	}
	t2 = std::chrono::high_resolution_clock::now();
	result.single_eval_rate=trialCount/
	(std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count()-rng_overhead);
	
	rng.seed(52);
	t1 = std::chrono::high_resolution_clock::now();
	for(size_t i=0; i<trialCount; i++){
		for(size_t j=0; j<ndim; j++)
			coords[j]=dists[j](rng);
		
		if(!searchcenters(coords.data(), centers.data()))
			throw std::logic_error("center lookup failure for point which should be in bounds");
		
		dummy=eval.ndsplineeval(coords.data(), centers.data(), 0);
		for(size_t j=0; j<ndim; j++)
			dummy=eval.ndsplineeval(coords.data(), centers.data(), 1u<<j);
	}
	t2 = std::chrono::high_resolution_clock::now();
	result.gradient_single_eval_rate=trialCount/
	(std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count()-rng_overhead);
	
	rng.seed(52);
	t1 = std::chrono::high_resolution_clock::now();
	for(size_t i=0; i<trialCount; i++){
		for(size_t j=0; j<ndim; j++)
			coords[j]=dists[j](rng);
		
		if(!searchcenters(coords.data(), centers.data()))
			throw std::logic_error("center lookup failure for point which should be in bounds");
		
		eval.ndsplineeval_gradient(coords.data(), centers.data(), gradeval.data());
	}
	t2 = std::chrono::high_resolution_clock::now();
	result.gradient_multi_eval_rate=trialCount/
	(std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count()-rng_overhead);
	
	if(verbose){
		printf("Benchmark results:\n");
		printf("\t%.2le single evaluations/second\n",result.single_eval_rate);
		printf("\t%.2le 'single' gradient evaluations/second\n",result.gradient_single_eval_rate);
		printf("\t%.2le 'multiple' gradient evaluations/second\n",result.gradient_multi_eval_rate);
		printf("\t(%zu trial evaluations)\n",trialCount);
	}
	(void)dummy;
	
	return(result);
}

} //namespace photospline

#endif
