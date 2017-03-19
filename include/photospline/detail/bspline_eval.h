#ifndef PHOTOSPLINE_BSPLINE_EVAL_H
#define PHOTOSPLINE_BSPLINE_EVAL_H

#include <random>
#include <chrono>

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
		uint32_t max = nknots[i]-2;
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
		if (centers[i] == naxes[i])
			centers[i]--;
	}
	
	return (true);
}

template<typename Alloc>
double splinetable<Alloc>::ndsplineeval_core(const int* centers, int maxdegree, detail::buffer2d<float> localbasis) const
{
	uint32_t n;
	float basis_tree[ndim+1];
	int decomposedposition[ndim];
	
	int64_t tablepos = 0;
	for (n = 0; n < ndim; n++) {
		decomposedposition[n] = 0;
		tablepos += (centers[n] - (int64_t)order[n])*(int64_t)strides[n];
	}
	
	basis_tree[0] = 1;
	for (n = 0; n < ndim; n++)
		basis_tree[n+1] = basis_tree[n]*localbasis[n][0];
	uint32_t nchunks = 1;
	for (n = 0; n < ndim - 1; n++)
		nchunks *= (order[n] + 1);
	
	float result = 0;
	n = 0;
	while (true) {
		for (uint32_t i = 0; __builtin_expect(i < order[ndim-1] + 1, 1); i++)
			result+=basis_tree[ndim-1]*localbasis[ndim-1][i]*coefficients[tablepos + i];
		
		if (__builtin_expect(++n == nchunks, 0))
			break;
		
		tablepos += strides[ndim-2];
		decomposedposition[ndim-2]++;
		
		// Carry to higher dimensions
		uint32_t i;
		for (i = ndim-2;
			decomposedposition[i] > order[i]; i--) {
			decomposedposition[i-1]++;
			tablepos += (strides[i-1] - decomposedposition[i]*strides[i]);
			decomposedposition[i] = 0;
		}
		for (uint32_t j = i; __builtin_expect(j < ndim-1, 1); j++)
			basis_tree[j+1] = basis_tree[j]*
			localbasis[j][decomposedposition[j]];
	}
	
	return result;
}

template<typename Alloc>
template<unsigned int D>
double splinetable<Alloc>::ndsplineeval_coreD(const int* centers, int maxdegree, detail::buffer2d<float> localbasis) const
{
	uint32_t n;
	float basis_tree[D+1];
	int decomposedposition[D];
	
	int64_t tablepos = 0;
	for (n = 0; n < D; n++) {
		decomposedposition[n] = 0;
		tablepos += (centers[n] - (int64_t)order[n])*(int64_t)strides[n];
	}
	
	basis_tree[0] = 1;
	for (n = 0; n < D; n++)
		basis_tree[n+1] = basis_tree[n]*localbasis[n][0];
	uint32_t nchunks = 1;
	for (n = 0; n < D - 1; n++)
		nchunks *= (order[n] + 1);
	
	float result = 0;
	for(uint32_t n=0; __builtin_expect(n<(nchunks-1),1); n++){
		for (uint32_t i = 0; __builtin_expect(i < order[D-1] + 1, 1); i++)
			result+=basis_tree[D-1]*localbasis[D-1][i]*coefficients[tablepos + i];
		
		tablepos += strides[D-2u];
#ifdef __clang__ //this code is unreachable if D<2
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Warray-bounds"
#endif
		decomposedposition[D-2u]++;
#ifdef __clang__
	#pragma clang diagnostic pop
#endif
		
		// Carry to higher dimensions
		uint32_t i;
		for (i = D-2; decomposedposition[i] > order[i]; i--) {
			decomposedposition[i-1]++;
			tablepos += (strides[i-1] - decomposedposition[i]*strides[i]);
			decomposedposition[i] = 0;
		}
		for (uint32_t j = i; __builtin_expect(j < D-1, 1); j++)
			basis_tree[j+1] = basis_tree[j]*localbasis[j][decomposedposition[j]];
	}
	for (uint32_t i = 0; i < (order[D-1] + 1); i++)
		result+=basis_tree[D-1]*localbasis[D-1][i]*coefficients[tablepos+i];
	
	return result;
}

template<typename Alloc>
template<unsigned int D, unsigned int O>
double splinetable<Alloc>::ndsplineeval_coreD_FixedOrder(const int* centers, int maxdegree, detail::buffer2d<float> localbasis) const
{
	uint32_t n;
	float basis_tree[D+1];
	int decomposedposition[D];
	
	int64_t tablepos = 0;
	for (n = 0; n < D; n++) {
		decomposedposition[n] = 0;
		tablepos += (centers[n] - (int64_t)O)*(int64_t)strides[n];
	}
	
	basis_tree[0] = 1;
	for (n = 0; n < D; n++)
		basis_tree[n+1] = basis_tree[n]*localbasis[n][0];
	uint32_t nchunks = 1;
	for (n = 0; n < D - 1; n++)
		nchunks *= (O + 1);
	
	float result = 0;
	for(uint32_t n=0; __builtin_expect(n<(nchunks-1),1); n++){
		for (uint32_t i = 0; i < (O + 1); i++)
			result+=basis_tree[D-1]*localbasis[D-1][i]*coefficients[tablepos + i];
		
		tablepos += strides[D-2u];
#ifdef __clang__ //this code is unreachable if D<2
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Warray-bounds"
#endif
		decomposedposition[D-2u]++;
#ifdef __clang__
	#pragma clang diagnostic pop
#endif
		
		// Carry to higher dimensions
		uint32_t i;
		for (i = D-2; decomposedposition[i] > O; i--) {
			decomposedposition[i-1]++;
			tablepos += (strides[i-1] - decomposedposition[i]*strides[i]);
			decomposedposition[i] = 0;
		}
		for (uint32_t j = i; __builtin_expect(j < D-1, 1); j++)
			basis_tree[j+1] = basis_tree[j]*localbasis[j][decomposedposition[j]];
	}
	for (uint32_t i = 0; i < (O + 1); i++)
		result+=basis_tree[D-1]*localbasis[D-1][i]*coefficients[tablepos+i];
	
	return result;
}

namespace detail {

template<unsigned int O1>
constexpr unsigned int nchunks()
{
	return 1u;
}

// product of order+1 of the up to last dimension
template<unsigned int O1, unsigned int O2, unsigned int ... Orders>
constexpr unsigned int nchunks()
{
	return (O1+1)*nchunks<O2, Orders...>();
}

template<unsigned int O1>
constexpr unsigned int chunk()
{
	return O1+1;
}

// order+1 in last dimension
template<unsigned int O1, unsigned int O2, unsigned int ... Orders>
constexpr unsigned int chunk()
{
	return chunk<O2, Orders...>();
}

template<typename Alloc>
bool orders_are(const splinetable<Alloc> &spline, const std::initializer_list<unsigned int> &orders)
{
	if (orders.size() != spline.get_ndim())
		return false;
	unsigned int i = 0;
	for(auto order : orders) {
		if (spline.get_order(i) != order)
			return false;
		i++;
	}
	return true;
}

}

template<typename Alloc>
template<unsigned int ... Orders>
double splinetable<Alloc>::ndsplineeval_core_KnownOrder(const int* centers, int maxdegree, detail::buffer2d<float> localbasis) const
{
	constexpr unsigned int D = sizeof...(Orders);
	uint32_t n;
	float basis_tree[D+1];
	int decomposedposition[D];
	
	// 
	int64_t tablepos = 0;
	for (n = 0; n < D; n++) {
		decomposedposition[n] = 0;
		tablepos += (centers[n] - (int64_t)order[n])*(int64_t)strides[n];
	}
	
	basis_tree[0] = 1;
	for (n = 0; n < D; n++)
		basis_tree[n+1] = basis_tree[n]*localbasis[n][0];
	
	constexpr uint32_t nchunks = detail::nchunks<Orders...>();
	constexpr uint32_t chunk = detail::chunk<Orders...>();

	float result = 0;
	for(uint32_t n=0; __builtin_expect(n<(nchunks-1),1); n++){
		for (uint32_t i = 0; __builtin_expect(i < chunk, 1); i++)
			result+=basis_tree[D-1]*localbasis[D-1][i]*coefficients[tablepos + i];
		
		tablepos += strides[D-2u];
#ifdef __clang__ //this code is unreachable if D<2
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Warray-bounds"
#endif
		decomposedposition[D-2u]++;
#ifdef __clang__
	#pragma clang diagnostic pop
#endif
		
		// Carry to higher dimensions
		uint32_t i;
		for (i = D-2; decomposedposition[i] > order[i]; i--) {
			decomposedposition[i-1]++;
			tablepos += (strides[i-1] - decomposedposition[i]*strides[i]);
			decomposedposition[i] = 0;
		}
		for (uint32_t j = i; __builtin_expect(j < D-1, 1); j++)
			basis_tree[j+1] = basis_tree[j]*localbasis[j][decomposedposition[j]];
	}
	for (uint32_t i = 0; i < chunk; i++)
		result+=basis_tree[D-1]*localbasis[D-1][i]*coefficients[tablepos+i];
	
	return result;
}

template<typename Alloc>
double splinetable<Alloc>::ndsplineeval(const double* x, const int* centers, int derivatives) const
{
	uint32_t maxdegree = *std::max_element(order,order+ndim) + 1;
	float localbasis_store[ndim*maxdegree];
	detail::buffer2d<float> localbasis{localbasis_store,maxdegree};
	
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
			for (int32_t i = 0; i <= order[n]; i++)
				localbasis[n][i] = bspline_deriv(
												   &knots[n][0], x[n],
												   centers[n] - order[n] + i, 
												   order[n], derivatives[n]);
		}
	}
	
	return ndsplineeval_core(centers, maxdegree, localbasis);
}
	
template<typename Alloc>
typename splinetable<Alloc>::evaluator
splinetable<Alloc>::get_evaluator() const{
	evaluator eval(*this);
	
	uint32_t constOrder = order[0];
	for (unsigned int j = 1; j < ndim; j++) {
		if (order[j] != constOrder) {
			constOrder = 0;
			break;
		}
	}
	
	switch(constOrder){
#ifdef PHOTOSPLINE_EVAL_TEMPLATES
		case 2:
			switch(ndim){
				case 1:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD_FixedOrder<1,2>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD_FixedOrder<1,2>;
					break;
				case 2:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD_FixedOrder<2,2>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD_FixedOrder<2,2>;
					break;
				case 3:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD_FixedOrder<3,2>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD_FixedOrder<3,2>;
					break;
				case 4:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD_FixedOrder<4,2>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD_FixedOrder<4,2>;
					break;
				case 5:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD_FixedOrder<5,2>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD_FixedOrder<5,2>;
					break;
				case 6:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD_FixedOrder<6,2>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD_FixedOrder<6,2>;
					break;
				case 7:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD_FixedOrder<7,2>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD_FixedOrder<7,2>;
					break;
				case 8:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD_FixedOrder<8,2>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD_FixedOrder<8,2>;
					break;
				default:
					eval.eval_ptr=&splinetable::ndsplineeval_core;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_core;
			}
			break;
		case 3:
			switch(ndim){
				case 1:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD_FixedOrder<1,3>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD_FixedOrder<1,3>;
					break;
				case 2:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD_FixedOrder<2,3>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD_FixedOrder<2,3>;
					break;
				case 3:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD_FixedOrder<3,3>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD_FixedOrder<3,3>;
					break;
				case 4:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD_FixedOrder<4,3>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD_FixedOrder<4,3>;
					break;
				case 5:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD_FixedOrder<5,3>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD_FixedOrder<5,3>;
					break;
				case 6:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD_FixedOrder<6,3>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD_FixedOrder<6,3>;
					break;
				case 7:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD_FixedOrder<7,3>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD_FixedOrder<7,3>;
					break;
				case 8:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD_FixedOrder<8,3>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD_FixedOrder<8,3>;
					break;
				default:
					eval.eval_ptr=&splinetable::ndsplineeval_core;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_core;
			}
			break;
#endif
		default:
			switch(ndim){
#ifdef PHOTOSPLINE_EVAL_TEMPLATES
				case 1:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD<1>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD<1>; break;
					break;
				case 2:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD<2>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD<2>; break;
					break;
				case 3:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD<3>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD<3>; break;
					break;
				case 4:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD<4>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD<4>; break;
					break;
				case 5:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD<5>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD<5>; break;
					break;
				case 6:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD<6>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD<6>; break;
					break;
				case 7:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD<7>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD<7>; break;
					break;
				case 8:
					eval.eval_ptr=&splinetable::ndsplineeval_coreD<8>;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD<8>; break;
					break;
#endif
				default:
					eval.eval_ptr=&splinetable::ndsplineeval_core;
					eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_core;
			}
	}

#ifdef PHOTOSPLINE_EVAL_TEMPLATES
	// Mixed orders known to exist in the wild
	if (detail::orders_are(*this, {2,2,2,3,2,2})) {
		eval.eval_ptr=&splinetable::ndsplineeval_core_KnownOrder<2,2,2,3,2,2>;
		eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD<6>;
		eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_core_KnownOrder<2,2,2,3,2,2>;
		
	} else if (detail::orders_are(*this, {2,2,2,5,2,2})) {
		eval.eval_ptr=&splinetable::ndsplineeval_core_KnownOrder<2,2,2,5,2,2>;
		eval.v_eval_ptr=&splinetable::ndsplineeval_multibasis_coreD<6>;
	}
#endif

	return(eval);
}

template<typename Alloc>
bool splinetable<Alloc>::evaluator::searchcenters(const double* x, int* centers) const{
	return(table.searchcenters(x,centers));
}
	
template<typename Alloc>
double splinetable<Alloc>::evaluator::ndsplineeval(const double* x, const int* centers, int derivatives) const{
	uint32_t maxdegree = *std::max_element(table.order,table.order+table.ndim) + 1;
	float localbasis_store[table.ndim*maxdegree];
	detail::buffer2d<float> localbasis{localbasis_store,maxdegree};
	
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
double splinetable<Alloc>::evaluator::operator()(const double* x, int derivatives) const{
	int centers[table.ndim];
	if(!table.searchcenters(x,centers))
		return(0);
	return(ndsplineeval(x,centers,derivatives));
}
	
template<typename Alloc>
double splinetable<Alloc>::evaluator::ndsplineeval_deriv(const double* x, const int* centers, const unsigned int *derivatives) const
{
	uint32_t maxdegree = *std::max_element(table.order,table.order+table.ndim) + 1;
	float localbasis_store[table.ndim*maxdegree];
	detail::buffer2d<float> localbasis{localbasis_store,maxdegree};
	
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
			for (int32_t i = 0; i <= table.order[n]; i++)
				localbasis[n][i] = bspline_deriv(
												   &table.knots[n][0], x[n],
												   centers[n] - table.order[n] + i, 
												   table.order[n], derivatives[n]);
		}
	}
	
	return((table.*(eval_ptr))(centers, maxdegree, localbasis));
}
	
template<typename Alloc>
typename splinetable<Alloc>::benchmark_results
splinetable<Alloc>::benchmark_evaluation(size_t trialCount, bool verbose){
	std::default_random_engine rng;
	
	volatile double dummy;
	benchmark_results result;
	evaluator eval=get_evaluator();
	
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
	
	return(result);
}

} //namespace photospline

#endif
