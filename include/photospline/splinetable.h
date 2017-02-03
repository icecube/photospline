
#ifndef PHOTOSPLINE_SPLINETABLE_H
#define PHOTOSPLINE_SPLINETABLE_H

#include <algorithm>
#include <cassert>
#include <memory>
#include <numeric>
#include <sstream>
#include <vector>

#include "photospline/bspline.h"
#include "photospline/detail/simd.h"

#include <string.h>
#include <fitsio.h>
#include <fitsio2.h>

#ifdef PHOTOSPLINE_INCLUDES_SPGLAM
#include "photospline/detail/splineutil.h"
#include "photospline/detail/glam.h"
#endif //PHOTOSPLINE_INCLUDES_SPGLAM

#include <iostream>

namespace photospline{
	
#ifdef PHOTOSPLINE_INCLUDES_SPGLAM
///A more user-friendly version of the C ndsparse
///Note that this does not mean entirely safe; as all of the internals are still
///exposed, leaving direct modification of the structure possible.
struct ndsparse : public ::ndsparse{
	size_t entriesInserted;
	///\brief 'Raw' constructor which leaves the object essentially uninitialized
	ndsparse(){
		ndim=0;
		rows=0;
		i=nullptr;
		ranges=nullptr;
		x=nullptr;
		entriesInserted=0;
	}
	///\brief Proper constructor which allocates memory
	///\param rows the number of nonzero entries for which space will be allocated
	///\param ndim the number of dimensions
	ndsparse(size_t rows, size_t ndim){
		if(!ndim)
			throw std::logic_error("Tried to allocate an ndsparse with dimension 0");
		if(!rows)
			throw std::logic_error("Tried to allocate an ndsparse with 0 entries");
		if(ndsparse_allocate(this,rows,ndim)!=0)
			throw std::bad_alloc();
		entriesInserted=0;
	}
	~ndsparse(){
		ndsparse_free(this);
	}
	///\brief A convenience interface for adding entries.
	///It is a bad idea to mix use of this interface with direct alterations to
	///the x and i arrays.
	void insertEntry(double value, unsigned int* indices){
		//assert(entriesInserted+1<rows && "Attempt to insert and entry into a full ndsparse");
		if(!(entriesInserted<rows))
			throw std::runtime_error("Attempt to insert an entry into a full ndsparse");
		x[entriesInserted]=value;
		for(size_t j=0; j<ndim; j++){
			i[j][entriesInserted] = indices[j];
			ranges[j] = std::max(ranges[j],i[j][entriesInserted]+1);
		}
		entriesInserted++;
	}
};
#endif //PHOTOSPLINE_INCLUDES_SPGLAM

namespace detail{
	///Work-around to replace passing 2D VLAs
	template<typename T>
	struct buffer2d{
		T* buf;
		size_t dim1;
		buffer2d(T* b, size_t d1):buf(b),dim1(d1){}
		
		T* operator[](size_t idx1){
			return(buf+dim1*idx1); }
		const T* operator[](size_t idx1) const{ return(buf+dim1*idx1); }
	};
	///A simple view into un-owned data.
	template<typename T>
	struct array_view{
	public:
		using value_type=const T;
		using reference=const T&;
		using const_reference=const T&;
		using pointer=const T*;
		using const_pointer=const T*;
		using iterator=const T*;
		using const_iterator=const T*;
		using difference_type=typename std::iterator_traits<iterator>::difference_type;
		using size_type=size_t;
	private:
		const_pointer d;
		size_type s;
	public:
		array_view():d(nullptr),s(0){}
		array_view(const_pointer d, size_type s):d(d),s(s){}
		void reset(const_pointer d, size_type s){
			this->d=d;
			this->s=s;
		}
		size_type size() const{ return(s); }
		size_type max_size() const{ return(s); }
		bool empty() const{ return(s==0); }
		const_pointer data() const{ return(d); }
		const_reference operator[](size_type i) const{ return(*(d+i)); }
		iterator begin() const{ return(d); }
		const_iterator cbegin() const{ return(d); }
		iterator end() const{ return(d+s); }
		const_iterator cend() const{ return(d+s); }
	};
}
	
template<typename Alloc = std::allocator<void> >
class splinetable{
public:
	typedef Alloc allocator_type;
	typedef std::allocator_traits<allocator_type> allocator_traits;
	
	typedef typename allocator_traits::template rebind_traits<uint32_t>::pointer uint32_t_ptr;
	typedef typename allocator_traits::template rebind_traits<uint64_t>::pointer uint64_t_ptr;
	typedef typename allocator_traits::template rebind_traits<float>::pointer float_ptr;
	typedef typename allocator_traits::template rebind_traits<double>::pointer double_ptr;
	typedef typename allocator_traits::template rebind_traits<double_ptr>::pointer double_ptr_ptr;
	typedef typename allocator_traits::template rebind_traits<char>::pointer char_ptr;
	typedef typename allocator_traits::template rebind_traits<char_ptr>::pointer char_ptr_ptr;
	typedef typename allocator_traits::template rebind_traits<char_ptr_ptr>::pointer char_ptr_ptr_ptr;
	
	explicit splinetable(allocator_type alloc=Alloc()):
	ndim(0),order(NULL),knots(NULL),nknots(NULL),extents(NULL),periods(NULL),
	coefficients(NULL),naxes(NULL),strides(NULL),naux(0),aux(NULL),allocator(alloc)
	{}
	
	explicit splinetable(const std::string& filePath, allocator_type alloc=Alloc()):
	ndim(0),order(NULL),knots(NULL),nknots(NULL),extents(NULL),periods(NULL),
	coefficients(NULL),naxes(NULL),strides(NULL),naux(0),aux(NULL),allocator(alloc)
	{
		read_fits(filePath);
	}
	
	~splinetable(){
		if(ndim){
			uint64_t ncoeffs=strides[0]*naxes[0];
			for(uint32_t i=0; i<ndim; i++)
				deallocate(knots[i]-order[i],nknots[i]+2*order[i]);
			deallocate(knots,ndim);
			deallocate(nknots,ndim);
			deallocate(order,ndim);
			if(extents){
				deallocate(extents[0],2*ndim);
				deallocate(extents,ndim);
			}
			if(periods)
				deallocate(periods,ndim);
			deallocate(coefficients,ncoeffs);
			deallocate(naxes,ndim);
			deallocate(strides,ndim);
			for(uint32_t i=0; i<naux; i++){
				deallocate(aux[i][0],strlen(&aux[i][0][0])+1);
				deallocate(aux[i][1],strlen(&aux[i][1][0])+1);
				deallocate(aux[i],2);
			}
			deallocate(aux,naux);
		}
	}
	
	static size_t estimateMemory(const std::string& filePath,
	                             uint32_t n_convolution_knots = 1,
	                             uint32_t convolution_dimension = 0);
	
	///Read from a file
	bool read_fits(const std::string& path);
	
	///Read from a fits 'file' in a memory buffer
	bool read_fits_mem(void* buffer, size_t buffer_size);
	
	///Write to a file
	void write_fits(const std::string& path) const;
	
	///Write to a fits memory 'file'.
	///\returns A pair containing a buffer allocated by malloc() which should be
	///         deallocated with free() and the size of that buffer.
	std::pair<void*,size_t> write_fits_mem() const;
	
	///Get the count of 'auxiliary' keys
	size_t get_naux_values() const{ return(naux); }
	///Directly get a particular auxiliary key
	const char* get_aux_key(size_t i) const{ return(aux[i][0]); }
	///Directly get a particular auxiliary value
	const char* get_aux_value(const char* key) const;
	///Delete an auxiliary key and associated value
	///\returns true if the key existed and was removed
	bool remove_key(const char* key);
	///Look up the value associated with an auxiliary key
	///\param key the name of the key to look up
	///\param result location to store the value if the key is found
	template<typename T>
	bool read_key(const char* key, T& result) const;
	///Insert or overwrite an auxiliary key,value pair
	///\param key the name of the key to store
	///\param value the value to store for the key
	template<typename T>
	bool write_key(const char* key, const T& value);
	
	///\brief Manages use of optimized spline evaluation routines
	///
	///Particular splines may have properies which can be exploited to evaulate
	///them more efficiently (e.g. having a constant, low order in all
	///dimensions). For technical reasons (primarily the inadvisability of
	///storing function pointers in shared memory) it is best to separate the
	///data about which internal functions have been selected for greatest
	///efficiency from the splinetable object itself. The evaluator type
	///safely encapsulates this information, and presents the same evaluation
	///interface as the splietable itself to users.
	///
	///Since the evaluator holds a reference to the actual splinetable from
	///which is want obtained it must be considered invalidated if that table is
	///altered or destroyed.
	struct evaluator{
	private:
		const splinetable<Alloc>& table;
		double (splinetable::*eval_ptr)(const int*, int, detail::buffer2d<float>) const;
		void (splinetable::*v_eval_ptr)(const int*, const v4sf***, v4sf*) const;
		friend class splinetable<Alloc>;
		evaluator(const splinetable<Alloc>& table):table(table){}
	public:
		///\brief Get the underlying splinetable
		const splinetable<Alloc>& get_table() const{ return(table); }
		///\brief same as splinetable::searchcenters
		bool searchcenters(const double* x, int* centers) const;
		///\brief same as splinetable::ndsplineeval
		double ndsplineeval(const double* x, const int* centers, int derivatives) const;
		///\brief same as splinetable::ndsplineeval_gradient
		void ndsplineeval_gradient(const double* x, const int* centers, double* evaluates) const;
		///\brief same as splinetable::ndsplineeval_deriv2
		double ndsplineeval_deriv2(const double* x, const int* centers, int derivatives) const;
	};
	friend struct evaluator;
	
	///Constructs an optimized evaluator object which will use the best
	///available internal routines to perform evaulations. The evaluator holds
	///a reference to this splinetable, so it must be considered invalidated if
	///this table altered or destroyed.
	evaluator get_evaluator() const;
	
	/*
	 * Spline table based hypersurface evaluation. ndsplineeval() takes a spline
	 * coefficient table, a vector at which to evaluate the surface, and a vector
	 * indicating the evaluation centers, as for splineeval().
	 *
	 * searchcenters() provides a method to acquire a centers vector
	 * for ndsplineeval() using a binary search. Depending on how the table
	 * was produced, a more efficient method may be available.
	 */
	bool searchcenters(const double* x, int* centers) const;
	
	double ndsplineeval(const double* x, const int* centers, int derivatives) const;
	
	double ndsplineeval_deriv2(const double* x, const int* centers, int derivatives) const;
	
	/* Evaluate a spline surface and all its derivatives at x */
	void ndsplineeval_gradient(const double* x, const int* centers, double* evaluates) const;
	
	struct benchmark_results{
		double single_eval_rate;
		double gradient_single_eval_rate;
		double gradient_multi_eval_rate;
	};
	///Evaluate the spline at random points within the extent
	benchmark_results benchmark_evaluation(size_t trialCount=1e4, bool verbose=false);
	
	/*
	 * Convolve a table with the spline defined on a set of knots along a given
	 * dimension and store the spline expansion of the convolved surface in the
	 * table. This will raise the order of the splines in the given dimension by
	 * (n_knots - 1), i.e. convolving with an order-0 spline (a box function,
	 * defined on two knots) will raise the order of the spline surface by 1.
	 */
	void convolve(const uint32_t dim, const double* knots, size_t n_knots);
	
	uint32_t get_ndim() const{ return(ndim); }
	uint32_t get_order(uint32_t dim) const{
		assert(dim<ndim);
		return(order[dim]);
	}
	uint64_t get_nknots(uint32_t dim) const{
		assert(dim<ndim);
		return(nknots[dim]);
	}
	const double* get_knots(uint32_t dim) const{
		assert(dim<ndim);
		return(knots[dim]);
	}
	double get_knot(uint32_t dim, uint64_t knot) const{
		assert(dim<ndim);
		assert(knot<nknots[dim]);
		return(knots[dim][knot]);
	}
	double lower_extent(uint32_t dim) const{
		assert(dim<ndim);
		return(extents[dim][0]);
	}
	double upper_extent(uint32_t dim) const{
		assert(dim<ndim);
		return(extents[dim][1]);
	}
	double get_period(uint32_t dim) const{
		assert(dim<ndim);
		return(periods[dim]);
	}
	///Get the total number of spline coefficients
	uint64_t get_ncoeffs() const{
		return(std::accumulate(naxes,naxes+ndim,1ULL,std::multiplies<uint64_t>()));
	}
	///Get the number of coefficients along a given dimension
	uint64_t get_ncoeffs(uint32_t dim) const{
		assert(dim<ndim);
		return(naxes[dim]);
	}
	uint64_t get_stride(uint32_t dim) const{
		assert(dim<ndim);
		return(strides[dim]);
	}
	///Raw access to the coefficients. Use with care.
	float* get_coefficients(){
		return(&coefficients[0]);
	}
	const float* get_coefficients() const{
		return(&coefficients[0]);
	}
	
#ifdef PHOTOSPLINE_INCLUDES_SPGLAM
	
	///constant used to specify that no dimension should be forced to be
	///monotonic when fitting.
	static constexpr uint32_t no_monodim=PHOTOSPLINE_GLAM_NO_MONODIM;
	
	///Construct an approximating spline for possibly sparse, rectangularly gridded data
	///
	///\param data the datapoints to be fit, whose locations are specifed as
	///            indices into the coordinates arrays
	///\param weights the relative weights associated with the data points
	///\param coordinates a set of arrays containing the values of the abscissas
	///                   describing the rectangular grid on which the data
	///                   points are located
	///\param splineOrder an array of the b-spline orders which the fitted spline
	///                   should have in each dimension
	///\param knots a set of arrays containing the positions for the knots of the
	///             fitted spline in each dimension
	///\param smoothing an array specifying the strength of the regularization
	///                 to apply in each fit dimension
	///\param penaltyOrder an array specifying the order of the regularization
	///                    to apply in each fit dimension
	///\param monodim optionally, a single dimension in which the fit should be
	///               constrained to be monotonic
	///\param verbose whether to print logging and progress messages to standard
	///               output
	template<typename DoubleCont, typename UInt32Cont, typename DoubleContCont>
	void fit(const ::ndsparse& data,
	         const DoubleCont& weights,
	         const DoubleContCont& coordinates,
	         const UInt32Cont& splineOrder,
	         const DoubleContCont& knots,
	         const DoubleCont& smoothing,
	         const UInt32Cont& penaltyOrder,
	         uint32_t monodim=no_monodim, bool verbose=true);
	
	///Evaluate the spline on a grid
	///
	///\param coords a set of arrays specifying the evaluation points in each
	///            dimension
	template<typename DoubleContCont>
	std::unique_ptr<ndsparse> grideval(const DoubleContCont& coords) const;
	
#endif //PHOTOSPLINE_INCLUDES_SPGLAM
	
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
	template<size_t N, typename Distribution, typename RNG, typename Transform>
	std::vector<std::array<double,N>> sample(
	  size_t nresults, size_t burnin, std::array<size_t,N> samplingDimensions,
	  int derivatives, std::vector<double> coordinates,
	  Distribution distribution, RNG& rng,
	  Transform transform) const;

	template<size_t N, typename Distribution, typename RNG>
	std::vector<std::array<double,N>> sample(
	  size_t nresults, size_t burnin, std::array<size_t,N> samplingDimensions,
	  int derivatives, std::vector<double> coordinates,
	  Distribution distribution, RNG& rng) const;

	template<size_t N, typename Distribution, typename RNG>
	std::vector<std::array<double,N>> sample(
	  size_t nresults, size_t burnin, std::array<size_t,N> samplingDimensions,
	  std::vector<double> coordinates,
	  Distribution distribution, RNG& rng) const;

	///Reorder the dimensions of the spline
	///\param permutation the new order in which the current spline dimensions
	///       should appear
	void permuteDimensions(const std::vector<size_t>& permutation);
private:
	
	uint32_t ndim;
	uint32_t_ptr order;
	
	double_ptr_ptr knots;
	uint64_t_ptr nknots;
	
	double_ptr_ptr extents;
	
	double_ptr periods;
	
	float_ptr coefficients;
	uint64_t_ptr naxes;
	uint64_t_ptr strides;
	
	uint32_t naux;
	char_ptr_ptr_ptr aux;
	
	allocator_type allocator;
	uint32_t constOrder;
	double (splinetable::*eval_ptr)(const int*, int, detail::buffer2d<float>) const;
	void (splinetable::*v_eval_ptr)(const int*, const v4sf***, v4sf*) const;
	
	splinetable(const splinetable&);
	splinetable& operator=(const splinetable& other);
	
	/*
	 * The N-Dimensional tensor product basis version of splineeval.
	 * Evaluates the results of a full spline basis given a set of knots,
	 * a position, an order, and a central spline for the position (or -1).
	 * The central spline should be the index of the 0th order basis spline
	 * that is non-zero at the position x.
	 *
	 * x is the vector at which we will evaluate the space
	 */
	double ndsplineeval_core(const int* centers, int maxdegree, detail::buffer2d<float> localbasis) const;
	template<unsigned int D>
	double ndsplineeval_coreD(const int* centers, int maxdegree, detail::buffer2d<float> localbasis) const;
	template<unsigned int D, unsigned int O>
	double ndsplineeval_coreD_FixedOrder(const int* centers, int maxdegree, detail::buffer2d<float> localbasis) const;
	template<unsigned int ... Orders>
	double ndsplineeval_core_KnownOrder(const int* centers, int maxdegree, detail::buffer2d<float> localbasis) const;
	
	void ndsplineeval_multibasis_core(const int *centers, const v4sf*** localbasis, v4sf* result) const;
	template<unsigned int D>
	void ndsplineeval_multibasis_coreD(const int *centers, const v4sf*** localbasis, v4sf* result) const;
	template<unsigned int D, unsigned int O>
	void ndsplineeval_multibasis_coreD_FixedOrder(const int *centers, const v4sf*** localbasis, v4sf* result) const;
	template<unsigned int ... Orders>
	void ndsplineeval_multibasis_core_KnownOrder(const int *centers, const v4sf*** localbasis, v4sf* result) const;
	
	template<typename T>
	typename allocator_traits::template rebind_traits<T>::pointer allocate(size_t n){
		typedef typename allocator_traits::template rebind_alloc<T> other_alloc_t;
		typedef std::allocator_traits<other_alloc_t> other_alloc_traits;
		other_alloc_t other_alloc(allocator);
		return(other_alloc_traits::allocate(other_alloc,n));
	}
	
	template<typename Ptr, typename T=typename std::pointer_traits<Ptr>::element_type>
	void deallocate(Ptr buf, size_t n){
		typedef typename allocator_traits::template rebind_alloc<T> other_alloc_t;
		typedef std::allocator_traits<other_alloc_t> other_alloc_traits;
		other_alloc_t other_alloc(allocator);
		other_alloc_traits::deallocate(other_alloc,buf,n);
	}
	
	///Read from a file
	bool read_fits_core(fitsfile*, const std::string& filePath="");
	
	///Write to a file
	void write_fits_core(fitsfile*) const;
};
	
} //namespace photospline

#include "photospline/detail/bspline_eval.h"
#include "photospline/detail/bspline_multi.h"
#include "photospline/detail/aux.h"
#include "photospline/detail/convolve.h"
#include "photospline/detail/fitsio.h"
#include "photospline/detail/sample.h"
#include "photospline/detail/permute.h"

#ifdef PHOTOSPLINE_INCLUDES_SPGLAM
#include "photospline/detail/fit.h"
#include "photospline/detail/grideval.h"
#endif

#include "photospline/detail/config.h"

#endif /* PHOTOSPLINE_SPLINETABLE_H */
