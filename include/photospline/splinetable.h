
#ifndef PHOTOSPLINE_SPLINETABLE_H
#define PHOTOSPLINE_SPLINETABLE_H

#include <algorithm>
#include <cassert>
#include <memory>
#include <numeric>
#include <sstream>
#include <vector>
#include <cstdlib>

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
	
	///Construct an empty splinetable.
	///The resulting object is useful only for calling read_fits, read_fits_mem, or fit.
	explicit splinetable(allocator_type alloc=Alloc()):
	ndim(0),order(NULL),knots(NULL),nknots(NULL),extents(NULL),periods(NULL),
	coefficients(NULL),naxes(NULL),strides(NULL),naux(0),aux(NULL),allocator(alloc)
	{}
	
	///Construct a splinetable from serialized data previously stored in a FITS file.
	///\param filePath the path to the input file
	explicit splinetable(const std::string& filePath, allocator_type alloc=Alloc()):
	ndim(0),order(NULL),knots(NULL),nknots(NULL),extents(NULL),periods(NULL),
	coefficients(NULL),naxes(NULL),strides(NULL),naux(0),aux(NULL),allocator(alloc)
	{
		read_fits(filePath);
	}

	///Construct a splinetable from stacking other splines.
	///\param splines the splines to stack
	///\param coordinates the coordiantes in the stacking dimension of the tables to be stacked
	///\param stackOrder the order of the spline in the stacking dimension
	explicit splinetable(std::vector<splinetable<Alloc>*> tables, std::vector<double> coordinates, int stackOrder=2, allocator_type alloc=Alloc()):
	ndim(0),order(NULL),knots(NULL),nknots(NULL),extents(NULL),periods(NULL),
	coefficients(NULL),naxes(NULL),strides(NULL),naux(0),aux(NULL),allocator(alloc)
	{
    assert(!tables.empty());
    assert(tables.size()==coordinates.size());
    int inputDim=tables.front()->get_ndim();
    for(auto table : tables){
      assert(table->get_ndim() == inputDim);
      assert(table->get_ncoeffs() && tables.front()->get_ncoeffs());
      for(unsigned int i=0; i<inputDim; i++){
        assert(table->get_order(i) && tables.front()->get_order(i));
      }
    }

    // add padding dimensions
    {
      auto extrapolateSpline=[](const splinetable<Alloc>* s1, const splinetable<Alloc>* s2)->splinetable<Alloc>*{
        splinetable<Alloc>* snew = new splinetable<Alloc>();

        snew->ndim = s2->ndim;

        snew->order = snew->allocate<uint32_t>(s2->ndim);
        std::copy_n(s2->order,s2->ndim,snew->order);

        snew->nknots = snew->allocate<uint64_t>(s2->ndim);
        std::copy_n(s2->nknots,s2->ndim,snew->nknots);

        snew->knots = snew->allocate<double_ptr>(s2->ndim);
        for(unsigned int i=0; i<s2->ndim; i++){
          snew->knots[i] = snew->allocate<double>(s2->nknots[i]+2*s2->order[i]) + s2->order[i];
          std::copy_n(s2->knots[i],s2->nknots[i],snew->knots[i]);
        }

        snew->naxes = snew->allocate<uint64_t>(s2->ndim);
        std::copy_n(s2->naxes,s2->ndim,snew->naxes);

        snew->strides = snew->allocate<uint64_t>(s2->ndim);
        std::copy_n(s2->strides,s2->ndim,snew->strides);

        snew->extents = snew->allocate<double_ptr>(s2->ndim);
        snew->extents[0] = snew->allocate<double>(2*s2->ndim);
        for(unsigned int i=0;i<s2->ndim; i++){
          snew->extents[i] = &snew->extents[0][2*i];
        }

        for(unsigned int i=0; i<s2->ndim; i++){
          snew->extents[i][0] = s2->extents[i][0];
          snew->extents[i][1] = s2->extents[i][1];
        }

        snew->periods = NULL;
        snew->naux = 0;
        snew->aux = NULL;

        unsigned long nCoeffs=snew->get_ncoeffs();
        snew->coefficients=snew->allocate<float>(nCoeffs);
        for(unsigned long i=0; i<nCoeffs; i++){
          auto c1=s1->get_coefficients()[i];
          auto c2=s2->get_coefficients()[i];
          snew->get_coefficients()[i]=2*c2-c1;
        }

        return(snew);
      };

      tables.insert(tables.begin(),extrapolateSpline(tables[1],tables[0]));
      coordinates.insert(coordinates.begin(),2*coordinates[0]-coordinates[1]);

      tables.push_back(extrapolateSpline(tables[tables.size()-2],tables[tables.size()-1]));
      coordinates.push_back(2*coordinates[coordinates.size()-1]-coordinates[coordinates.size()-2]);
    }

    //set dimensions
    ndim=inputDim+1;
    //copy/set spline orders and knots
    order=allocate<uint32_t>(ndim);
    nknots=allocate<uint64_t>(ndim);

    for(unsigned int i=0; i<inputDim; i++){
      order[i] = tables.front()->get_order(i);
      nknots[i] = tables.front()->get_nknots(i);
    }
    order[inputDim]=stackOrder;
    nknots[inputDim]=tables.size()+stackOrder+1;

    knots=allocate<double_ptr>(ndim);
    //copy existing knots
    for(unsigned int i=0; i<inputDim; i++){
      knots[i]=allocate<double>(nknots[i]+2*order[i]) + order[i];
      std::copy_n(tables.front()->get_knots(i),nknots[i],knots[i]);
    }
    //figure out knots for the new dimension
    {
      knots[inputDim]=allocate<double>(nknots[inputDim]+2*order[inputDim]) + order[inputDim];
      double_ptr lastKnots=knots[inputDim];
      //copy input positions
      std::copy_n(coordinates.begin(),tables.size(),lastKnots+stackOrder);

      //shift knots
      double knotShift=(stackOrder-1)*(lastKnots[stackOrder+tables.size()-1]-lastKnots[stackOrder])/(2*tables.size());
      for(unsigned int i=0; i<tables.size(); i++)
        lastKnots[stackOrder+i]+=knotShift;

      //add stackOrder padding knots before
      double knotStep=lastKnots[stackOrder+1]-lastKnots[stackOrder];
      for(int i=0; i<stackOrder; i++)
        lastKnots[i]=lastKnots[stackOrder]+(i-stackOrder)*knotStep;
      //add one padding knot after
      lastKnots[nknots[inputDim]-1]=2*lastKnots[nknots[inputDim]-2]-lastKnots[nknots[inputDim]-3];
    }

    //set naxes
    naxes=allocate<uint64_t>(ndim);
    for(unsigned int i=0; i<inputDim; i++)
      naxes[i] = tables.front()->get_ncoeffs(i);
    naxes[inputDim]=tables.size();

    //copy coefficients
    unsigned long nCoeffs=std::accumulate(naxes, naxes+ndim, 1UL, std::multiplies<uint64_t>());
    unsigned long nInputCoeffs=std::accumulate(naxes, naxes+ndim-1, 1UL, std::multiplies<uint64_t>());
    coefficients=allocate<float>(nCoeffs);
    unsigned int step=naxes[ndim-1];
    for(unsigned int i=0; i<tables.size(); i++){
      for(unsigned int j=0; j<nInputCoeffs; j++)
        coefficients[i+j*step]=tables[i]->get_coefficients()[j];
    }

    //set strides
    strides = allocate<uint64_t>(ndim);
    uint64_t arraysize;
    strides[ndim-1] = arraysize = 1;
    for(int i=ndim-1; i >= 0; i--){
      arraysize *= naxes[i];
      if(i>0)
        strides[i-1] = arraysize;
    }
	}

	splinetable(splinetable&& other):
	ndim(other.ndim),order(std::move(other.order)),knots(std::move(other.knots)),
	nknots(other.nknots),extents(std::move(other.extents)),
	periods(std::move(other.periods)),coefficients(std::move(other.coefficients)),
	naxes(std::move(other.naxes)),strides(std::move(other.strides)),
	naux(other.naux),aux(std::move(other.aux)),
	allocator(std::move(other.allocator))
	{
		other.ndim=0;
		other.order=NULL;
		other.knots=NULL;
		other.nknots=NULL;
		other.extents=NULL;
		other.periods=NULL;
		other.coefficients=NULL;
		other.naxes=NULL;
		other.strides=NULL;
		other.naux=0;
		other.aux=NULL;
		other.allocator=Alloc();
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
	
	splinetable& operator=(splinetable&& other){
		if(&other==this)
			return(*this);
		using std::swap;
		swap(ndim,other.ndim);
		swap(order,other.order);
		swap(knots,other.knots);
		swap(nknots,other.nknots);
		swap(extents,other.extents);
		swap(periods,other.periods);
		swap(coefficients,other.coefficients);
		swap(naxes,other.naxes);
		swap(strides,other.strides);
		swap(naux,other.naux);
		swap(aux,other.aux);
		swap(allocator,other.allocator);
		return(*this);
	}
	
	///Compare splines for equality.
	///Splines are considered equal if evaluation would return the same result
	///at all possible coordinate values. Auxiliary FITS header entries are not
	///considered.
	bool operator==(const splinetable& other) const{
		if (ndim != other.ndim)
			return false;
		if (!std::equal(order,order+ndim,other.order))
			return false;
		if (!std::equal(naxes,naxes+ndim,other.naxes))
			return false;
		if (!std::equal(nknots,nknots+ndim,other.nknots))
			return false;
		for (uint32_t i=0; i<ndim; i++)
			if (!std::equal(knots[i],knots[i]+nknots[i],other.knots[i]))
				return false;
		if (get_ncoeffs() != other.get_ncoeffs())
			return false;
		if (!std::equal(coefficients,coefficients+get_ncoeffs(),other.coefficients))
			return false;
		return true;
	}
	
	///Compare splines for inequality.
	///Splines are considered equal if evaluation would return the same result
	///at all possible coordinate values. Auxiliary FITS header entries are not
	///considered.
	bool operator!=(const splinetable& other) const{
		//there doesn't appear to be any particular performance disadvantage to 
		//doing this, we mostly just supply this operator for greater user 
		//convenience.
		return(!operator==(other));
	}

	///Estimate the memory needed to load an existing spline.
	///This function is intended for users using specialized allocators for which
	///it is useful to know the total memory required before constructing the
	///allocator and splinetable. If a convolution will be applied to the spline
	///after it is loaded, doing so will increase its memory requirements, which
	///this function can also take into account.
	///\param filePath the path to the input FITS file
	///\param n_convolution_knots the number of knots possessed by the other
	///       spline with which one dimension of this spline will be convolved.
	///       A value of one is equivalent to no convolution.
	///\param convolution_dimension the dimension of this spline which will have
	///       a convolution applied.
	static size_t estimateMemory(const std::string& filePath,
	                             uint32_t n_convolution_knots = 1,
	                             uint32_t convolution_dimension = 0);
	
	///Read from a FITS file
	///\param path the path to the input file
	bool read_fits(const std::string& path);
	
	///Read from a FITS 'file' in a memory buffer
	///\param the input data buffer
	///\param buffer_size the length of the input buffer
	bool read_fits_mem(void* buffer, size_t buffer_size);
	
	///Write to a FITS file
	///\param path the path to the output file
	void write_fits(const std::string& path) const;
	
	///Write to a FITS memory 'file'.
	///\returns A pair containing a buffer allocated by malloc(), which should be
	///         deallocated with free(), and the size of that buffer.
	std::pair<void*,size_t> write_fits_mem() const;
	
	///Get the count of 'auxiliary' keys
	size_t get_naux_values() const{ return(naux); }
	///Directly get a particular auxiliary key
	const char* get_aux_key(size_t i) const{ return(aux[i][0]); }
	///Directly get a particular auxiliary value
	///\param key the key whose value should be fetched
	///\return the value if key exists, otherwise NULL
	const char* get_aux_value(const char* key) const;
	///Delete an auxiliary key and associated value
	///\returns true if the key existed and was removed
	bool remove_key(const char* key);
	///Look up the value associated with an auxiliary key
	///\param key the name of the key to look up
	///\param result location to store the value if the key is found
	///\return true if the key was found and the corresponding value was
	///        sucessfully parsed into result, otherwise false
	template<typename T>
	bool read_key(const char* key, T& result) const;
	///Look up the value associated with an auxiliary key
	///\note The data extracted may be padded with extra whitespace if it has
	///      been read from a FITS file. 
	///\param key the name of the key to look up
	///\param result location to store the value if the key is found
	///\return true if the key was found
	bool read_key(const char* key, std::string& result) const;
	///Insert or overwrite an auxiliary key,value pair
	///\param key the name of the key to store
	///\param value the value to store for the key
	///\return whether the value was successfully stored
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
	template <typename Float=float>
	struct evaluator_type{
	private:
		const splinetable<Alloc>& table;
		double (splinetable::*eval_ptr)(const int*, int, detail::buffer2d<Float>) const;
		void (splinetable::*v_eval_ptr)(const int*, const typename detail::simd_vector<Float>::type***, typename detail::simd_vector<Float>::type*) const;
		friend class splinetable<Alloc>;
		evaluator_type(const splinetable<Alloc>& table):table(table){}
	public:
		///\brief Get the underlying splinetable
		const splinetable<Alloc>& get_table() const{ return(table); }
		///\brief same as splinetable::searchcenters
		bool searchcenters(const double* x, int* centers) const;
		///\brief same as splinetable::ndsplineeval
		double ndsplineeval(const double* x, const int* centers, int derivatives=0) const;
		///\brief Convenince short-cut for ndsplineeval
		double operator()(const double* x, int derivatives=0) const;
		///\brief same as splinetable::ndsplineeval_gradient
		void ndsplineeval_gradient(const double* x, const int* centers, double* evaluates) const;
		///\brief same as splinetable::ndsplineeval_deriv
		double ndsplineeval_deriv(const double* x, const int* centers, const unsigned int *derivatives) const;
	};
	template <typename Float> friend struct evaluator_type;
	
	///Constructs an optimized evaluator object which will use the best
	///available internal routines to perform evaulations, at the requested
	///precision. The evaluator holds a reference to this splinetable, so it
	///must be considered invalidated if this table altered or destroyed.
	template <typename Float=float>
	evaluator_type<Float> get_evaluator() const;

	typedef evaluator_type<float> evaluator;
	
	/*
	 * Spline table based hypersurface evaluation. ndsplineeval() takes a spline
	 * coefficient table, a vector at which to evaluate the surface, and a vector
	 * indicating the evaluation centers, as for splineeval().
	 *
	 * searchcenters() provides a method to acquire a centers vector
	 * for ndsplineeval() using a binary search. Depending on how the table
	 * was produced, a more efficient method may be available.
	 */
	///Acquire a centers vector for use with the ndsplineeval functions.
	///This performs a binary search over the knots in each spline dimension, so
	///depending on how the spline was produced a more efficient method may be exist.
	///\param x a vector of coordinates at which the spline is to be evaluated
	///\param centers a vector of indices which will be populated by this function
	///\return whether centers was sucessfully populated
	///\pre x and centers must both have lengths matching the spline's dimension
	bool searchcenters(const double* x, int* centers) const;
	
	///Evaluate the spline hypersurface.
	///\param x a vector of coordinates at which the spline is to be evaluated
	///\param centers a vector of knot indices derived from x, constructed using
	///       searchcenters
	///\param derivatives a bitmask indicating in which dimensions the spline
	///       should be differentiated rather than directly evaluated. Mixed
	///       partial derivatives are supported.
	///\return the spline value or derivative value
	template <typename Float=float>
	double ndsplineeval(const double* x, const int* centers, int derivatives) const;
	
	///Evaluate the spline hypersurface.
	///This convenience interface for evaluation simply performs searchcenters
	///and ndsplineeval. It assumes that no derivative is desired, and yields
	///zero when searchcenters fails. This is suitable for simple uses where the
	///user knows that the coordinates will be inside the table bounds.
	///\param x a vector of coordinates at which the spline is to be evaluated
	///\return the spline value or zero if center lookup fails
	double operator()(const double* x) const;
	
	///Evaluate arbitrary derivative of the spline hypersurface.
	///\param x a vector of coordinates at which the spline is to be evaluated
	///\param centers a vector of knot indices derived from x, constructed using
	///       searchcenters
	///\param derivatives a vector giving the order of derivative to compute
	///       for each dimension
	///\return the derivative of the spline
	double ndsplineeval_deriv(const double* x, const int* centers, const unsigned int *derivatives) const;
	
	///Evaluate the spline hypersurface and its gradient.
	///If the full gradient is needed along with the spline value this function
	///can be more efficient than repeated calls to ndsplineeval.
	///\param x a vector of coordinates at which the spline is to be evaluated
	///\param centers a vector of knot indices derived from x, constructed using
	///       searchcenters
	///\param evaluates a vector which will be populated by this function with
	///       the spline value, followed by the components of the gradient.
	///\pre evaluates must have a length one greater than the dimension of the spline
	template <typename Float=float>
	void ndsplineeval_gradient(const double* x, const int* centers, double* evaluates) const;
	
	///A container for results obtained from benchmark_evaluation
	struct benchmark_results{
		///The rate at which ndsplineeval can evaluate the value of the spline
		double single_eval_rate;
		///The rate at which ndsplineeval can evaluate the gradient of the spline
		///(requiring a number of calls to ndsplineeval equal to the spline
		///dimension per evaluation)
		double gradient_single_eval_rate;
		///The rate at which ndsplineeval_gradient can evaluate the value and
		///gradient of the spline
		double gradient_multi_eval_rate;
	};
	///Evaluate the spline at random points within the extent
	///\param trialCountthe number of times each type of evaluation should be
	///       performed. This should be large enough to get a stable result.
	///\param verbose whether the results should be printed to stdout.
	template <typename Float=float>
	benchmark_results benchmark_evaluation(size_t trialCount=1e4, bool verbose=false);
	
	///Convolve a single dimension of this spline with another spline.
	///This rewrites the data of this spline in place. The order of the spline
	///will be raised by (n_knots - 1) in the given dimension, i.e. convolving
	///with an order-0 spline (a box function, defined on two knots) will raise
	///the order of the spline surface by 1.
	///\param dim the dimension in which the convolution should be done
	///\param knots a vector of positions of the knots of the one dimensional
	///       spline used for the convolution
	///\param n_knots the length of knots
	void convolve(const uint32_t dim, const double* knots, size_t n_knots);
	
	///Get the dimension of the spline
	uint32_t get_ndim() const{ return(ndim); }
	///Get the order of the spline in a given dimension
	uint32_t get_order(uint32_t dim) const{
		assert(dim<ndim);
		return(order[dim]);
	}
	///Get the number of knots in a given dimension
	uint64_t get_nknots(uint32_t dim) const{
		assert(dim<ndim);
		return(nknots[dim]);
	}
	///Get the knot vector for a given dimension
	const double* get_knots(uint32_t dim) const{
		assert(dim<ndim);
		return(knots[dim]);
	}
	///Get a particular knot from a given dimension
	///\param dim the dimension
	///\param knot the index of the knot to fetch
	double get_knot(uint32_t dim, uint64_t knot) const{
		assert(dim<ndim);
		assert(knot<nknots[dim]);
		return(knots[dim][knot]);
	}
	///Get the left boundary of the spline in a given dimension
	double lower_extent(uint32_t dim) const{
		assert(dim<ndim);
		return(extents[dim][0]);
	}
	///Get the right boundary of the spline in a given dimension
	double upper_extent(uint32_t dim) const{
		assert(dim<ndim);
		return(extents[dim][1]);
	}
	///Get the period of the spline in a given dimension
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
	///Get the stride through the coefficient array corresponding to a given dimension
	uint64_t get_stride(uint32_t dim) const{
		assert(dim<ndim);
		return(strides[dim]);
	}
	///Raw access to the coefficients. Use with care.
	float* get_coefficients(){
		return(&coefficients[0]);
	}
	///Raw access to the coefficients.
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
	
	///Draw a number of random samples from the distribution given by a slice
	///through the spline.
	///\tparam N the number of dimensions in which to sample
	///\tparam Distribution a type which can both evaluate a proposal pdf
	///        via double operator()(const std::vector<double>& coordinates) const
	///        and can draw samples from that same pdf
	///        via std::vector<double> sample(RNG rng)
	///\tparam RNG a random number generator, conforming to the standard interface
	///\tparam Transform a callable type which can transform values of the spline,
	///        given the vector of coordinates at which the spline was evaluated
	///        and the resulting value of the spline surface.
	///\param nresults the number of samples to draw
	///\param burnin the length of the burn in period to use
	///\param samplingDimensions the dimensions in which sampling is to be done
	///\param coordinates vector of coordinate values at which the spline is to
	///       be evaluated in the dimensions not being sampled. The length of this
	///       vector must be the same as spline dimension, but entries for the
	///       sampling dimensions will be ignored.
	///\param distribution the proposal distribution
	///\param rng the source of randomness used for sampling
	///\param derivatives a bitmask indicating in which dimensions the derivative
	///       of the spline surface is to be taken rather
	///\param transform the transformation function to apply to the spline suface
	///       before sampling. This is useful when one wishes to sample from a
	///       distribution, but has created a spline representing a transformed
	///       version of the pdf (for example its logarithm), in which case
	///       transform should be a function which undoes the transformation
	///       which was applied to the create the spline. Alternatively, one may
	///       wish to sample from a distribution in the space of some set of
	///       variables xi, but have created a spline in the space of some
	///       transformed variables xi', in which case transform should account
	///       for the jacobian of the change of variables. A combination of
	///       these cases is also possible.
	///\pre The entries of coordinates for dimensions which do not appear in 
	///        samplingDimensions must be with in the extents of spline
	template<size_t N, typename Distribution, typename RNG, typename Transform>
	std::vector<std::array<double,N>> sample(
	  size_t nresults, size_t burnin, std::array<size_t,N> samplingDimensions,
	  std::vector<double> coordinates, Distribution distribution, RNG& rng,
	  int derivatives, Transform transform) const;

	///A version of sample which does not apply a transformation.
	///See sample for more details.
	template<size_t N, typename Distribution, typename RNG>
	std::vector<std::array<double,N>> sample(
	  size_t nresults, size_t burnin, std::array<size_t,N> samplingDimensions,
	  std::vector<double> coordinates, Distribution distribution, RNG& rng,
	  int derivatives=0) const;

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
	template <typename Float>
	double ndsplineeval_core(const int* centers, int maxdegree, detail::buffer2d<Float> localbasis) const;
	template<typename Float, unsigned int D>
	double ndsplineeval_coreD(const int* centers, int maxdegree, detail::buffer2d<Float> localbasis) const;
	template<typename Float, unsigned int D, unsigned int O>
	double ndsplineeval_coreD_FixedOrder(const int* centers, int maxdegree, detail::buffer2d<Float> localbasis) const;
	template<typename Float, unsigned int ... Orders>
	double ndsplineeval_core_KnownOrder(const int* centers, int maxdegree, detail::buffer2d<Float> localbasis) const;
	
	template <typename Float>
	void ndsplineeval_multibasis_core(const int *centers, const typename detail::simd_vector<Float>::type*** localbasis, typename detail::simd_vector<Float>::type* result) const;
	template<typename Float, unsigned int D>
	void ndsplineeval_multibasis_coreD(const int *centers, const typename detail::simd_vector<Float>::type*** localbasis, typename detail::simd_vector<Float>::type* result) const;
	template<typename Float, unsigned int D, unsigned int O>
	void ndsplineeval_multibasis_coreD_FixedOrder(const int *centers, const typename detail::simd_vector<Float>::type*** localbasis, typename detail::simd_vector<Float>::type* result) const;
	template<typename Float, unsigned int ... Orders>
	void ndsplineeval_multibasis_core_KnownOrder(const int *centers, const typename detail::simd_vector<Float>::type*** localbasis, typename detail::simd_vector<Float>::type* result) const;
	
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
#include "photospline/detail/auxiliary.h"
#include "photospline/detail/convolve.h"
#include "photospline/detail/fitsio.h"
#include "photospline/detail/sample.h"
#include "photospline/detail/permute.h"

#ifdef PHOTOSPLINE_INCLUDES_SPGLAM
#include "photospline/detail/fit.h"
#include "photospline/detail/grideval.h"
#endif

#endif /* PHOTOSPLINE_SPLINETABLE_H */
