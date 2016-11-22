#include "photospline/cinter/splinetable.h"
#include "photospline/splinetable.h"

#ifdef __cplusplus
extern "C" {
#endif
	
int splinetable_init(struct splinetable* table){
	if(!table)
		return(1);
	try{
		table->data=new photospline::splinetable<>();
	}catch(std::exception& ex){
		fprintf(stderr,"%s\n",ex.what());
		return(1);
	}catch(...){
		return(1);
	}
	return(0);
}
	
void splinetable_free(struct splinetable* table){
	if(!table)
		return;
	auto real_table=static_cast<photospline::splinetable<>*>(table->data);
	delete real_table;
	table->data=NULL;
}

int readsplinefitstable(const char* path, struct splinetable* table){
	if(!path || !table)
		return(1);
	if(table->data)
		splinetable_free(table);
	try{
		table->data=new photospline::splinetable<>(path);
	}catch(std::exception& ex){
		fprintf(stderr,"%s\n",ex.what());
		return(1);
	}catch(...){
		return(1);
	}
	return(0);
}

int writesplinefitstable(const char* path, const struct splinetable* table){
	if(!path || !table)
		return(1);
	try{
		const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
		real_table.write_fits(path);
	}catch(std::exception& ex){
		fprintf(stderr,"%s\n",ex.what());
		return(1);
	}catch(...){
		return(1);
	}
	return(0);
}

const char* splinetable_get_key(const struct splinetable* table, const char* key){
	if(!table || !table->data || !key)
		return(NULL);
	try{
		const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
		return(real_table.get_aux_value(key));
	}catch(std::exception& ex){
		fprintf(stderr,"%s\n",ex.what());
		return(NULL);
	}catch(...){
		return(NULL);
	}
	return(NULL);
}

int splinetable_read_key(const struct splinetable* table, splinetable_dtype type,
						 const char* key, void* result){
	if(!table || !table->data || !key || !result)
		return(1);
	try{
		const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
		switch(type){
			case SPLINETABLE_INT:
				real_table.read_key(key,*static_cast<int*>(result));
				break;
			case SPLINETABLE_DOUBLE:
				real_table.read_key(key,*static_cast<double*>(result));
				break;
		}
	}catch(std::exception& ex){
		fprintf(stderr,"%s\n",ex.what());
		return(1);
	}catch(...){
		return(1);
	}
	return(0);
}

int splinetable_write_key(struct splinetable* table, splinetable_dtype type,
                          const char* key, const void* value){
	if(!table || !table->data || !key || !value)
		return(1);
	try{
		auto& real_table=*static_cast<photospline::splinetable<>*>(table->data);
		switch(type){
			case SPLINETABLE_INT:
				real_table.write_key(key,*static_cast<const int*>(value));
				break;
			case SPLINETABLE_DOUBLE:
				real_table.write_key(key,*static_cast<const double*>(value));
				break;
		}
	}catch(std::exception& ex){
		fprintf(stderr,"%s\n",ex.what());
		return(1);
	}catch(...){
		return(1);
	}
	return(0);
}
	
uint32_t splinetable_ndim(const struct splinetable* table){
	const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
	return(real_table.get_ndim());
}
uint32_t splinetable_order(const struct splinetable* table, uint32_t dim){
	const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
	return(real_table.get_order(dim));
}
uint64_t splinetable_nknots(const struct splinetable* table, uint32_t dim){
	const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
	return(real_table.get_nknots(dim));
}
const double* splinetable_knots(const struct splinetable* table, uint32_t dim){
	const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
	return(real_table.get_knots(dim));
}
double splinetable_knot(const struct splinetable* table, uint32_t dim,
                        uint64_t knot){
	const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
	return(real_table.get_knot(dim,knot));
}
double splinetable_lower_extent(const struct splinetable* table, uint32_t dim){
	const auto& real_table=*static_cast<photospline::splinetable<>*>(table->data);
	return(real_table.lower_extent(dim));
}
double splinetable_upper_extent(const struct splinetable* table, uint32_t dim){
	const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
	return(real_table.upper_extent(dim));
}
double splinetable_period(const struct splinetable* table, uint32_t dim){
	const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
	return(real_table.get_period(dim));
}
uint64_t splinetable_ncoeffs(const struct splinetable* table, uint32_t dim){
	const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
	return(real_table.get_ncoeffs(dim));
}
uint64_t splinetable_total_ncoeffs(const struct splinetable* table){
	const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
	return(real_table.get_ncoeffs());
}
uint64_t splinetable_stride(const struct splinetable* table, uint32_t dim){
	const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
	return(real_table.get_stride(dim));
}
const float* splinetable_coefficients(const struct splinetable* table){
	const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
	return(real_table.get_coefficients());
}
	
int tablesearchcenters(const struct splinetable* table, const double* x,
                       int* centers){
	const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
	return(real_table.searchcenters(x,centers));
}
	
double ndsplineeval(const struct splinetable* table, const double* x,
                    const int* centers, int derivatives){
	const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
	return(real_table.ndsplineeval(x,centers,derivatives));
}
	
void ndsplineeval_gradient(const struct splinetable* table, const double* x,
                           const int* centers, double* evaluates){
	const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
	real_table.ndsplineeval_gradient(x,centers,evaluates);
}
	
double ndsplineeval_deriv2(const struct splinetable* table, const double* x,
                           const int* centers, int derivatives){
	const auto& real_table=*static_cast<const photospline::splinetable<>*>(table->data);
	return(real_table.ndsplineeval_deriv2(x,centers,derivatives));
}
	
int splinetable_convolve(struct splinetable* table, const int dim,
                         const double* knots, size_t n_knots){
	auto& real_table=*static_cast<photospline::splinetable<>*>(table->data);
	real_table.convolve(dim, knots, n_knots);
	return(0);
}
	
int readsplinefitstable_mem(const struct splinetable_buffer* buffer,
                            struct splinetable* table){
	if(!buffer || !buffer->data || !table)
		return(1);
	try{
		if(!table->data)
			table->data=new photospline::splinetable<>();
		auto& real_table=*static_cast<photospline::splinetable<>*>(table->data);
		real_table.read_fits_mem(buffer->data, buffer->size);
	}catch(std::exception& ex){
		fprintf(stderr,"%s\n",ex.what());
		return(1);
	}catch(...){
		return(1);
	}
	return(0);
}
	
int writesplinefitstable_mem(struct splinetable_buffer* buffer,
                             const struct splinetable* table){
	if(!buffer || buffer->data || !table)
		return(1);
	try{
		auto& real_table=*static_cast<photospline::splinetable<>*>(table->data);
		auto result=real_table.write_fits_mem();
		buffer->data=result.first;
		buffer->size=result.second;
	}catch(std::exception& ex){
		fprintf(stderr,"%s\n",ex.what());
		return(1);
	}catch(...){
		return(1);
	}
	return(0);
}

#ifdef PHOTOSPLINE_INCLUDES_SPGLAM
int splinetable_glamfit(struct splinetable* table, const struct ndsparse* data,
						const double* weights, const double* const* coords,
						const uint32_t* splineOrder, const double* const* knots,
						const uint64_t* nknots,
						const double* smoothing, const uint32_t* penaltyOrder,
						uint32_t monodim, bool verbose){
	if(!table || !table->data || !data)
		return(1);
	try{
		auto& real_table=*static_cast<photospline::splinetable<>*>(table->data);
		
		//TODO: this is rather suboptimal: the C++ interface wants vectors,
		//which want to own their data, so we have to copy everything.
		//Template the C++ fit function on container type, and provide a shim
		//type with an STL interface but which doesn't own it storage?
		std::vector<double> weightsv(data->rows);
		std::copy_n(weights,data->rows,weightsv.begin());
		
		std::vector<std::vector<double>> coordsv(data->ndim);
		for(size_t i=0; i<data->ndim; i++){
			coordsv[i].resize(data->ranges[i]);
			std::copy_n(coords[i],data->ranges[i],coordsv[i].begin());
		}
		
		std::vector<uint32_t> splineOrderv(data->ndim);
		std::copy_n(splineOrder,data->ndim,splineOrderv.begin());
		
		std::vector<std::vector<double>> knotsv(data->ndim);
		for(size_t i=0; i<data->ndim; i++){
			knotsv[i].resize(nknots[i]);
			std::copy_n(knots[i],nknots[i],knotsv[i].begin());
		}
		
		std::vector<double> smoothingv(data->ndim);
		std::copy_n(smoothing,data->ndim,smoothingv.begin());
		
		std::vector<uint32_t> penaltyOrderv;
		std::copy_n(penaltyOrder,data->ndim,penaltyOrderv.begin());
		
		real_table.fit(*data,weightsv,coordsv,splineOrderv,knotsv,smoothingv,
					   penaltyOrderv,monodim,verbose);
	}catch(std::exception& ex){
		fprintf(stderr,"%s\n",ex.what());
		return(1);
	}catch(...){
		return(1);
	}
	return(0);
}
#endif
	
int splinetable_permute(struct splinetable* table, size_t* permutation){
	try{
		auto& real_table=*static_cast<photospline::splinetable<>*>(table->data);
		std::vector<size_t> permutationv(real_table.get_ndim());
		//the user had better have supplied the right number of entries; we have no way to check
		std::copy(permutation,permutation+real_table.get_ndim(),permutationv.begin());
		real_table.permuteDimensions(permutationv);
	}catch(std::exception& ex){
		fprintf(stderr,"%s\n",ex.what());
		return(1);
	}catch(...){
		return(1);
	}
	return(0);
}
	
#ifdef __cplusplus
} //extern "C"
#endif