#ifndef PHOTOSPLINE_FITSIO_H
#define PHOTOSPLINE_FITSIO_H

#include <string.h>

namespace photospline{
	
std::vector<uint32_t> readOrder(fitsfile* fits, uint32_t ndim);
bool reservedFitsKeyword(const char* key);
uint32_t countAuxKeywords(fitsfile* fits);

template<typename Alloc>
size_t splinetable<Alloc>::estimateMemory(const std::string& filePath,
                                          uint32_t n_convolution_knots,
                                          uint32_t convolution_dimension){
	fitsfile* fits;
	int error = 0;

	fits_open_file(&fits, filePath.c_str(), READONLY, &error);
	
	struct fits_cleanup{
		fitsfile* fits;
		fits_cleanup(fitsfile* fits):fits(fits){}
		~fits_cleanup(){
			int error=0;
			fits_close_file(fits, &error);
			fits_report_error(stderr, error);
		}
	} cleanup(fits);
	
	if (error != 0)
		throw std::runtime_error("Unable to open "+filePath);
	
	{
		int hdus, type;
		fits_get_num_hdus(fits, &hdus, &error);
		fits_movabs_hdu(fits, 1, &type, &error);
		if (error != 0)
			throw std::runtime_error("Failed to move to HDU 1 in "+filePath);
		
		if (type != IMAGE_HDU)
			throw std::runtime_error("First HDU in "+filePath+" is not an image");
	}
	
	int dim;
	fits_get_img_dim(fits, &dim, &error);
	if (error != 0)
		throw std::runtime_error("Unable to read table dimension from "+filePath);
	
	//Find out how many coefficients there are
	std::vector<long> naxes(dim);
	fits_get_img_size(fits, dim, naxes.data(), &error);
	if (error != 0)
		throw std::runtime_error("Unable to read coefficient array 'image' size: Error "+std::to_string(error));
	std::reverse(naxes.begin(),naxes.end()); //pull necessary switcheroo
	
	std::vector<uint32_t> order = readOrder(fits,dim);
	order[convolution_dimension] += n_convolution_knots-1;
	
	size_t size = sizeof(splinetable<Alloc>); //main object
	
	//count knots
	for (int i = 0; i < dim; i++) {
		std::ostringstream hduname;
		hduname << "KNOTS" << i;
		fits_movnam_hdu(fits, IMAGE_HDU, const_cast<char*>(hduname.str().c_str()), 0, &error);
		long nknots;
		fits_get_img_size(fits, 1, &nknots, &error);
		
		if (error != 0) {
			throw std::runtime_error("Error reading knot vector "+std::to_string(i));
		}
		
		if (unsigned(i) == convolution_dimension){
			nknots *= n_convolution_knots;
			naxes[i] = nknots - order[i] - 1;
		}
		size += (nknots+2*order[i])*sizeof(double);
	}
		
	int64_t ncoeffs = std::accumulate(naxes.begin(),naxes.end(),(int64_t)1,std::multiplies<int64_t>());
	
	//count up size
	size += dim*sizeof(uint32_t); //order
	size += dim*sizeof(double_ptr); //knot pointers, knots themselves accounted for above
	size += dim*sizeof(uint64_t); //nknots
	size += 2*dim*sizeof(double)+dim*sizeof(double_ptr); //extents
	size += dim*sizeof(double); //periods
	size += ncoeffs*sizeof(float); //coefficients
	size += dim*sizeof(uint64_t); //naxes
	size += dim*sizeof(uint64_t); //strides
	size += dim*sizeof(double); //rmin_sep
	size += dim*sizeof(double); //rmax_sep
	
	uint32_t naux = countAuxKeywords(fits);
	//pessimistically assume all keys and values are maximal length
	size += naux*(FLEN_KEYWORD+FLEN_VALUE)*sizeof(char);
	
	const size_t KB=1ULL<<10;
	//round up to the nearest KB, and add one more,
	//to allow for a little overhead
	size += (KB-size%KB)+KB;
	
	return(size);
}

template<typename Alloc>
bool splinetable<Alloc>::read_fits(const std::string& filePath){
	if(ndim!=0)
		throw std::runtime_error("splinetable already contains data, cannot read from file");
	
	fitsfile* fits;
	int error = 0;

	fits_open_diskfile(&fits, filePath.c_str(), READONLY, &error);
	if (error != 0)
		throw std::runtime_error(("CFITSIO failed to open "+filePath+" for reading").c_str());
	
	struct fits_cleanup{
		fitsfile* fits;
		fits_cleanup(fitsfile* fits):fits(fits){}
		~fits_cleanup(){
			int error=0;
			fits_close_file(fits, &error);
			fits_report_error(stderr, error);
		}
	} cleanup(fits);
	return(read_fits_core(fits, filePath));
}
	
template<typename Alloc>
bool splinetable<Alloc>::read_fits_mem(void* buffer, size_t buffer_size){
	if(ndim!=0)
		throw std::runtime_error("splinetable already contains data, cannot read from (memory) file");
	
	fitsfile* fits;
	int error = 0;
	
	fits_open_memfile(&fits, "", READONLY, &buffer, &buffer_size, 0, NULL, &error);
	if (error != 0){
		fits_report_error(stderr, error);
		throw std::runtime_error("CFITSIO failed to open memory 'file' for reading");
	}
	
	struct fits_cleanup{
		fitsfile* fits;
		fits_cleanup(fitsfile* fits):fits(fits){}
		~fits_cleanup(){
			int error=0;
			fits_close_file(fits, &error);
			fits_report_error(stderr, error);
		}
	} cleanup(fits);
	return(read_fits_core(fits, "memory 'file'"));
}
	
template<typename Alloc>
bool splinetable<Alloc>::read_fits_core(fitsfile* fits, const std::string& filePath){
	int error = 0;
	//if (error != 0)
	//	throw std::runtime_error("Failed to move to HDU 1 in "+filePath);
	
	//Set the HDU and check its type
	{
		int hdus, type;
		fits_get_num_hdus(fits, &hdus, &error);
		fits_movabs_hdu(fits, 1, &type, &error);
		if (error != 0)
			throw std::runtime_error("Unable to move to first HDU: Error "+std::to_string(error));
		if (type != IMAGE_HDU)
			throw std::runtime_error("First HDU in "+filePath+" is not an image");
	}
	
	
	//Read header information
	{
		int temp_dim;
		fits_get_img_dim(fits, &temp_dim, &error);
		if (error != 0)
			throw std::runtime_error("Unable to read table dimension from "+filePath);
		if (temp_dim < 1)
			throw std::runtime_error("Invalid table dimension "+std::to_string(temp_dim));
		ndim = temp_dim;
	}
	
	//Read in any auxiliary keywords.
	{
		int nkeys = 0;
		fits_get_hdrspace(fits, &nkeys, NULL, &error);
		if (nkeys > 0) {
			char key[FLEN_KEYWORD], value[FLEN_VALUE];
			int keylen, valuelen;
			
			// figure out how many generic keys there are
			naux = 0;
			for (int j = 1 ; j-1 < nkeys; j++) {
				error = 0;
				fits_read_keyn(fits, j, key, value, NULL, &error);
				if (error != 0)
					continue;
				if (reservedFitsKeyword(key))
					continue;
				naux++;
			}
			
			aux = allocate<char_ptr_ptr>(naux);
			std::fill(aux,aux+naux,nullptr);
			
			for (unsigned i = 0, j = 1 ; (i < naux) && (j-1 < unsigned(nkeys)); j++) {
				error = 0;
				fits_read_keyn(fits, j, key, value, NULL, &error);
				if (error != 0)
					continue;
				if (reservedFitsKeyword(key))
					continue;
				
				keylen = strlen(key) + 1;
				valuelen = strlen(value) + 1;
				aux[i] = allocate<char_ptr>(2);
				aux[i][0] = aux[i][1] = NULL;
				aux[i][0] = allocate<char>(keylen);
				aux[i][1] = allocate<char>(valuelen);
				std::copy(key,key+keylen,aux[i][0]);
				//remove stupid quotes mandated by FITS, but not removed by cfitsio on reading
				//Note that we do not attempt to remove whitespace, because we cannot 
				//distinguish whitespace included by the user and whitespace pointlessly
				//added by FITS.
				if(valuelen>1 && value[0]=='\''){
					if(valuelen>2 && value[valuelen-2]=='\''){ //remove a trailing quote also
						std::copy(value+1,value+valuelen-2,aux[i][1]);
						aux[i][1][valuelen-3]='\0';
					}
					else{ //just remove an opening quote
						std::copy(value+1,value+valuelen-1,aux[i][1]);
						aux[i][1][valuelen-2]='\0';
					}
				}
				else{
					std::copy(value,value+valuelen,aux[i][1]);
					aux[i][1][valuelen-1]='\0';
				}
				i++;
			}
		} else {
			aux = NULL;
			naux = 0;
		}
	}
	
	//Read the spline orders
	order = allocate<uint32_t>(ndim);
	//See if there is a single order value
	fits_read_key(fits, TINT, "ORDER", &order[0], NULL, &error);
	if (error != 0) {
		error = 0;
		
		//There is not, so look for a separate order in each dimension
		for (unsigned i = 0; i < ndim; i++) {
			std::ostringstream ss;
			ss << "ORDER" << i;
			fits_read_key(fits, TUINT, ss.str().c_str(), &order[i], NULL, &error);
			if (error != 0)
				throw std::runtime_error("Unable to read order for dimension "+std::to_string(i));
		}
	} else {
		//all orders are the same
		std::fill(order+1,order+ndim,order[0]);
	}
	
	if (error != 0)
		return (error);
	
	//read the table periods
	periods = allocate<double>(ndim);
	for (unsigned i = 0; i < ndim; i++) {
		std::ostringstream ss;
		ss << "PERIOD" << i;
		fits_read_key(fits, TDOUBLE, ss.str().c_str(), &periods[i], NULL, &error);
		//If the PERIOD keys cannot be read, just interpret this as
		//a non-periodic table.
		if (error != 0) {
			periods[i] = 0;
			error = 0;
		}
	}
	
	//We won't read these things until later, but it's useful to allocate all
	//arrays which don't depend on the orders or numbers of knots before the
	//ones which do
	rmin_sep = allocate<double>(ndim);
	rmax_sep = allocate<double>(ndim);
	knots = allocate<double_ptr>(ndim);
	nknots = allocate<uint64_t>(ndim);
	extents = allocate<double_ptr>(ndim);
	extents[0] = allocate<double>(2*ndim);
	
	//Read the coefficient table
	std::vector<long> naxes_temp(ndim);
	fits_get_img_size(fits, ndim, naxes_temp.data(), &error);
	if (error != 0)
		throw std::runtime_error("Unable to read coefficient array 'image' size: Error "+std::to_string(error));
	for(size_t i=0; i<ndim; i++){
		if(naxes_temp[i]<0)
			throw std::runtime_error("Invalid size in dimension "+std::to_string(i));
	}
	naxes = allocate<uint64_t>(ndim);
	
	
	//FITS multidimensional arrays are stored as FORTRAN arrays,
	//not C arrays, so we need to swizzle the matrix into being
	//a C array. Or we should. Instead, PyFITS, which writes these
	//files, writes a C array, but with the axis counts transposed.
	//Fix it while copying to the actual array.
	std::copy(naxes_temp.rbegin(),naxes_temp.rend(),naxes);
	
	// Compute the total array size and the strides into each dimension
	strides = allocate<uint64_t>(ndim);
	strides[0]=1;
	std::partial_sum(naxes_temp.begin(),naxes_temp.end()-1,strides+1,std::multiplies<uint64_t>());
	std::reverse(strides,strides+ndim);
	uint64_t ncoeffs=strides[0]*naxes[0];
	coefficients = allocate<float>(ncoeffs);
	
	std::vector<long> fpixel(ndim,1);
	fits_read_pix(fits, TFLOAT, fpixel.data(), ncoeffs, NULL,
				  &coefficients[0], NULL, &error);
	
	if (error != 0){
		//destroy
		throw std::runtime_error("Error reading table coefficients");
	}
	
	//Read the knot vectors, which are stored one each in extension HDUs
	for (unsigned i = 0; i < ndim; i++) {
		std::ostringstream hduname;
		hduname << "KNOTS" << i;
		fits_movnam_hdu(fits, IMAGE_HDU, const_cast<char*>(hduname.str().c_str()), 0, &error);
		long nknots_temp;
		fits_get_img_size(fits, 1, &nknots_temp, &error);
		
		if (error != 0)
			throw std::runtime_error("Error reading size of knot vector "+std::to_string(i));
		if(nknots_temp<=0)
			throw std::runtime_error("Invalid number of knots ("+std::to_string(nknots_temp)+") in dimension "+std::to_string(i));
		nknots[i]=nknots_temp;
		
		//Allow spline evaluations to run off the ends of the
		//knot field without segfaulting.
		knots[i] = allocate<double>(nknots[i]+2*order[i]) + order[i];
		//TODO: should the 'off the end' entries of knots be set to zero?
		
		long fpix = 1;
		fits_read_pix(fits, TDOUBLE, &fpix, nknots[i], NULL, &knots[i][0], NULL, &error);
		if (error != 0)
			throw std::runtime_error("Error reading knot vector "+std::to_string(i)+" data");
	}
	
	//Read the axes extents, stored in a single extension HDU.
	{
		//TODO: is this safe when in a shared memory/relocatable mode?
		for (unsigned i = 1; i < ndim; i++)
			extents[i] = &extents[0][2*i];
		
		long n_extents = 0;
		long fpix = 1;
		int ext_error = 0;
		fits_movnam_hdu(fits, IMAGE_HDU, const_cast<char*>("EXTENTS"), 0, &ext_error);
		fits_get_img_size(fits, 1, &n_extents, &ext_error);
		if (n_extents != 2*ndim)
			ext_error = 1;
		
		if (ext_error != 0) { // No extents. Make up some reasonable ones.
			for (unsigned i = 0; i < ndim; i++) {
				extents[i][0] = knots[i][order[i]];
				extents[i][1] = knots[i][nknots[i] - order[i] - 1];
			}
		} else {
			fits_read_pix(fits, TDOUBLE, &fpix, n_extents, NULL,
						  &extents[0][0], NULL, &ext_error);
			if (ext_error!=0)
				throw std::runtime_error("Error reading extent data");
		}
	}

	fill_knot_spacing_bounds();
	
	if(error!=0)
		throw std::runtime_error("Error reading "+filePath+": Error "+std::to_string(error));
	
	return (error==0);
}

template<typename Alloc>
void splinetable<Alloc>::write_fits(const std::string& filePath) const{
	if(ndim==0)
		throw std::runtime_error("splinetable contains no data, cannot write to file");
	
	fitsfile* fits;
	int error = 0;
	
	fits_create_file(&fits, ("!"+filePath).c_str(), &error);
	if (error != 0)
		throw std::runtime_error(("CFITSIO failed to open "+filePath+" for writing").c_str());
	
	struct fits_cleanup{
		fitsfile* fits;
		fits_cleanup(fitsfile* f):fits(f){}
		~fits_cleanup(){
			int error=0;
			fits_close_file(fits, &error);
			fits_report_error(stderr, error);
		}
	} cleanup(fits);
	
	write_fits_core(fits);
}
	
template<typename Alloc>
std::pair<void*,size_t> splinetable<Alloc>::write_fits_mem() const{
	if(ndim==0)
		throw std::runtime_error("splinetable contains no data, cannot write to (memory) file");
	
	fitsfile* fits;
	int error = 0;
	
	const size_t FITS_blocksize=2880;
	size_t memsize=FITS_blocksize;
	void* buf=malloc(memsize);
	
	try{
		fits_create_memfile(&fits, &buf, &memsize, FITS_blocksize, realloc, &error);
		
		struct fits_cleanup{
			fitsfile* fits;
			fits_cleanup(fitsfile* f):fits(f){}
			~fits_cleanup(){
				int error=0;
				fits_close_file(fits, &error);
				fits_report_error(stderr, error);
			}
		} cleanup(fits);
		
		write_fits_core(fits);
	}catch(std::exception& ex){
		throw std::runtime_error("Failed to write FITS memory 'file': \n"+std::string(ex.what()));
	}
	
	return(std::make_pair(buf,memsize));
}
	
template<typename Alloc>
void splinetable<Alloc>::write_fits_core(fitsfile* fits) const{
	int error = 0;
	/*
	 * Write the coefficients
	 * Fits stores arrays in a sort-of Fortran-like way,
	 * so we need to write the axes in reverse order.
	 * Note that the strides will not need to be written explicitly,
	 * as they can be reconstructed from naxes.
	 */
	{
		std::unique_ptr<long[]> naxes(new long[ndim]);
		uint64_t nelements=1;
		for(uint32_t i=0; i<ndim; i++) {
			naxes[i] = this->naxes[ndim - i - 1];
			nelements *= naxes[i];
		}
		fits_create_img(fits, FLOAT_IMG, ndim, naxes.get(), &error);
		if (error != 0)
			throw std::runtime_error("Failed to create FITS image for spline coefficients");
	
		std::unique_ptr<long[]> fpixel(new long[ndim]);
		std::fill_n(fpixel.get(),ndim,1L);
		fits_write_pix(fits, TFLOAT, fpixel.get(), nelements, &coefficients[0], &error);
		if (error != 0)
			throw std::runtime_error("Failed to write coefficients to FITS image");
	}
	
	// Write out header information
	const char typeString[]="Spline Coefficient Table";
	fits_write_key(fits, TSTRING, "TYPE", (void*)&typeString, NULL, &error);
	if (error != 0)
		throw std::runtime_error("Failed to write TYPE key");
	
	char nameBuffer[64];
	// Write spline orders
	for(uint32_t i=0; i<ndim; i++) {
		int chars=snprintf(nameBuffer,sizeof(nameBuffer),"ORDER%d",i);
		if (chars < 0 || unsigned(chars)>=sizeof(nameBuffer))
			throw std::runtime_error("ORDER key too long");
		fits_write_key(fits, TINT, nameBuffer, &order[i], "B-Spline Order", &error);
		if (error != 0)
			throw std::runtime_error("Failed to write ORDER");
	}
	
	// Write periods
	if (periods) {
		for(uint32_t i=0; i<ndim; i++) {
			int chars=snprintf(nameBuffer,sizeof(nameBuffer),"PERIOD%d",i);
			if (chars < 0 || unsigned(chars)>=sizeof(nameBuffer))
				throw std::runtime_error("PERIOD key too long");
			fits_write_key(fits, TDOUBLE, nameBuffer, &periods[i], NULL, &error);
			if (error != 0)
				throw std::runtime_error("Failed to write PERIOD");
		}
	}
	
	// Write 'aux' things, whatever they may be
	// TODO: error checking that these don't collide with anything else?
	for(uint32_t i=0; i<naux; i++) {
		fits_write_key(fits, TSTRING, aux[i][0], aux[i][1], NULL, &error);
		if (error != 0)
			throw std::runtime_error("Failed to write aux entry");
	}
	// done with headers
	
	// Write knot vectors
	for(uint32_t i=0; i<ndim; i++) {
		if(nknots[i]>(uint64_t)std::numeric_limits<long>::max())
			throw std::runtime_error("Too many knots to store in FITS format");
		long axis=nknots[i];
		fits_create_img(fits, DOUBLE_IMG, 1, &axis, &error);
		if (error != 0)
			throw std::runtime_error("Failed to create FITS image for knot vector");
		
		int chars=snprintf(nameBuffer,sizeof(nameBuffer),"KNOTS%d",i);
		if (chars < 0 || unsigned(chars)>=sizeof(nameBuffer))
			throw std::runtime_error("Knot vector name too long");
		fits_update_key(fits, TSTRING, "EXTNAME", nameBuffer, NULL, &error);
		if (error != 0)
			throw std::runtime_error("Failed to set knot vector EXTNAME");
		
		long pixel=1;
		fits_write_pix(fits, TDOUBLE, &pixel, axis, knots[i], &error);
		if (error != 0)
			throw std::runtime_error("Failed to write knot vector");
	}
	
	// Write extents, if they exist
	if (extents) {
		long axis=ndim*2;
		fits_create_img(fits, DOUBLE_IMG, 1, &axis, &error);
		if (error != 0)
			throw std::runtime_error("Failed to create FITS image for extents");
		
		const char extName[]="EXTENTS";
		fits_update_key(fits, TSTRING, "EXTNAME", (void*)&extName, NULL, &error);
		if (error != 0)
			throw std::runtime_error("Failed to set extents EXTNAME");
		
		long pixel=1;
		fits_write_pix(fits, TDOUBLE, &pixel, axis, extents[0], &error);
		if (error != 0)
			throw std::runtime_error("Failed to write extents");
	}
}

} //namespace photospline

#endif //PHOTOSPLINE_FITSIO_H
