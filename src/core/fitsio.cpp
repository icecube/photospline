#include "../include/photospline/splinetable.h"

namespace photospline{

std::vector<uint32_t> readOrder(fitsfile* fits, uint32_t ndim){
	int error = 0;
	std::vector<uint32_t> order(ndim);
	//See if there is a single order value
	fits_read_key(fits, TINT, "ORDER", &order[0], NULL, &error);
	if (error != 0) {
		error = 0;
		
		//There is not, so look for a separate order in each dimension
		for (int i = 0; i < ndim; i++) {
			std::ostringstream ss;
			ss << "ORDER" << i;
			fits_read_key(fits, TUINT, ss.str().c_str(), &order[i], NULL, &error);
			if (error != 0) {
				throw std::runtime_error("Needs real error message 6");
			}
		}
	} else {
		//all orders are the same
		std::fill(order.begin()+1,order.end(),order[0]);
	}
	return (order);
}
	
bool reservedFitsKeyword(const char* key){
	return(strncmp("BITPIX", key, 6) == 0 ||
	       strncmp("SIMPLE", key, 6) == 0 ||
	       strncmp("TYPE", key, 4) == 0 ||
	       strncmp("ORDER", key, 5) == 0 ||
	       strncmp("NAXIS", key, 5) == 0 ||
	       strncmp("PERIOD", key, 6) == 0 ||
	       strncmp("EXTEND", key, 6) == 0 ||
	       strncmp("COMMENT", key, 7) == 0);
}

uint32_t countAuxKeywords(fitsfile* fits){
	int nkeys = 0, error = 0;
	fits_get_hdrspace(fits, &nkeys, NULL, &error);
	if (nkeys == 0)
		return (0);
	char key[FLEN_KEYWORD], value[FLEN_VALUE];
	int keylen, valuelen;
	
	// figure out how many generic keys there are
	uint32_t naux = 0;
	for (int j = 1 ; j-1 < nkeys; j++) {
		error = 0;
		fits_read_keyn(fits, j, key, value, NULL, &error);
		if (error != 0)
			continue;
		if (reservedFitsKeyword(key))
			continue;
		naux++;
	}
	return (naux);
}

} //namespace photospline
