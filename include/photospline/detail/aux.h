#ifndef PHOTOSPLINE_DETAIL_AUX_H
#define PHOTOSPLINE_DETAIL_AUX_H

#include "photospline/detail/fitsio.h"

namespace photospline{

template<typename Alloc>
const char* splinetable<Alloc>::get_aux_value(const char* key) const{
	const char* value=nullptr;
	for (uint32_t i=0; i < naux; i++) {
		// NB: aux[i][0] may be a smart pointer
		if (strcmp(key, &*aux[i][0]) == 0) {
			value = &*aux[i][1];
			break;
		}
	}
	return (value);
}

template<typename Alloc>
bool splinetable<Alloc>::remove_key(const char* key){
	uint32_t i;
	for (i=0; i < naux; i++) {
		if (strcmp(key, &*aux[i][0]) == 0)
			break;
	}
	if (i==naux) //key was never here
		return (false);
	
	//remove the key
	char_ptr_ptr tmp_aux=nullptr;
	try{
		//first, shuffle all of the remaining keys and values into a temporary buffer
		tmp_aux = new char_ptr[naux-1];
		for (uint32_t j=0, k=0; j<naux; j++) {
			if (j!=i)
				tmp_aux[k++]=aux[j];
		}
		//eliminate the selected key and value
		deallocate(aux[i][0],strlen(&aux[i][0][0])+1);
		deallocate(aux[i][1],strlen(&aux[i][1][0])+1);
		deallocate(aux[i],2);
		//deallocate the old aux
		deallocate(aux,naux);
		//allocate new aux
		naux--;
		//this should be able to fit in the space vacated by the previous version,
		//even if nowhere else is available, so it should not fail under sane
		//circumstances
		aux = allocate<char_ptr_ptr>(naux);
		//copy back remaining keys and values
		std::copy_n(&tmp_aux[0],naux,&aux[0]);
	}catch(...){
		delete[] tmp_aux;
		throw;
	}
	return (true);
}
	
template<typename Alloc>
template<typename T>
bool splinetable<Alloc>::read_key(const char* key, T& result) const{
	const char* value = get_aux_value(key);
	if(!value)
		return (false);
	std::istringstream ss(&*value);
	ss >> result;
	return (!ss.fail());
}
	
template<typename Alloc>
bool splinetable<Alloc>::read_key(const char* key, std::string& result) const{
	const char* value = get_aux_value(key);
	if(!value)
		return (false);
	result=&*value;
	return (true);
}
	
template<typename Alloc>
template<typename T>
bool splinetable<Alloc>::write_key(const char* key, const T& value){
	//check if the key is allowed
	if (reservedFitsKeyword(key))
		throw std::runtime_error("Cannot set key with reserved name "+std::string(key));
	size_t keylen = strlen(key) + 1;
	size_t maxdatalen=68; //valid for short keys
	if(keylen<=9){ //up to 8 bytes of data
		for(size_t i=0; i<keylen-1; i++){
			if(!(std::isupper(key[i]) || std::isdigit(key[i])) || key[i]=='-' || key[i]=='_')
				throw std::runtime_error("Standard (short) FITS header keywords are forbidden "
										 "to contain characters other than uppercase letters, "
										 "digits, dashes, and underscores (key was '"+
										 std::string(key)+"')");
		}
	}
	else{
		//it is unclear what the contraints on the format of long keyword names 
		//are, since the 'HIERARCH Keyword Convention' document refers to "the 
		//rules for free-format keywords, as defined in the FITS Standard 
		//document", when no such rules appear to exist. If this was intended to 
		//refer to section 4.1.2.1 then cfitsio's behavior of allowing long 
		//keywords (not split by spaces or periods) at all is non-conforming anyway. 
		for(size_t i=0; i<keylen-1; i++){
			if(key[i]=='=')
				throw std::runtime_error("Standard (short) FITS header keywords must not "
										 "contain '=' characters (key was '"+
										 std::string(key)+"')");
		}
		maxdatalen=80-(13+keylen-1); //14 characters for "HIERARCH ", "= '", and "'"
	}
	std::ostringstream ss;
	ss << value;
	if(ss.fail())
		return(false);
	std::string valuedata=ss.str();
	size_t valuelen = valuedata.size() + 1;
	//For normal (short) keys, we get up to 68 bytes of storage, but for longer keywords
	//the 'HIERARCH Keyword Convention' kicks in and limits us further
	if(valuelen-1>maxdatalen){
		throw std::runtime_error("Value is too long to be stored as a FITS keyword ('"
								 +valuedata+"' has length "+std::to_string(valuelen-1)
								 +", but a maximum of "+std::to_string(maxdatalen)+
								 " characters will fit with this key since continued "
								 "string keywords are not currently implemented.)");
	}
	//check if the key already exists and we should update it
	size_t i;
	for (i=0; i < naux; i++) {
		if (strcmp(key, &*aux[i][0]) == 0)
			break;
	}
	if (i!=naux) { //the key was found, so update it
		//try to allocate the correct amount of space for the new value
		try{
			char_ptr new_value=allocate<char>(valuelen);
			std::copy(valuedata.begin(),valuedata.end(),new_value);
			*(new_value+valuelen-1)=0;
			deallocate(aux[i][1],strlen(&aux[i][1][0])+1);
			aux[i][1]=new_value;
		}catch(...){
			throw std::runtime_error("Unable to allocate storage for additional aux key");
		}
		return (false);
	}
	//otherwise, allocate more space and append
	else {
		//see if we can get space for the new data before we touch any existing things
		char_ptr new_key=nullptr, new_value=nullptr;
		char_ptr_ptr new_entry=nullptr;
		char_ptr_ptr_ptr new_aux=nullptr;
		try{
			new_aux=allocate<char_ptr_ptr>(naux+1);
			new_entry=allocate<char_ptr>(2);
			new_key=allocate<char>(keylen);
			new_value=allocate<char>(valuelen);
		}catch(...){
			deallocate(new_aux,naux+1);
			deallocate(new_entry,2);
			deallocate(new_key,keylen);
			deallocate(new_value,valuelen);
			throw std::runtime_error("Unable to allocate storage for additional aux key");
		}
		//copy over existing data
		for(size_t j=0; j<naux; j++)
			new_aux[j] = aux[j];
		new_aux[naux] = new_entry;
		new_aux[naux][0] = new_key;
		new_aux[naux][1] = new_value;
		std::copy(key,key+keylen,new_aux[naux][0]);
		std::copy(valuedata.begin(),valuedata.end(),new_value);
		*(new_value+valuelen-1)=0;
		deallocate(aux,naux);
		aux = new_aux;
		naux++;
		return (true);
	}
}
	
} //namespace photospline

#endif
