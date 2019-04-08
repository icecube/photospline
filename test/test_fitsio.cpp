#include "test.h"
#include "photospline/splinetable.h"
#include <unistd.h>

TEST(read_fits_spline){
	photospline::splinetable<> spline("test_data/test_spline_4d.fits");
	
	ENSURE_EQUAL(spline.get_ndim(),4u,"Test spline should have dimension 3");
	
	for(uint32_t i=0; i<spline.get_ndim(); i++)
		ENSURE_EQUAL(spline.get_order(i),2u,"Test spline should be order 2 in all dimensions");
	
	for(uint32_t i=0; i<spline.get_ndim(); i++)
		ENSURE_EQUAL(spline.get_nknots(i),uint64_t(16)-2*i);
}

void compare_splines(const photospline::splinetable<>& spline1,
					 const photospline::splinetable<>& spline2){
	ENSURE_EQUAL(spline2.get_ndim(),spline1.get_ndim(),"Splines shold have same dimension");
	
	for(uint32_t i=0; i<spline1.get_ndim(); i++)
		ENSURE_EQUAL(spline2.get_order(i),spline1.get_order(i),"Splines should have same orders");
	
	for(uint32_t i=0; i<spline1.get_ndim(); i++){
		ENSURE_EQUAL(spline2.get_nknots(i),spline1.get_nknots(i),"Splines should have same number of knots");
		for(uint64_t j=0; j<spline1.get_nknots(i); j++)
			ENSURE_EQUAL(spline2.get_knot(i,j),spline1.get_knot(i,j),"Splines should have same knots");
	}
	
	for(uint32_t i=0; i<spline1.get_ndim(); i++){
		ENSURE_EQUAL(spline2.lower_extent(i),spline1.lower_extent(i),"Splines should have same extents");
		ENSURE_EQUAL(spline2.upper_extent(i),spline1.upper_extent(i),"Splines should have same extents");
	}
	
	for(uint32_t i=0; i<spline1.get_ndim(); i++)
		ENSURE_EQUAL(spline2.get_period(i),spline1.get_period(i),"Splines should have same periods");
	
	for(uint32_t i=0; i<spline1.get_ndim(); i++)
		ENSURE_EQUAL(spline2.get_ncoeffs(i),spline1.get_ncoeffs(i),"Splines should have same number of coefficients");
	
	ENSURE_EQUAL(spline2.get_ncoeffs(),spline1.get_ncoeffs(),"Splines should have same total number of coefficients");
	for(uint64_t i=0; i<spline1.get_ncoeffs(); i++)
		ENSURE_EQUAL(spline2.get_coefficients()[i],spline1.get_coefficients()[i],"Splines should have same coefficients");
	
	for(uint32_t i=0; i<spline1.get_ndim(); i++)
		ENSURE_EQUAL(spline2.get_stride(i),spline1.get_stride(i),"Splines should have same strideas");
}

TEST(write_fits_spline){
	photospline::splinetable<> spline("test_data/test_spline_4d.fits");
	spline.write_key("SHORTKEY",123);
	try {
		// lowercased names should throw
		spline.write_key("ALongerKey", "a string of text");
		throw std::logic_error("Should have thrown");
	} catch (std::runtime_error &) {}
	spline.write_key("ALONGERKEY","a string of text");
	
	spline.write_fits("write_test_spline.fits");
	photospline::splinetable<> spline2("write_test_spline.fits");
	
	compare_splines(spline,spline2);
	
	int extra_key=0;
	ENSURE(spline2.read_key("SHORTKEY",extra_key));
	ENSURE_EQUAL(extra_key,123,"Additional keys must survive FITS serialization");
	std::string extra_key2;
	ENSURE(spline2.read_key("ALONGERKEY",extra_key2));
	ENSURE_EQUAL(extra_key2,"a string of text","Additional keys must survive FITS serialization");
	
	unlink("write_test_spline.fits");
}

TEST(fits_mem_spline){
	photospline::splinetable<> spline("test_data/test_spline_4d.fits");
	spline.write_key("SHORTKEY",123);
	spline.write_key("ALONGERKEY","a string of text");
	
	auto buffer=spline.write_fits_mem();
	std::unique_ptr<void,void(*)(void*)> data(buffer.first,&free);
	size_t buffer_size=buffer.second;
	photospline::splinetable<> spline2;
	spline2.read_fits_mem(data.get(),buffer_size);
	
	compare_splines(spline,spline2);
	
	int extra_key=0;
	ENSURE(spline2.read_key("SHORTKEY",extra_key));
	ENSURE_EQUAL(extra_key,123,"Additional keys must survive FITS serialization");
	std::string extra_key2;
	ENSURE(spline2.read_key("ALONGERKEY",extra_key2));
	ENSURE_EQUAL(extra_key2,"a string of text","Additional keys must survive FITS serialization");
}