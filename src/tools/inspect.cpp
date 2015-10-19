#include <iostream>

#include <photospline/splinetable-mod.h>

int main(int argc, char* argv[]){
	if(argc!=2){
		std::cerr << "Usage: photospline-inspect spline_file" << std::endl;
		return(1);
	}
	photospline::splinetable<> spline(argv[1]);
	std::cout << spline.get_ndim() << " dimensional spline" << std::endl;
	std::cout << "Spline orders: ";
	for(size_t i=0; i<spline.get_ndim(); i++)
		std::cout << spline.get_order(i) << ' ';
	std::cout << std::endl;
	std::cout << "Extents: ";
	for(size_t i=0; i<spline.get_ndim(); i++)
		std::cout << '[' << spline.lower_extent(i) << ',' << spline.upper_extent(i) << "] ";
	std::cout << std::endl;
	std::cout << "Knot counts: ";
	for(size_t i=0; i<spline.get_ndim(); i++)
		std::cout << spline.get_nknots(i) << ' ';
	std::cout << std::endl;
	if(spline.get_naux_values()){
		std::cout << "Auxilliary keys:\n";
		for(size_t i=0; i<spline.get_naux_values(); i++)
			std::cout << "  " << spline.get_aux_key(i) << " = " << spline.get_aux_value(spline.get_aux_key(i)) << std::endl;
	}
}