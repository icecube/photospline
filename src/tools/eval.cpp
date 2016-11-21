#include <iostream>
#include <sstream>

#include <photospline/splinetable-mod.h>

int main(int argc, char* argv[]){
	if(argc<3){
		std::cerr << "Usage: photospline-eval spline_file coordinates..." << std::endl;
		return(1);
	}

	photospline::splinetable<> spline(argv[1]);
	if(argc!=spline.get_ndim()+2){
		std::cerr << "Wrong number of coordinates (" << argc-2
		<< ") for evaluation of " << spline.get_ndim() << "-dimensional spline" << std::endl;
		return(1);
	}
	
	std::vector<double> coords(spline.get_ndim());
	std::vector<int> centers(spline.get_ndim());
	{
		std::istringstream ss;
		for(unsigned int i=0; i<spline.get_ndim(); i++){
			ss.clear();
			ss.str(argv[i+2]);
			ss >> coords[i];
			if(ss.fail()){
				std::cerr << "Failed to interpret \"" << argv[i+2]
				<< "\" as a floating-point coordinate value" << std::endl;
				return(1);
			}
		}
	}
	
	bool okay=spline.searchcenters(coords.data(), centers.data());
	if(!okay){
		std::cerr << "Failed to look up knot indices for coordinates" << std::endl;
		return(1);
	}
	
	double evaluate=spline.ndsplineeval(coords.data(), centers.data(), 0);
	
	std::cout << evaluate << std::endl;
	return(0);
}