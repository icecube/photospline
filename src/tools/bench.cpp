#include <iostream>
#include "photospline/splinetable.h"

void help(){
	std::cerr
	    << "usage: photospline-bench [-h] [--double] spline_file" << std::endl
	    << std::endl
	    << "positional arguments:" << std::endl
	    << " spline_file   path the FITS file" << std::endl
	    << std::endl
	    << "optional arguments:" << std::endl
	    << " -h, --help    show this message and exit" << std::endl
	    << " --double      evaluate in double precision" << std::endl
	;
}

int main(int argc, char* argv[]){
	bool double_precision=false;
	std::string fname;
	if(argc<2){
		help();
		return(1);
	}
	for (int i=1; i < argc; i++) {
		std::string arg(argv[i]);
		if (arg == "-h" || arg == "--help") {
			help();
			return(0);
		} else if (arg == "--double") {
			double_precision=true;
			continue;
		} else if (fname.empty()) {
			fname = arg;
		} else {
			std::cerr << "got unrecognized argument " << arg << std::endl;
			help();
			return(1);
		}
	}
	photospline::splinetable<> spline(fname);
	unsigned dim=spline.get_ndim();
	unsigned int iterations=(8.e7*exp(-(double)dim/1.4427));
	if (double_precision) {
		spline.benchmark_evaluation<double>(iterations,true);
	} else {
		spline.benchmark_evaluation<float>(iterations,true);
	}
}