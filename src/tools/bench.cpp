#include <iostream>
#include "photospline/splinetable.h"

int main(int argc, char* argv[]){
	if(argc<2){
		std::cerr << "Usage: photospline-bench spline_file" << std::endl;
		return(1);
	}
	photospline::splinetable<> spline(argv[1]);
	unsigned dim=spline.get_ndim();
	unsigned int iterations=(8.e7*exp(-(double)dim/1.4427));
	spline.benchmark_evaluation(iterations,true);
}