#include <cassert>
#include <chrono>
#include <iostream>
#include <random>

#include <photospline/splinetable-mod.h>

#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/allocators/allocator.hpp>

volatile double sink;

int main(){
	const std::string tablePath="../../photon_tables/ems_mie_z20_a10.prob.fits";
	
	size_t memEst;
	bool convolve=false;
	if(convolve)
		memEst = photospline::splinetable<>::estimateMemory(tablePath,3/*nknots*/,0/*dim*/);
	else
		memEst = photospline::splinetable<>::estimateMemory(tablePath);
	std::string labels[]={"bytes","KB","MB","GB","TB","PB"};
	std::cout << "Estimated memory required: " << memEst/pow(1024.,floor(log(memEst)/log(1024.))) << ' ' << labels[(int)floor(log(memEst)/log(1024.))] << std::endl;
	
	using namespace boost::interprocess;
	struct shm_remove{
		//shm_remove() { shared_memory_object::remove("TestSharedMemory"); }
		/*~shm_remove(){
			bool removed=shared_memory_object::remove("TestSharedMemory");
			if(removed)
				std::cout << "Removed shared memory space" << std::endl;
		}*/
	} remover;
	managed_shared_memory shm(open_or_create, "TestSharedMemory", memEst);
	
	typedef allocator<void, managed_shared_memory::segment_manager> SharedAllocator;
	SharedAllocator allocator(shm.get_segment_manager());
	photospline::splinetable<SharedAllocator>* shared_table = shm.find_or_construct<photospline::splinetable<SharedAllocator>>("testSplinetable")(tablePath,allocator);
	
	if(convolve){
		double conv_knots[3]={-.01,0,.01};
		std::cout << "Convolving. . . " << std::endl;
		shared_table->convolve(0/*dim*/,conv_knots,3/*nknots*/);
		std::cout << " Done" << std::endl;
	}
	
	//'create' the spline again
	photospline::splinetable<SharedAllocator>* shared_table2 = shm.find_or_construct<photospline::splinetable<SharedAllocator>>("testSplinetable")(tablePath,allocator);
	
	photospline::splinetable<SharedAllocator>::fast_evaluation_token
	tok1=shared_table->get_evaluation_token(),
	tok2=shared_table->get_evaluation_token();
	
	//make sure tables do the same thing
	std::cout << "Comparing tables" << std::endl;
	std::mt19937 rng;
	rng.seed(29);
	
	std::vector<std::uniform_real_distribution<>> dists;
	for(size_t i=0; i<shared_table->get_ndim(); i++)
		dists.push_back(std::uniform_real_distribution<>(shared_table->lower_extent(i),shared_table->upper_extent(i)));
	
	std::vector<double> coords(shared_table->get_ndim());
	std::vector<int> centers1(shared_table->get_ndim()), centers2(shared_table2->get_ndim());
	for(size_t i=0; i<1e4; i++){
		for(size_t j=0; j<shared_table->get_ndim(); j++)
			coords[j]=dists[j](rng);
		
		if(!shared_table->searchcenters(coords.data(), centers1.data())){
			std::cout << "center lookup failure" << std::endl;
			continue;
		}
		
		if(!shared_table2->searchcenters(coords.data(), centers2.data())){
			std::cout << "center lookup failure" << std::endl;
			continue;
		}
		
		assert(std::equal(centers1.begin(),centers1.end(),centers2.begin()));
		
		double e1=shared_table->ndsplineeval(coords.data(), centers1.data(), 0, tok1);
		double e2=shared_table2->ndsplineeval(coords.data(), centers1.data(), 0, tok2);
		assert(e1==e2);
	}
	
	std::cout << "Waiting 30 seconds" << std::endl;
	sleep(30);
	std::cout << "Done" << std::endl;
}