#include <iostream>
#include <stdexcept>
#include <unistd.h>
#include <vector>

#include "test.h"

struct test_exception : public std::runtime_error{
	test_exception(const std::string& msg):std::runtime_error(msg){}
};

void emit_error(const std::string& file, size_t line,
				const std::string& criterion, const std::string& message){
	std::ostringstream ss;
	ss << file << ':' << line << "\n\t";
	if(message.empty())
		ss << "Assertion failed: \n";
	else
		ss << message << ": \n";
	ss << '\t' << criterion << std::endl;
	throw test_exception(ss.str());
}

std::map<std::string,void(*)()>&
test_registry()
{
	static std::map<std::string,void(*)()> *registry = new std::map<std::string,void(*)()>;
	return *registry;
}

int main(int argc, char* argv[]){
	std::vector<std::string> test_filters;
	for(int i=1; i<argc; i++){
		std::string arg=argv[i];
		if(arg=="WORKING_DIRECTORY"){
			if(i+1>=argc){
				std::cerr << "WORKING_DIRECTORY not specified" << std::endl;
				return(1);
			}
			if(chdir(argv[i+1])!=0){
				std::cerr << "Failed to change working directory to " << argv[i+1] << std::endl;
				return(1);
			}
			i++;
		} else if (arg=="-k"){
			if(i+1>=argc){
				std::cerr << "-k requires an argument" << std::endl;
				return(1);
			}
			test_filters.push_back(argv[i+1]);
			i++;
		}
	}
	
	std::cout << "Running " << test_registry().size() << " tests" << std::endl;
	bool all_pass=true;
	size_t passes=0, failures=0;
	for(std::map<std::string,void(*)()>::const_iterator test=test_registry().begin();
		test!=test_registry().end(); test++){
		if (!test_filters.empty()) {
			bool target=false;
			for(auto &f : test_filters) {
				if (test->first.find(f) != std::string::npos) {
					target=true;
					break;
				}
			}
			if (!target)
				continue;
		}
		bool pass=false;
		std::cout << test->first << ": ";
		std::cout.flush();
		try{
			(test->second)();
			pass=true;
		}catch(test_exception& ex){
			std::cout << "FAIL\n " << ex.what() << std::endl;
		}catch(std::exception& ex){
			std::cout << "FAIL\n Exception: " << ex.what() << std::endl;
		}catch(...){
			std::cout << "FAIL\n Unknown object thrown" << std::endl;
		}
		if(pass)
			std::cout << "PASS" << std::endl;
		(pass?passes:failures)++;
		all_pass &= pass;
	}
	std::cout << passes << " test" << (passes!=1?"s":"") << " pass"
	<< (passes!=1?"":"es") << ", "
	<< failures << " fail" << (failures!=1?"":"s") << std::endl;
	return(all_pass ? 0 : 1);
}