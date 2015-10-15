#include <map>
#include <string>
#include <sstream>

void emit_error(const std::string& file, size_t line,
				const std::string& criterion, const std::string& message="");

#define ENSURE(cond,...) \
	do{ \
		if(!(cond)) \
			emit_error(__FILE__,__LINE__,#cond,##__VA_ARGS__); \
	}while(0)

template<typename T>
std::string express_comparison(const std::string& e1, const T& v1,
                               const std::string& e2, const T& v2,
                               const T& tolerance=0){
	std::ostringstream ss;
	ss << v1 << " (" << e1 << ") != " << v2 << " (" << e2 << ")";
	if(tolerance!=0)
		ss << " to within " << tolerance;
	return(ss.str());
}

#define ENSURE_EQUAL(first,second,...) \
	do{ \
		if(!((first)==(second) || \
             (std::isnan(first) && std::isnan(second)))) \
			emit_error(__FILE__,__LINE__, \
			  express_comparison(#first,first,#second,second),##__VA_ARGS__); \
	}while(0)

#define ENSURE_DISTANCE(first,second,tolerance,...) \
do{ \
	if(!(std::abs((first)-(second))<(tolerance))) \
		emit_error(__FILE__,__LINE__, \
		  express_comparison(#first,first,#second,second,tolerance),##__VA_ARGS__); \
}while(0)

#define FAIL(...) \
	emit_error(__FILE__,__LINE__,"FAIL",##__VA_ARGS__)

extern std::map<std::string,void(*)()> test_registry;

struct register_test{
	register_test(const std::string& test_name, void(*test)()){
		test_registry.insert(std::make_pair(test_name,test));
	}
};

#define TEST(name) \
	void test_func ## name (); \
	static register_test register_ ## name (#name,&test_func ## name); \
	void test_func ## name ()
