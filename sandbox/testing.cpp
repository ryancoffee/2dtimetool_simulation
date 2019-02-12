#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <iterator>
#include <DataOps.hpp>
#include <cstdint> 
#include <limits>

using namespace DataOps;

int main(void){
        std::vector<double> v = {-1.,3.,4.,10.,-100.,101.};
	std::cout << v;
	auto bounds = std::minmax_element(v.begin(),v.end());
	std::cout << "min = " << *bounds.first << " and max = " << *bounds.second << "\n" << std::flush;
	double scale;
	long int max = std::numeric_limits<int8_t>::max();
	long int min = std::numeric_limits<int8_t>::min();
	if (*bounds.second>std::abs(*bounds.first)){
		scale = max/ *bounds.second;
	} else {
		scale = min/ *bounds.first;
	}
	std::vector<int> result(v.size());
	std::transform(v.begin(),v.end(),result.begin(),[scale](double x){return int(scale*x);});
	std::cout << result;
	std::cout << "max = " << max << " and min = " << min << std::endl;


        return 0;
}

