#ifndef FIBERBUNDLE

#define FIBERBUNDLE

#define _USE_MATH_DEFINES
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <complex>
#include <iostream>
#include <fstream>

#include "DataOps.hpp"

class CalibMat {
	public:
		CalibMat(size_t n);
		bool print_mapping(std::ofstream & out);
		inline void set_fsWindow(const double x = 3333){fsWindow = x;}

		inline double delay(const size_t i){ return fsWindow * double(i - ndelays/2); }
		inline size_t get_ndelays(void){return ndelays;}
		inline size_t set_ndelays(const size_t n){ndelays = n; return ndelays;}

	private:
		size_t ndelays;
		double fsWindow;
};

#endif

