#ifndef CALIBMAT

#define CALIBMAT

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
		CalibMat(size_t n,double window);
		bool print_delays(std::ofstream & out);
		inline void set_fsWindow(const double x = 3333){fsWindow = x;}

		inline double get_delay(const size_t i){ return fsCenter + fsWindow * (double(i)/double(ndelays)-0.5); }
		inline size_t get_ndelays(void){return ndelays;}
		inline size_t set_ndelays(const size_t n){ndelays = n; return ndelays;}
		inline void set_center(const double x = -3000){fsCenter = x;}

	private:
		size_t ndelays;
		double fsWindow;
		double fsCenter;
};

#endif

