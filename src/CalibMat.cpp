#include "CalibMat.hpp"

#include "Constants.hpp"
using Constants::pi;
using Constants::pi_3;


CalibMat::CalibMat(size_t n = 109,double win = 3000)
: ndelays(n)
, fsWindow(win)
{
	std::cerr << "In constructor CalibMat() " << std::endl;
}
bool CalibMat::print_delays(ofstream & out)
{
	if (!out.is_open()){
		std::cerr << "Passing unopened ofstream reference" << std::endl;
		return false;
	}
	out << "# Calibration Matrix\n";
	for (size_t i =0;i<ndelays;++i){
		out << delay(i) << "\n";
	}
	out << std::flush;
	return true;
}
