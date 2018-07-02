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
