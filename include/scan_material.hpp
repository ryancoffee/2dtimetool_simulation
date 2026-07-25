#ifndef SCAN_MATERIAL

#define SCAN_MATERIAL

// standard includes
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <cstddef> // for nullptr
#include <cstdint> // for writing as an int32_t and int16_t
#include <math.h>
#include <string>
#include <ctime>
#include <complex>
#include <algorithm>
#include <vector>
#include <memory>
#include <random>
#include <chrono>

/* Commenting out HDF5 for sake of flight, need to build hdf5 build environment for laptop
*/
// HDF5 includes
#include "hdf5.h"
#include "H5Cpp.h"
#include "H5Attribute.h"

/*
H5::IntType h5uint16( H5::PredType::NATIVE_USHORT );
H5::IntType h5int16( H5::PredType::NATIVE_SHORT );
H5::IntType h5uint32( H5::PredType::NATIVE_UINT );
H5::IntType h5int32( H5::PredType::NATIVE_INT );
H5::StrType h5string(0, H5T_VARIABLE);
h5uint16.setOrder( H5T_ORDER_LE );
h5int16.setOrder( H5T_ORDER_LE );
h5uint32.setOrder( H5T_ORDER_LE );
h5int32.setOrder( H5T_ORDER_LE );
*/

// FFTW
#include <fftw3.h>

// my includes
#include "Constants.hpp"
#include "Pulse.hpp"
#include "MatResponse.hpp"
#include "FiberBundle.hpp"
#include "CalibMat.hpp"

// parallel includes
#include <omp.h>



#endif
