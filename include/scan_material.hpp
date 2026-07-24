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
// HDF5 includes
#include "hdf5.h"
#include "H5Cpp.h"
#include "H5Attribute.h"
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
