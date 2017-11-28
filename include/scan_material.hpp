#ifndef SCAN_MATERIAL

#define SCAN_MATERIAL

// standard includes
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <ctime>
#include <complex>
#include <algorithm>
#include <vector>

// OpenCV includes
//#include <opencv/core.h>

// gsl includes
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_fft_complex.h>

// my includes
#include "Constants.hpp"
#include "Pulse.hpp"
#include "MatResponse.hpp"
#include "FiberBundle.hpp"

// parallel includes
#include <omp.h>

/* Defining global variables */
/* Prototyping functions */
/*
void setstepvec_amp(gsl_vector * modamp,double stepwidth,double steptime,double attenuation,double dt);
void setstepvec_phase(gsl_vector * modphase,double stepwidth,double steptime,double maxphase,double dt);
void setstepvec_marco(gsl_vector * modamp,double steptime,double dt);
void setstepvec_fastrise(gsl_vector * modamp,double stepwidth,double steptime,double attenuation,double dt);
void setstepvec_blip(gsl_vector * modamp,double stepwidth,double steptime,double attenuation,double dt);
*/

#endif
