// Noise class definition

#ifndef NOISE_H
#define NOISE_H


// standard includes
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

// gsl includes
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_fft_complex.h>


// my headers

// my definitions

class Noise {

public:
	Noise(unsigned lengthin) :
	length(lengthin) 
	{
		doublevec = new double[length];
		cvec = gsl_vector_complex_alloc(length);
	}
	~Noise()
	{
		delete doublevec;
		gsl_vector_complex_free(cvec);
	}

	void setcvec(const double meanamp,const double sigmaamp, const double phasep2p);
	void setdouble(const double mean,const double sigma);

	void addnoise(double * doublevec,const unsigned length);
	void addnoise(gsl_complex_vector cvec);
	inline void addnoise(double * value){
		value *=
	}

private:
	double * doublevec;
	gsl_vector_complex *cvec;
	const unsigned length;
	double amp;
};

#endif 
