// MatResponse class definition

#ifndef MATRESPONSE_H
#define MATRESPONSE_H


// standard includes
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

// gsl includes
#include <gsl/gsl_const_num.h> 
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_fft_complex.h>

// my headers
#include "Constants.hpp"
#include "Pulse.hpp"
        
// my definitions

using namespace Constants;
                
class MatResponse {
        
public: 
	MatResponse(double t0_in=0.0,double width_in=10.0,double atten_in = 0.95,double phase_in = 0.03) :
	t0(t0_in / fsPau<float>()),
	twidth(width_in * M_SQRTPI / fsPau<float>() / 2.0),
	attenuation(atten_in),
	phase(phase_in)
	{     
		a=0.75;
		b=0.25;
		alpha=5.0e-4*fsPau<float>();
		beta=0.0;
	}        

	~MatResponse(void){}
	
	
	void aalphabbeta(double ain,double alphain,double bin,double betain){
		a=ain/(ain+bin);
		b=bin/(ain+bin);
		alpha=fsPau<float>()*alphain;
		beta=fsPau<float>()*betain;
	}

	void setstepvec_full(PulseFreq & pulse);
	void setstepvec_full(PulseFreq * pulse);
	void setstepvec_amp(gsl_vector * modamp,double dt);
	void setstepvec_phase(gsl_vector * modphase,double dt);
	void addstepvec_amp(gsl_vector * modamp,double dt,double delay);
	void addstepvec_phase(gsl_vector * modphase,double dt,double delay);
	void setstepvec_amp(PulseFreq & pulse);
	void setstepvec_amp(PulseFreq * pulse);
	void setstepvec_phase(PulseFreq & pulse);
	void setstepvec_phase(PulseFreq * pulse);
	void addstepvec_amp(PulseFreq & pulse,double delay);
	void addstepvec_amp(PulseFreq * pulse,double delay);
	void addstepvec_phase(PulseFreq & pulse,double delay);
	void addstepvec_phase(PulseFreq * pulse,double delay);
	
	void buffervectors(gsl_vector * modamp,gsl_vector * modphase,double dt);
	void buffervectors(PulseFreq & pulse);
	void buffervectors(PulseFreq * pulse);

	void setstepvec_marco(gsl_vector * modamp,double dt);
	
	inline void setdelay(double tin){ t0 = tin/fsPau<float>(); }
	inline void setatten(double attenin){ attenuation = attenin; }
	inline void setphase(double phasein){ phase = phasein; }


        inline void setreflectance(const double reflectancein){reflectance = reflectancein;}
        inline void setetalondelay(const double etalondelayin){etalondelay = etalondelayin;}

	inline double getetalondelay(void){return etalondelay;}
	inline double getetalondelayinau(void){return etalondelay/fsPau<float>();}
	inline double getreflectance(void){return reflectance;}



private:
	
	double t0,twidth,attenuation,phase,a,alpha,b,beta;
	double etalondelay, reflectance;


};

#endif
