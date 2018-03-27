// MatResponse class definition

#ifndef MATRESPONSE_H
#define MATRESPONSE_H


// standard includes
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

// my headers
#include "Constants.hpp"
#include "Pulse.hpp"
        
// my definitions

using namespace Constants;
class PulseFreq;
                
class MatResponse {

	friend class PulseFreq;
        
public: 
	MatResponse(double t0_in,double width_in,double atten_in,double phase_in); 
	MatResponse(MatResponse & rhs);
	~MatResponse(void){}
	
	
	void aalphabbeta(double ain,double alphain,double bin,double betain){
		a=ain/(ain+bin);
		b=bin/(ain+bin);
		alpha=fsPau<double>()*alphain;
		beta=fsPau<double>()*betain;
	}

	void setstepvec_full(PulseFreq & pulse);
	void setstepvec_full(PulseFreq * pulse);
	void setstepvec_amp(PulseFreq & pulse);
	void setstepvec_amp(PulseFreq * pulse);
	void setstepvec_phase(PulseFreq & pulse);
	void setstepvec_phase(PulseFreq * pulse);
	void addstepvec_amp(PulseFreq & pulse,double delay);
	void addstepvec_amp(PulseFreq * pulse,double delay);
	void addstepvec_phase(PulseFreq & pulse,double delay);
	void addstepvec_phase(PulseFreq * pulse,double delay);
	
	void buffervectors(PulseFreq & pulse);
	void buffervectors(PulseFreq * pulse);

	
	inline void set_delay(double in) { setdelay(in); }
	inline void setdelay(double in){ t0 = in/fsPau<float>(); }
	inline void set_attenuation(double in) { setatten(in); }
	inline void setatten(double in){ attenuation = in; }
	inline void set_phase(double in) { setphase(in); }
	inline void setphase(double in){ phase = in; }


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
