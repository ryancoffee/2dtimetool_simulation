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
#include <Constants.hpp>
#include <Pulse.hpp>
        
// my definitions

using namespace Constants;
class PulseFreq;

class CarriersNotSet : public exception
{
	virtual const char* what() const throw()
	{
		return "Carriers were not yet set";
	}
};
                
class MatResponse {

	friend class PulseFreq;
        
public: 
	MatResponse(double t0_in,double width_in,double atten_in,double phase_in); 
	MatResponse(MatResponse & rhs);
	~MatResponse(void){}
	
	MatResponse & operator=(const MatResponse & rhs);
	
	void aalphabbeta(double ain,double alphain,double bin,double betain){
		a=ain/(ain+bin);
		b=bin/(ain+bin);
		alpha=fsPau<double>()*alphain;
		beta=fsPau<double>()*betain;
	}

	inline double thermaletalondelay(const double TinK = 300.){
		return n_refractive(TinK)/C_umPfs<double>()*l_thickness_um * 2;
	}
	inline double n_refractive(const double TinK = 300.){
		return std::exp( (std::log(TinK) - n_b ) / n_a ) + n_0;
	}
	void set_n_refractive(const double n_0in, const double n_ain, const double n_bin);
	void set_thickness(const double lin);
	
 	inline double bandgap(double in=4.47){bandgap_eV = in; return bandgap_eV;}
	inline double bandgap(void){return bandgap_eV;}

	void setstepvec_amp(PulseFreq & pulse,double delay_in = 0.,double scale_in = 1.);
	void setstepvec_amp(PulseFreq * pulse,double delay_in = 0.,double scale_in = 1.);
	void setstepvec_phase(PulseFreq & pulse,double delay_in = 0.,double scale_in = 1.);
        void setstepvec_phase(PulseFreq * pulse,double delay_in = 0.,double scale_in = 1.);
        void addstepvec_amp(PulseFreq & pulse,double delay_in = 0.,double scale_in = 1.);
        void addstepvec_amp(PulseFreq * pulse,double delay_in = 0.,double scale_in = 1.);
	void addstepvec_phase(PulseFreq & pulse,double delay_in = 0.,double scale_in = 1.);
	void addstepvec_phase(PulseFreq * pulse,double delay_in = 0.,double scale_in = 1.);
	void setstepvec_both(PulseFreq * pulse,double delay_in = 0.,double scale_in = 1.);
	void setstepvec_both(PulseFreq & pulse,double delay_in = 0.,double scale_in = 1.);
	void addstepvec_both(PulseFreq * pulse,double delay_in = 0.,double scale_in = 1.);
	void addstepvec_both(PulseFreq & pulse,double delay_in = 0.,double scale_in = 1.);


	bool fill_carriersvec(PulseFreq * pulse,std::ifstream * instream);
	bool fill_carriersvec(PulseFreq & pulse,std::ifstream & instream);
	bool fill_carriersvec(PulseFreq * pulse,double energy_keV);
	bool fill_carriersvec(PulseFreq & pulse,double energy_keV);
	bool setstepvec_both_carriers(PulseFreq * pulse,double delay_in = 0.,double scale_in = 1.);
	bool setstepvec_both_carriers(PulseFreq & pulse,double delay_in = 0.,double scale_in = 1.);
	bool addstepvec_both_carriers(PulseFreq * pulse,double delay_in = 0.,double scale_in = 1.);
	bool addstepvec_both_carriers(PulseFreq & pulse,double delay_in = 0.,double scale_in = 1.);
	
	void buffervectors(PulseFreq & pulse);
	void buffervectors(PulseFreq * pulse);

	
	inline void set_scale(const double in) { setscale(in); }
	inline void setscale(const double in) { scale = in; }
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
	
	bool carriers_set;
	double t0,twidth,attenuation,phase,a,alpha,b,beta,scale;
	double etalondelay, reflectance;
	double bandgap_eV;
	double n_0,n_a,n_b,l_thickness_um;
        std::vector<double> carriers;
};


#endif
