// Pulse class definition

#ifndef PULSE_H
#define PULSE_H


// standard includes
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <random>
#include <fftw3.h>


// my headers
#include <MatResponse.hpp>
#include <Constants.hpp>
#include <ScanParams.hpp>
#include <DataOps.hpp>

// my definitions
#define SAMPLEROUND 1000

#define NU_LOW 0.3
#define NU_HIGH 0.7
#define LAMBDA_LOW 400
#define LAMBDA_HIGH 900

using namespace Constants;

class ScanParams;

class PulseTime {

	public:
		PulseTime(double strength_in = 1e-3 * 0.696, double width_in = 50, double t0_in = 0.0) : 
			strength(strength_in * auenergy<double>()/Eh<double>() * std::pow(aufor10PW<double>(),int(2))), 
			Ctau(width_in * root_pi<double>() / fsPau<double>() / 2.0),
			t0(t0_in / fsPau<double>())
	{
		//    std::clog << "Creating Pulse " << this << std::endl;
	}
		~PulseTime()
		{
			//    std::clog << "Destroying Pulse " << this << std::endl;
		}

		void setstrength(const double in);
		void setwidth(const double in);
		void sett0(const double in);

		inline double getstrength() { return strength; }
		inline double getCtau() { return Ctau; }
		inline double gett0() { return t0; }

		inline bool getenvelope(const double t,double *FF,double *dFFdt) 
		{
			if ( inpulse(t) ){
				*FF = strength * ( std::pow( cos(half_pi<double>()*(t-t0)/Ctau) , int(2)) );
				*dFFdt = -strength/2 * ( pi<double>()/Ctau * sin(pi<double>()*(t-t0)/Ctau));
				return true;
			} else {
				*FF = 0.0;
				*dFFdt = 0.0;
				return false;
			}
		}
		inline bool getenvelope(const double t,double *FF) 
		{
			if ( inpulse(t) ){
				*FF = strength * ( std::pow( cos(half_pi<double>()*(t-t0)/Ctau), int(2) ) );
				return true;
			} else {
				*FF = 0.0;
				return false;
			}
		}

	private:
		double strength, Ctau, t0;

		inline bool inpulse(const double t) 
		{
			if (t >= -Ctau && t <= Ctau){
				return true;
			} else {
				return false;
			}
		}
};




class PulseFreq {

	friend class MatResponse;

	public:
		PulseFreq & operator=(PulseFreq const & rhs); // assignment
		PulseFreq & operator+=(const PulseFreq &rhs); // function definitions must follow the PulseFreq definition
		PulseFreq & operator-=(const PulseFreq &rhs); // function definitions must follow the PulseFreq definition
		PulseFreq & operator*=(const PulseFreq &rhs); // function definitions must follow the PulseFreq definition
		PulseFreq & operator/=(const PulseFreq &rhs); // function definitions must follow the PulseFreq definition
		PulseFreq & operator*=(const double scale);
		PulseFreq & diffamps(const PulseFreq &rhs);
		PulseFreq & normamps(const PulseFreq &rhs);
		PulseFreq & interfere(const PulseFreq &rhs);

	public:
		PulseFreq(const double omcenter_in,const double omwidth_in,const double omonoff_in, double tspan_in);
		PulseFreq(PulseFreq &rhs); // copy constructor
		~PulseFreq(void);

		void set_lamsamples(size_t in);
		void set_gain(int in);
		void set_saturate(int in);

		bool addrandomphase();

		inline unsigned getsamples(void) {
			return samples;
		}
		inline unsigned getdt(void) {
			return dtime;
		}
		inline void fft_totime(void) {

			fftw_execute(FTplan_backward);
			DataOps::mul(cvec,1./std::sqrt(samples),samples);
			cvec2rhophi();
			infreq = false;
			intime = true;
		}
		inline void fft_tofreq(void) {
			if (infreq)
				std::cerr << "died here at fft_tofreq()" << std::endl;
			fftw_execute(FTplan_forward);
			DataOps::mul(cvec,1./std::sqrt(samples),samples);
			cvec2rhophi();
			infreq=true;
			intime=false;
		}

		int addchirp(double chirp_in) {
			if (intime){
				std::cerr << "whoops, trying to add phase in the time domain" << std::endl;
				return 1;
			}
			phase_GDD = chirp_in;
			addGDDtoindex(0,1);
			addGDDtoindex(samples/2,-1);
			for (unsigned i = 1; i<samples/2;i++){
				addGDDtoindex(i,1);
				addGDDtoindex(samples-i,-1);
			}
			return 0;
		}

		int addchirp(double* chirp_in) {
			if (intime){
				std::cerr << "whoops, trying to add phase in the time domain" << std::endl;
				return 1;
			}
			phase_GDD = chirp_in[0];
			phase_TOD = chirp_in[1];
			phase_4th = chirp_in[2];
			phase_5th = chirp_in[3];
			addGDDtoindex(0,1);
			addGDDtoindex(samples/2,-1);
			addTODtoindex(0,1);
			addTODtoindex(samples/2,-1);
			add4thtoindex(0,1);
			add4thtoindex(samples/2,-1);
			add5thtoindex(0,1);
			add5thtoindex(samples/2,-1);
			for (unsigned i = 1; i<samples/2;i++){
				addGDDtoindex(i,1);
				addGDDtoindex(samples-i,-1);
				addTODtoindex(i,1);
				addTODtoindex(samples-i,-1);
				add4thtoindex(i,1);
				add4thtoindex(samples-i,-1);
				add5thtoindex(i,1);
				add5thtoindex(samples-i,-1);
			}
			return 0;
		}

		int addchirp(std::vector<double> & chirp_in) {
			if (intime){
				std::cerr << "whoops, trying to add phase in the time domain" << std::endl;
				return 1;
			}
			assert(chirp_in.size()==4);
			phase_GDD = chirp_in[0];
			phase_TOD = chirp_in[1];
			phase_4th = chirp_in[2];
			phase_5th = chirp_in[3];
			addGDDtoindex(0,1);
			addGDDtoindex(samples/2,-1);
			addTODtoindex(0,1);
			addTODtoindex(samples/2,-1);
			add4thtoindex(0,1);
			add4thtoindex(samples/2,-1);
			add5thtoindex(0,1);
			add5thtoindex(samples/2,-1);
			for (unsigned i = 1; i<samples/2;i++){
				addGDDtoindex(i,1);
				addGDDtoindex(samples-i,-1);
				addTODtoindex(i,1);
				addTODtoindex(samples-i,-1);
				add4thtoindex(i,1);
				add4thtoindex(samples-i,-1);
				add5thtoindex(i,1);
				add5thtoindex(samples-i,-1);
			}
			return 0;
		}


		void attenuate(double attenfactor);
		void delay(double delayin); // expects delay in fs
		void phase(double phasein); // expects delay in units of pi , i.e. 1.0 = pi phase flip 

		inline int modulateamp_time(const std::vector<double> & modulation) {
			if (infreq){
				std::cerr << "whoops, trying time modulation but in frequency domain" << std::endl;
				return 1; }
			if (modulation.size() < samples){
				std::cerr << "size mismatch, out of range in modulateamp_time()" << std::endl;
				return 2;
			}
			modampatindx(0,modulation);
			modampatindx(samples/2,modulation);
			for (unsigned i = 1;i<samples/2;i++){
				modampatindx(i,modulation);
				modampatindx(samples-i,modulation);
			}
			return 0;
		}
		inline int modulateamp_time(void) {
			if (infreq){
				std::cerr << "whoops, trying time modulation but in frequency domain" << std::endl;
				return 1; }
			modampatindx(0);
			modampatindx(samples/2);
			for (unsigned i = 1;i<samples/2;i++){
				modampatindx(i);
				modampatindx(samples-i);
			}
			return 0;
		}

		inline int modulatephase_time(const std::vector<double> & modulation) {
			if (infreq){
				std::cerr << "whoops, trying time modulation but in frequency domain" << std::endl;
				return 1;
			}
			if (modulation.size() < samples){
				std::cerr << "size mismatch, out of range in modulateamp_time()" << std::endl;
				return 2;
			}
			modphaseatindx(0,modulation);
			modphaseatindx(samples/2,modulation);
			for (unsigned i = 1;i<samples/2;i++){
				modphaseatindx(i,modulation);
				modphaseatindx(samples-i,modulation);
			}
			return 0;
		}
		inline int modulatephase_time(void) {
			if (infreq){
				std::cerr << "whoops, trying time modulation but in frequency domain" << std::endl;
				return 1;
			}
			modphaseatindx(0);
			modphaseatindx(samples/2);
			for (unsigned i = 1;i<samples/2;i++){
				modphaseatindx(i);
				modphaseatindx(samples-i);
			}
			return 0;
		}
		void print_amp(std::ofstream & outfile);
		void print_phase(std::ofstream & outfile);
		void print_phase_powerspectrum(std::ofstream & outfile);
		void printwavelengthbins(ofstream * outfile);
		void appendwavelength(ofstream * outfile);
		void appendfrequency(ofstream * outfile);
		void appendnoisy(ofstream * outfile);
		void appendfrequencybins(ofstream * outfile);
		void printfrequencybins(ofstream * outfile);
		void printfrequency(ofstream * outfile);
		void printfrequencydelay(ofstream * outfile, const double *delay);
		void printfrequencydelaychirp(ofstream * outfile, const double *delay,const double *chirp);
		void printtime(ofstream * outfile);
		void printwavelength(ofstream * outfile,const double *delay);
		inline double gettime(unsigned ind){return (time[ind]*fsPau<float>());}

		gsl_vector * modamp;
		gsl_vector * modphase;

	private:

		float m_noisescale;
		size_t m_sampleinterval;
		size_t m_lamsamples;
		unsigned long m_gain;
		unsigned m_saturate;

		const bool parent,child;
		bool intime,infreq;
		unsigned samples;
		unsigned startind,stopind,onwidth,offwidth;
		double tspan;
		double domega,lambda_center,lambda_width,omega_center,omega_width,omega_high;
		double omega_onwidth,omega_offwidth;
		double phase_GDD,phase_TOD,phase_4th,phase_5th;
		double dtime,time_center,time_wdith;
		gsl_fft_complex_wavetable * cwave;
		gsl_fft_complex_workspace * cspace;
		gsl_vector_complex * cvec;
		gsl_vector * rhovec;
		gsl_vector * phivec;
		double * omega, *time;
		double nu0;

		const unsigned i_low, i_high;


		void allocatetables(void);
		void buildvectors(void);
		void factorization(void);
		void killvectors(void);
		void killtheseonly(void);

		inline void addGDDtoindex(const unsigned indx,const int omega_sign) {
			static gsl_complex z;
			phivec->data[indx] += omega_sign*phase_GDD*gsl_pow_2(omega[indx]-(double(omega_sign)*omega_center));
			z = gsl_complex_polar(rhovec->data[indx],phivec->data[indx]);
			gsl_vector_complex_set(cvec,indx,z);
		}	
		inline void addTODtoindex(const unsigned indx,const int omega_sign) {
			static gsl_complex z;
			phivec->data[indx] += omega_sign*phase_TOD*gsl_pow_3(omega[indx]-(double(omega_sign)*omega_center));
			z = gsl_complex_polar(rhovec->data[indx],phivec->data[indx]);
			gsl_vector_complex_set(cvec,indx,z);
		}	
		inline void add4thtoindex(const unsigned indx,const int omega_sign) {
			static gsl_complex z;
			phivec->data[indx] += omega_sign*phase_4th*gsl_pow_4(omega[indx]-(double(omega_sign)*omega_center));
			z = gsl_complex_polar(rhovec->data[indx],phivec->data[indx]);
			gsl_vector_complex_set(cvec,indx,z);
		}	
		inline void add5thtoindex(const unsigned indx,const int omega_sign) {
			static gsl_complex z;
			phivec->data[indx] += omega_sign*phase_5th*gsl_pow_5(omega[indx]-(double(omega_sign)*omega_center));
			z = gsl_complex_polar(rhovec->data[indx],phivec->data[indx]);
			gsl_vector_complex_set(cvec,indx,z);
		}	

		void print_amp(std::ofstream & outfile);
		void print_phase(std::ofstream & outfile);
		void print_phase_powerspectrum(std::ofstream & outfile);
		void printwavelengthbins(std::ofstream * outfile);
		void appendwavelength(std::ofstream * outfile);
		void appendfrequency(std::ofstream * outfile);
		void appendnoisy(std::ofstream * outfile);
		void appendfrequencybins(std::ofstream * outfile);
		void printfrequencybins(std::ofstream * outfile);
		void printfrequency(std::ofstream * outfile);
		void printfrequencydelay(std::ofstream * outfile, const double *delay);
		void printfrequencydelaychirp(std::ofstream * outfile, const double *delay,const double *chirp);
		void printtime(std::ofstream * outfile);
		void printwavelength(std::ofstream * outfile,const double *delay);
		inline double gettime(unsigned ind){return (time[ind]*fsPau<double>());}


	private:

		double m_noisescale;
		size_t m_sampleinterval;
		size_t m_lamsamples;
		unsigned long m_gain;
		unsigned m_saturate;

		size_t sampleround;

		bool intime,infreq;
		unsigned samples;
		unsigned startind,stopind,onwidth,offwidth;
		double tspan;
		double domega,lambda_center,lambda_width,omega_center,omega_width,omega_high;
		double omega_onwidth,omega_offwidth;
		double phase_GDD,phase_TOD,phase_4th,phase_5th;
		double dtime,time_center,time_wdith;

		// FFTW variables //
		fftw_plan FTplan_forward;
		fftw_plan FTplan_backward;
		std::complex<double> * cvec; // this is still fftw_malloc() for sake of fftw memory alignment optimization

		std::vector<double> rhovec;
		std::vector<double> modamp;
		std::vector<double> omega;
		std::vector<double> time;
		std::vector<double> modphase;
		std::vector<double> phivec;

		double nu0;

		unsigned i_low, i_high;

		inline void rhophi2cvec(const size_t i){ cvec[i] = std::polar(rhovec[i],phivec[i]);}
		inline void cvec2rhophi(const size_t i){ rhovec[i] = std::abs(cvec[i]); phivec[i] = std::arg(cvec[i]); }
		void rhophi2cvec(void);
		void cvec2rhophi(void);

		std::random_device rng;

		void buildvectors(void);
		void factorization(void);
		void killvectors(void);
		void killtheseonly(void);

		inline void addGDDtoindex(const unsigned indx,const int omega_sign) {
			phivec[indx] += omega_sign*phase_GDD*std::pow(omega[indx]-(double(omega_sign)*omega_center),int(2));
			rhophi2cvec(indx);
		}	
		inline void addTODtoindex(const unsigned indx,const int omega_sign) {
			phivec[indx] += omega_sign*phase_TOD*std::pow(omega[indx]-(double(omega_sign)*omega_center),int(3));
			rhophi2cvec(indx);
		}	
		inline void add4thtoindex(const unsigned indx,const int omega_sign) {
			phivec[indx] += omega_sign*phase_4th*std::pow(omega[indx]-(double(omega_sign)*omega_center),int(4));
			rhophi2cvec(indx);
		}	
		inline void add5thtoindex(const unsigned indx,const int omega_sign) {
			phivec[indx] += omega_sign*phase_5th*std::pow(omega[indx]-(double(omega_sign)*omega_center),int(5));
			rhophi2cvec(indx);
		}	

		inline void modampatindx(const unsigned indx,const std::vector<double> & modvec) {
			rhovec[indx] *= modvec[indx];
			rhophi2cvec(indx);
		}
		inline void modampatindx(const unsigned indx) {
			rhovec[indx] *= modamp[indx];
			rhophi2cvec(indx);
		}
		/*
		   E(t) = E0exp(iwt)=E0 exp(i w(t) (t-t0))
		   = E0 exp(i (w0+GDD*t) (t-t0) )
		   = E0 exp(i ((w0 * t) - (GDD*t*t0) + (GDD*t**2) - (w0*t0)) so adding a delay means adding an arbitrary phase of w0*t0 and subtracting GDD*t*t0 ( a linear term in phase which we know only dials the pulse position in z, or likewise phase in time... unless t0 is actually t0(t), then you have to include this one as well as w0*t0(t).
		   = E0 exp(i ((w0 * t) + (GDD*t**2) - (w0*t0) - (GDD*t*t0) ) )
		   = E0 exp(i w(t)*t) exp(i -w(t)*t0(t))
		   = E(t) exp(i -(w0 + 1/GDD*t + 1/TOD*t**2 + 1/FOD*t**3)*t0(t))
		 */
		inline void modphaseatindx(const unsigned indx,const std::vector<double> & modvec) {
			if (modvec[indx]!=0){
				phivec[indx] += modvec[indx];
				rhophi2cvec(indx);
			}
		}
		inline void modphaseatindx(const unsigned indx) {
			if (modphase[indx] != 0){
				phivec[indx] += modphase[indx];
				rhophi2cvec(indx);
			}
		}
		inline double rising(const unsigned indx) {
			return std::pow(sin(double(Constants::half_pi<double>()*(indx-startind)/onwidth)),int(2));
		}
		inline double falling(const unsigned indx) {
			return std::pow(sin(double(Constants::half_pi<double>()*(stopind-indx)/offwidth)),int(2));
		}
};

#endif
