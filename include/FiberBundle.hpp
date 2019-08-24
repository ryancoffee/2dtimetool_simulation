#ifndef FIBERBUNDLE

#define FIBERBUNDLE

#define _USE_MATH_DEFINES
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <complex>
#include <iostream>
#include <fstream>

#include "DataOps.hpp"

class FiberBundle {
	public:
		FiberBundle(size_t n);
		bool print_mapping(std::ofstream & out,double t0);
		bool shuffle_output(void);
		void scalePolarCoords(void);
		inline void set_fsPmm(const double x = 3333){fsPmm = x;}

		inline void fiberdiameter(const double x){ fiberdiam = x; for (size_t i=0;i<ovals.size();++i){ovals[i] = fiberdiam * double(i);} }
		inline double laserdiameter(const double x){ laserdiam = x; return laserdiam; }
		inline double xraydiameter(const double x){ xraydiam = x; return xraydiam; }
		inline double thermaldiameter(const double x){ thermaldiam = x; return thermaldiam; }
		inline double fiberdiameter(void){ return fiberdiam; }
		inline double laserdiameter(void){ return laserdiam; }
		inline double xraydiameter(void){ return xraydiam; }
		inline double thermaldiameter(void){ return thermaldiam; }
		inline double Ilaser(const double x){ilaser=x; return ilaser; }
		inline double Ixray(const double x){ixray=x; return ixray; }
		inline double Ilaser(void){return ilaser;}
		inline double Ixray(void){return ixray;}
		inline double Ilaser(const size_t i){return ilaser*std::exp(-std::pow(std::abs(zvals[ids[i]]-laser_center)/laserdiam,int(2)));}
		void setTmax_Tbase(const double Tmax,const double Tbase){Tmax_K = Tmax; Tbase_K = Tbase; }
		inline double TinK(const size_t i){Tmax_K*std::exp(-std::pow(std::abs(zvals[ids[i]]-thermalcenter)/thermaldiam,int(2)));}

		void print_zvals(void);

		inline double Ixray(const size_t i){return std::exp(-1.0*std::pow(std::abs(zvals[ids[i]]-xray_center)/xraydiam,int(2)));}
		//inline double Ixray(const size_t i){return ixray*std::exp(-1.0*std::pow(std::abs(zvals[ids[i]]-xray_center)/xraydiam,int(2)));}
		inline std::complex<double> center_Ilaser(const double dx,const double dy){laser_center = std::complex<double>(dx,dy); return laser_center; }
		inline std::complex<double> center_Ixray(const double dx,const double dy){xray_center = std::complex<double>(dx,dy); return xray_center; }
		inline std::complex<double> center_thermal(const double dx,const double dy){thermalcenter = std::complex<double>(dx,dy); return thermalcenter; }
		inline std::complex<double> center_Ilaser(void){return laser_center;}
		inline std::complex<double> center_Ixray(void){return xray_center;}
		inline std::complex<double> center_thermal(void){return thermalcenter;}
		inline double delay_angle(const double a){alpha = a; return alpha; }
		inline double delay_angle(void){return alpha;}
		inline double delay(const size_t i){ return fsPmm * (std::cos(alpha)*x(i) + std::sin(alpha)*y(i)); }
		inline double x(const size_t i){return zvals[ids[i]].real();}
		inline double y(const size_t i){return zvals[ids[i]].imag();}
		inline double o(const size_t i){return ovals[ids[i]];}
		inline double r(const size_t i){return std::abs(zvals[ids[i]]);}
		inline double t(const size_t i){return std::arg(zvals[ids[i]]);}
		inline size_t get_nfibers(void){return nfibers;}
		inline size_t set_nfibers(const size_t n){nfibers = n; return nfibers;}

	private:
		size_t nfibers;
		double fsPmm,fiberdiam,laserdiam,xraydiam,thermaldiam;
		double ilaser,ixray,alpha,Tmax_K,Tbase_K;
		std::complex<double> laser_center;
		std::complex<double> xray_center;
		std::complex<double> thermalcenter;
		std::vector<size_t> ids;
		std::vector<double> ovals;
		std::vector<std::complex<double> > zvals;

		void setnfibers(size_t n);
		bool set_polarcoords(size_t n);

	protected:
		double c,s;
};

#endif

