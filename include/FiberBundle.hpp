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
		bool print_mapping(std::ofstream & out);
		bool shuffle_output(void);
		void scalePolarCoords(void);
		inline void set_fsPmm(const float x = 0.5){fsPmm = x;}

		inline void fiberdiameter(const float x){ fiberdiam = x; for (size_t i=0;i<ovals.size();++i){ovals[i] = fiberdiam * float(i);}}
		inline void laserdiameter(const float x){ laserdiam = x; }
		inline void xraydiameter(const float x){ xraydiam = x; }
		inline float fiberdiameter(void){ return fiberdiam; }
		inline float laserdiameter(void){ return laserdiam; }
		inline float xraydiameter(void){ return xraydiam; }
		inline void Ilaser(const float x){ilaser=x;}
		inline void Ixray(const float x){ixray=x;}
		inline float Ilaser(void){return ilaser;}
		inline float Ixray(void){return ixray;}
		inline float Ilaser(const size_t i){return ilaser*std::exp(-std::pow(std::abs(zvals[ids[i]]-laser_center)/laserdiam,int(2)));}
		inline float Ixray(const size_t i){return ixray*std::exp(-std::pow(std::abs(zvals[ids[i]]-xray_center)/xraydiam,int(2)));}
		inline void center_Ilaser(const float dx,const float dy){laser_center = std::complex<float>(dx,dy);}
		inline void center_Ixray(const float dx,const float dy){xray_center = std::complex<float>(dx,dy);}
		inline std::complex<float> center_Ilaser(void){return laser_center;}
		inline std::complex<float> center_Ixray(void){return xray_center;}
		inline void delay_angle(const float a){alpha = a;}
		inline float delay_angle(void){return alpha;}
		inline float delay(const size_t i)
		{
			return fsPmm * std::abs(zvals[ids[i]])*std::cos(std::arg(zvals[ids[i]]) + alpha);
		}
		inline float x(const size_t i){return zvals[ids[i]].real();}
		inline float y(const size_t i){return zvals[ids[i]].imag();}
		inline float o(const size_t i){return ovals[ids[i]];}
		inline float r(const size_t i){return std::abs(zvals[ids[i]]);}
		inline float t(const size_t i){return std::arg(zvals[ids[i]]);}
		inline size_t get_nfibers(void){return nfibers;}
		inline size_t set_nfibers(const size_t n){nfibers = n; return nfibers;}

	private:
		size_t nfibers;
		size_t nrows;
		float fsPmm,fiberdiam,laserdiam,xraydiam;
		float ilaser,ixray,alpha;
		std::complex<float> laser_center;
		std::complex<float> xray_center;
		std::vector<size_t> ids;
		std::vector<float> ovals;
		std::vector<std::complex<float> > zvals;

		void setnrows(size_t n);
		bool set_polarcoords(size_t n);

	protected:
		float c,s;
};

#endif

