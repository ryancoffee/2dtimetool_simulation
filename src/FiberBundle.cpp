#include "FiberBundle.hpp"

FiberBundle::FiberBundle(size_t n = 109)
: nfibers(n)
, ixray(1.0)
, ilaser(1.0)
, alpha(0.0)
{
	ids.resize(nfibers);
	ovals.resize(nfibers);
	zvals.resize(nfibers,std::complex<float>(0));
	c=std::cos(float(M_PI)/3);
	s=std::sin(float(M_PI)/3);
	fiberdiam = 0.11;
	laserdiam = 0.7;
	xraydiam = 0.5;
	xray_center = std::complex<float>(0.,0.);
	laser_center = std::complex<float>(0.,0.);
	for (size_t i=0;i<ovals.size();++i){
		ovals[i] = fiberdiam * float(i);
		ids[i] = i;
	}
}

bool FiberBundle::shuffle_output(void)
{
	std::random_device rng;
	std::seed_seq seed{rng(), rng(), rng(), rng(), rng(), rng(), rng(), rng()};
	std::mt19937 e(seed);
	std::shuffle(ids.begin(),ids.end(),e);
}

bool FiberBundle::print_mapping(std::ofstream & out)
{
	if (!out.is_open() )
		return false;
	out << "#r\ttheta\tx\ty\to\tdelay\tIlas\tIxray\n"; 
	for (size_t i=0;i<zvals.size();++i){
		out << std::abs(zvals[i]) << "\t" 
			<< std::arg(zvals[ids[i]]) << "\t" 
			<< zvals[ids[i]].real() << "\t" 
			<< zvals[ids[i]].imag() << "\t" 
			<< ovals[ids[i]] << "\t" 
			<< delay(ids[i]) << "\t" 
			<< Ilaser(ids[i]) << "\t" 
			<< Ixray(ids[i]) << "\n";
	}
	out << std::flush;
	return true;
}


void FiberBundle::set_polarcoords(void)
{
	std::complex<float> z(0.0f,0.0f);
	if (nfibers!=109) // for now enforce 109 fiber bundle
		nfibers = 109;
		zvals.resize(nfibers);
	zvals[0] = std::complex<float>(0.f,0.f);
	float theta;
	for (size_t i;i<6;++i){
		theta = float(i)*M_PI/3.0f;
		zvals[0*6+i+1] = std::polar(1.0f, theta); 
		zvals[1*6+i+1] = std::polar(2.0f, theta); 
		z = std::complex<float>(1+c,s);
		zvals[2*6+i+1] = std::polar(std::abs(z),float(2*i+1)*std::arg(z));
		zvals[3*6+i+1] = std::polar(3.0f, theta);
		z = std::complex<float>(2+c,s);
		zvals[4*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
		zvals[5*6+i+1] = std::polar(std::abs(z), theta - std::arg(z));
		zvals[6*6+i+1] = std::polar(4.0f, theta);
		z = std::complex<float>(3+c,s);
		zvals[7*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
		zvals[8*6+i+1] = std::polar(std::abs(z), theta - std::arg(z));
		z = std::complex<float>(2+2*c,2*s);
		zvals[9*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
		zvals[10*6+i+1] = std::polar(5.0f, theta);
		z = std::complex<float>(4+c,s);
		zvals[11*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
		zvals[12*6+i+1] = std::polar(std::abs(z), theta - std::arg(z));
		z = std::complex<float>(3+c,3*s);
		zvals[13*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
		zvals[14*6+i+1] = std::polar(std::abs(z), theta - std::arg(z));
		z = std::complex<float>(4+c,3*s);
		zvals[15*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
		z = std::complex<float>(4+2*c,2*s);
		zvals[16*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
		zvals[17*6+i+1] = std::polar(std::abs(z), theta - std::arg(z));
	}
	//zvals *= fiberdiam; // Don't know why this isn't handeled in DataOps.hpp
	std::transform(zvals.begin(), zvals.end(), zvals.begin(), std::bind2nd(std::multiplies< std::complex<float> >(),fiberdiam));
}

