#include "FiberBundle.hpp"

#include "Constants.hpp"
using Constants::pi;
using Constants::pi_3;


FiberBundle::FiberBundle(size_t n = 109)
: ixray(1.0)
, ilaser(1.0)
, alpha(0.0)
, nfibers(n)
{
	//std::cerr << "In constructor FiberBundle() " << std::endl;

	if (!set_polarcoords(nfibers)) {
		std::cerr << "\n\n\t\t=========================" 
			<< "\n\t\t========================="
			<< "\n\t\t=========AAAAHHH========="
			<< "\n\t\t===alles ist scheisse===="
			<< "\n\t\t======check nfibers======"
			<< "\n\t\t========================="
			<< "\n\t\t========================="
			<< "\n\t\t========================="
			<< std::endl;
	}

        alpha = (double)atof(getenv("alpha"));

	ovals.resize(nfibers);
	zvals.resize(nfibers,std::complex<double>(0));
	ixray_vec.resize(nfibers);
	inds.resize(nfibers);
	DataOps::ramp(inds);
	fiberdiam = 0.11;
	laserdiam = 0.7;
	xraydiam = 0.5;
	thermaldiam = 0.5;
	xray_center = std::complex<double>(0.,0.);
	laser_center = std::complex<double>(0.,0.);
	thermalcenter = std::complex<double>(0.,0.);
	for (size_t i=0;i<ovals.size();++i){
		ovals[i] = fiberdiam * double(i);
	}
	fillIxray();
}

FiberBundle::FiberBundle(FiberBundle & rhs)
: ixray(rhs.ixray)
, ilaser(rhs.ilaser)
, alpha(rhs.alpha)
, nfibers(rhs.nfibers)
, fiberdiam(rhs.fiberdiam)
, laserdiam(rhs.laserdiam)
, xraydiam(rhs.xraydiam)
, thermaldiam(rhs.thermaldiam)
, xray_center(rhs.xray_center)
, laser_center(rhs.laser_center)
, thermalcenter(rhs.thermalcenter)
, fsPmm(rhs.fsPmm)
{
	ovals.resize(nfibers);
	zvals.resize(nfibers);
	ixray_vec.resize(nfibers);
	inds.resize(nfibers);
	std::copy(rhs.ovals.begin(),rhs.ovals.end(),ovals.begin());
	std::copy(rhs.zvals.begin(),rhs.zvals.end(),zvals.begin());
	std::copy(rhs.ixray_vec.begin(),rhs.ixray_vec.end(),ixray_vec.begin());
	std::copy(rhs.inds.begin(),rhs.inds.end(),inds.begin());
}

FiberBundle & FiberBundle::operator=(FiberBundle & rhs)
{
	ixray=rhs.ixray;
	ilaser=rhs.ilaser;
	alpha=rhs.alpha;
	nfibers=rhs.nfibers;
	fiberdiam=rhs.fiberdiam;
	laserdiam=rhs.laserdiam;
	xraydiam=rhs.xraydiam;
	thermaldiam=rhs.thermaldiam;
	xray_center=rhs.xray_center;
	laser_center=rhs.laser_center;
	thermalcenter = rhs.thermalcenter;
	fsPmm = rhs.fsPmm;
	if (ovals.size() != nfibers)
		ovals.resize(nfibers);
	if (zvals.size() != nfibers)
		zvals.resize(nfibers);
	if (ixray_vec.size() != nfibers)
		ixray_vec.resize(nfibers);
	if (inds.size() != nfibers)
		inds.resize(nfibers);
	std::copy(rhs.ovals.begin(),rhs.ovals.end(),ovals.begin());
	std::copy(rhs.zvals.begin(),rhs.zvals.end(),zvals.begin());
	std::copy(rhs.ixray_vec.begin(),rhs.ixray_vec.end(),ixray_vec.begin());
	std::copy(rhs.inds.begin(),rhs.inds.end(),inds.begin());
	return *this;
}

bool FiberBundle::set_inds(std::vector<uint16_t> in)
{
	if (in.size() != inds.size()){
		std::cerr << "setting FiberBundle inds with not the same length vector of uint8_t\n" << std::flush;
		return false;
	}
	std::copy(in.begin(),in.end(),inds.begin());
	return true;
}
bool FiberBundle::shuffle_inds(void)
{
	std::random_device rng;
	std::seed_seq seed{rng(), rng(), rng(), rng(), rng(), rng(), rng(), rng()};
	std::mt19937 e(seed);
	std::shuffle(inds.begin(),inds.end(),e);
	return true;
}

bool FiberBundle::shadow_xrays(const double xin, const double yin)
{
	fillIxray();
	for (size_t i = 0;i<ixray_vec.size(); ++i){ 
		if (zvals[i].real() > xin && zvals[i].imag() > yin)
			ixray_vec[i] *= 0.5;
	}
	return true;
}
double FiberBundle::fillIxray(void)
{
	if (ixray_vec.size() != zvals.size())
		ixray_vec.resize(zvals.size(),double(0.0));
	for (size_t i = 0;i<ixray_vec.size(); ++i){ 
		ixray_vec[i] = ixray * double(std::exp(-1.0*std::pow(std::abs(zvals[i]-xray_center)/xraydiam,int(2))));
	}
	return ixray;
}

bool FiberBundle::print_mapping(std::ofstream & out,double t0 = 0.)
{
	if (!out.is_open() )
		return false;
	out << "#delay = " << t0 << "\n"; 
	out << "#key\tr\ttheta\tx\ty\to\tdelay\tIlas\tIxray\tTinK\n"; 
	for (size_t i=0;i<zvals.size();++i){
			out << inds[i] << "\t"
			<< std::abs(zvals[i]) << "\t" 
			<< std::arg(zvals[i]) << "\t" 
			<< zvals[i].real() << "\t" 
			<< zvals[i].imag() << "\t" 
			<< ovals[i] << "\t" 
			<< (t0+delay(i)) << "\t" 
			<< Ilaser(i) << "\t" 
			<< ixray_vec[i] << "\t"
			<< TinK(i) << "\n";
	}
	out << std::flush;
	return true;
}


bool FiberBundle::set_polarcoords(size_t n)
{
	/*
	   if (nfibers!=109) // for now enforce 109 fiber bundle
	   nfibers = 109;
	 */
	if (n < 1) {
		std::cerr << "nfibers must be larger than 1, and follow (6*r + 1) for 6 fold symmetry" << std::endl;
		return false;
	}
	c=std::cos(pi_3<double>());
	s=std::sin(pi_3<double>());

	setnfibers(n);
	std::cout << "\n========== running " << nfibers << " fibers ============\n" << std::flush;
	zvals.resize(nfibers);

	std::complex<double> z(0.0,0.0);
	zvals[0] = std::complex<double>(0.0,0.0);
	std::fill(zvals.begin(),zvals.end(),std::complex<double>(0.,0.));
	double theta;
	for (size_t i=0;i<6;++i){
		theta = double(i)*pi_3<double>();
		if (nfibers < 2){continue;}
			zvals[1+i] = std::polar(1.,theta);
		if (nfibers < 8){continue;}
			zvals[6+1+i] = std::polar(2.,theta);
			z = std::complex<double>(1.+c,s);
			zvals[2*6+i+1] = std::polar(std::abs(z),double(2*i+1)*std::arg(z));
		if (nfibers < 20){continue;}
			zvals[3*6+i+1] = std::polar(3.0, theta);
			z = std::complex<double>(2+c,s);
			zvals[4*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
			zvals[5*6+i+1] = std::polar(std::abs(z), theta - std::arg(z));
		if (nfibers < 38){continue;}
			zvals[6*6+i+1] = std::polar(4.0, theta);
			z = std::complex<double>(3+c,s);
			zvals[7*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
			zvals[8*6+i+1] = std::polar(std::abs(z), theta - std::arg(z));
			z = std::complex<double>(2+2*c,2*s);
			zvals[9*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
		if (nfibers <62){ continue;}
			zvals[10*6+i+1] = std::polar(5.0, theta);
			z = std::complex<double>(4+c,s);
			zvals[11*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
			zvals[12*6+i+1] = std::polar(std::abs(z), theta - std::arg(z));
			z = std::complex<double>(3+c,3*s);
			zvals[13*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
			zvals[14*6+i+1] = std::polar(std::abs(z), theta - std::arg(z));
			z = std::complex<double>(4+c,3*s);
			zvals[15*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
			z = std::complex<double>(4+2*c,2*s);
			zvals[16*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
			zvals[17*6+i+1] = std::polar(std::abs(z), theta - std::arg(z));
		if (nfibers < 110){ continue; }
	}
	return true;
}

void FiberBundle::scalePolarCoords(void)
{
	std::transform(zvals.begin(), zvals.end(), zvals.begin(), std::bind2nd(std::multiplies< std::complex<double> >(),fiberdiam));
}
void FiberBundle::print_zvals(void)
{
	for (size_t i = 0 ; i< zvals.size(); ++i){
		std::cerr << zvals[i].real() << "," << zvals[i].imag() << std::endl;
	}
	std::cerr << std::flush;
}

void FiberBundle::setnfibers(size_t n)
{
	if (n<2){ // r=0
		nfibers = 1;
		return; // this works... don't know why
	}
	if (n<8){ // r=1
		nfibers = 7;
		return; // this is failing... don't know why
	}
	if (n<20){ // r=2
		nfibers = 19;
		return; // this works... don't know why
	}
	if (n<38){ // r=3
		nfibers = 37;
		return;
	}
	if (n<62){ // r=4
		nfibers = 61;
		return;
	}
	if (n<110){ // r=5
		nfibers = 109;
		return;
	}
}

