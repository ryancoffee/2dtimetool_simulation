#include "FiberBundle.hpp"

FiberBundle::FiberBundle(size_t n = 109)
: ixray(1.0)
, ilaser(1.0)
, alpha(0.0)
, nrows(1)
{

	if (!set_polarcoords(n)) {
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

	setnrows(n);
	std::cout << "\n========== running " << nfibers << " fibers ============\n\n" << std::flush;
	zvals.resize(nfibers);

	std::complex<float> z(0.0f,0.0f);
	zvals[0] = std::complex<float>(0.f,0.f);
	float theta;
	for (size_t i=0;i<6;++i){
		theta = float(i)*M_PI/3.0f;
		for (size_t r = 0;r<nrows;++r){
			if (nfibers < 2) continue; 
			zvals[r*6+i+1] = std::polar(1.0f, theta); 
			if (nfibers < 8) continue; 
			zvals[r*6+i+1] = std::polar(2.0f, theta); 
			z = std::complex<float>(1+c,s);
			zvals[r*6+i+1] = std::polar(std::abs(z),float(2*i+1)*std::arg(z));
			if (nfibers < 20) continue; 
			zvals[r*6+i+1] = std::polar(3.0f, theta);
			z = std::complex<float>(2+c,s);
			zvals[r*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
			zvals[r*6+i+1] = std::polar(std::abs(z), theta - std::arg(z));
			if (nfibers < 38) continue; 
			zvals[r*6+i+1] = std::polar(4.0f, theta);
			z = std::complex<float>(3+c,s);
			zvals[r*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
			zvals[r*6+i+1] = std::polar(std::abs(z), theta - std::arg(z));
			z = std::complex<float>(2+2*c,2*s);
			zvals[r*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
			if (nfibers < 62) continue; 
			zvals[r*6+i+1] = std::polar(5.0f, theta);
			z = std::complex<float>(4+c,s);
			zvals[r*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
			zvals[r*6+i+1] = std::polar(std::abs(z), theta - std::arg(z));
			z = std::complex<float>(3+c,3*s);
			zvals[r*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
			zvals[r*6+i+1] = std::polar(std::abs(z), theta - std::arg(z));
			z = std::complex<float>(4+c,3*s);
			zvals[r*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
			z = std::complex<float>(4+2*c,2*s);
			zvals[r*6+i+1] = std::polar(std::abs(z), theta + std::arg(z));
			zvals[r*6+i+1] = std::polar(std::abs(z), theta - std::arg(z));
			if (nfibers < 110) continue; 
		}
	}
	return true;
}

void FiberBundle::scalePolarCoords(void)
{
	std::transform(zvals.begin(), zvals.end(), zvals.begin(), std::bind2nd(std::multiplies< std::complex<float> >(),fiberdiam));
}

void FiberBundle::setnrows(size_t n)
{
	if (n<2){ // r=0
		nfibers = 1;
		nrows = 0;
		return; // this works... don't know why
	}
	if (n<8){ // r=1
		nfibers = 7;
		nrows = 1;
		return; // this is failing... don't know why
	}
	if (n<20){ // r=2
		nfibers = 19;
		nrows = 3;
		return; // this works... don't know why
	}
	if (n<38){ // r=3
		nfibers = 37;
		nrows = 6;
		return;
	}
	if (n<62){ // r=4
		nfibers = 61;
		nrows = 10;
		return;
	}
	if (n<110){ // r=5
		nfibers = 109;
		nrows = 18;
		return;
	}
}

