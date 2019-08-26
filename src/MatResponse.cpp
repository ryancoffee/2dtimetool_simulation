#include <MatResponse.hpp>
#include <Constants.hpp>
#include <DataOps.hpp>
#include <exception>
#include <algorithm>
#include <boost/math/interpolators/barycentric_rational.hpp>
#include <limits>

MatResponse::MatResponse(MatResponse & rhs) // copy constructor
: carriers_set(rhs.carriers_set)
, t0(rhs.t0)
, twidth(rhs.twidth)
, scale(rhs.scale)
, attenuation(rhs.attenuation)
, phase(rhs.phase)
, a(rhs.a)
, b(rhs.b)
, alpha(rhs.alpha)
, beta(rhs.beta)
, etalondelay(rhs.etalondelay)
, reflectance(rhs.reflectance)
, bandgap_eV(rhs.bandgap_eV)
, n_0(rhs.n_0)
, n_a(rhs.n_a)
, n_b(rhs.n_b)
, l_thickness_um(rhs.l_thickness_um)
{
	carriers.resize(rhs.carriers.size(),0.);
	std::copy(rhs.carriers.begin(),rhs.carriers.end(),carriers.begin());
}        

MatResponse::MatResponse(double t0_in=0.0,double width_in=10.0,double atten_in = 0.05,double phase_in = 0.03)
: t0(t0_in / fsPau<double>())
, twidth(width_in * root_pi<double>()/ fsPau<double>() / 2.0)
, attenuation(atten_in)
, phase(phase_in)
, carriers_set(false)
, n_0(2.38)
, n_a(1.5)
, n_b(-13.9931)
, l_thickness_um(5.5)
{     
	carriers.resize(10,0.);
	std::cerr << "In constructor MatRespons() " << std::endl;
	a=0.75;
	b=0.25;
	alpha=5.0e-4*fsPau<double>();
	beta=0.0;
	scale=1.0;
	bandgap_eV = 5.47; // this is 5.47 for diamond.
}        

MatResponse & MatResponse::operator=(const MatResponse & rhs)
{
	carriers_set = rhs.carriers_set;
	t0 = rhs.t0;
	twidth = rhs.twidth;
	scale = rhs.scale;
	attenuation = rhs.attenuation;
 	phase = rhs.phase;
	a = rhs.a;
	b = rhs.b;
	alpha = rhs.alpha;
	beta = rhs.beta;
	etalondelay = rhs.etalondelay;
	reflectance = rhs.reflectance;
	bandgap_eV = rhs.bandgap_eV;
	n_0 = rhs.n_0;
	n_a = rhs.n_a;
	n_b = rhs.n_b;
	l_thickness_um = rhs.l_thickness_um;
	carriers.resize(rhs.carriers.size(),0.);
	std::copy(rhs.carriers.begin(),rhs.carriers.end(),carriers.begin());
	return *this;
}

/*
 * total(x)=f(x)*fall(x)+e(x)*rise(x)+y0
 * f(x)=a*x**2+c*x**p
 * fall(x)=0.5*(1-erf((x-xfall)/wfall)
 * e(x)=aa*exp(-(x-xfall)/ww)
 * rise(x)=0.5*(1+erf((x-xfall)/wfall))
 *
#f(x)=a*x**2+c*x**p
#e(x)=aa*exp(-(x-xfall)/ww)
#fall(x)=0.5*(1-erf((x-xfall)/wfall))
#rise(x)=0.5*(1+erf((x-xfall)/wfall))
#y0=1e-4
#total(x) = f(x)*fall(x)+e(x)*rise(x)+y0
 * def diamond_cascade_params_improved(x,energy=9500):
    # lookout, energy needs to be in eV not keV as before
    #still basically hand fit, but these are more complete including long tail as additional exponential
    # all params seemed to fit a nice power law
    #f(x)=a*x**2+c*x**p
    #e(x)=aa*exp(-(x-xfall)/ww)
    #fall(x)=0.5*(1-erf((x-xfall)/wfall))
    #rise(x)=0.5*(1+erf((x-xfall)/wfall))
    #total(x) = f(x)*fall(x)+e(x)*rise(x)+y0
    ##
    #
    y0=1e-4
    a = 7.e7 * np.power(float(energy),float(-2.5))
    c = 128. * np.power(energy-531,-1./3)
    p = 10 * np.power(energy,-0.5)
    aa = 0.499467* np.exp(-np.power(energy/45.e3,int(2)))
    ww = 10*np.power(energy,2.e-6*energy)
    xfall = 8.24e-05 * np.power(energy,4./3)
    wfall = 0.000373 * energy
    y=np.zeros(x.shape)
    inds = np.where(x>0)
    f = a*np.power(x[inds],int(2))+c*np.power(x[inds],p)
    fall = 0.5*(1-erf((x[inds] - xfall)/wfall))
    e=aa*np.exp(-(x[inds]-xfall)/ww)
    rise=0.5*(1+erf((x[inds]-xfall)/wfall))
    y[inds] = f*fall+e*rise+y0
    return y

*/

void MatResponse::set_n_refractive(const double n_0in = 2.38, const double n_ain = 1.5, const double n_bin = -13.9931){
	n_0 = n_0in;
	n_a = n_ain;
	n_b = n_bin;
}
void MatResponse::set_thickness(const double lin){
	l_thickness_um = lin;
}

bool MatResponse::fill_carriersvec(PulseFreq * pulse,std::ifstream * instream){return fill_carriersvec(*pulse,*instream);}
bool MatResponse::fill_carriersvec(PulseFreq & pulse,std::ifstream & instream)
{
	using namespace DataOps;
	size_t column = (size_t)atoi(getenv("carriercolumn"));
	carriers.resize(pulse.getsamples(),double(0.));
	std::vector<double> times(carriers.size(),double(0.));
	std::vector< std::vector<double> > data;
	instream >> data;
	std::cerr << "data.size() = " << data.size() << "\t" << data.front().size() << "\n" << std::flush;
	if(column>=data.front().size()){
		std::cerr << "failed data column in fill_carrersvec() method\n" << std::flush;
		return fill_carriersvec(pulse,9.5);
	}
	double tstep = data[1][0] - data[0][0];
	std::vector <double> c_times(data.size()+2);
	std::vector <double> c_vholes(data.size()+2,double(0));
	for (size_t i=0;i<data.size();++i){
		if (i<10 || i>data.size()-10){
			std::cerr << data[i][0] << "\t" << data[i][column] << "\t...\t" << std::flush;
		}
		c_times[i+2] = data[i][0];
		c_vholes[i+2] = data[i][column];
		if (i<10 || i>data.size()-10){
			std::cerr << c_times[i] << "\t" << c_vholes[i] << "\n" << std::flush;
		}
	}
	c_times[0] = (-double(carriers.size())*pulse.getdt()*Constants::fsPau<double>());
	double smallestdouble = std::numeric_limits<double>::min();
	c_times[1] = -1.*(smallestdouble);
	c_vholes[0] = 0.;
	c_vholes[1] = 0.;
	double ncarriers_final = 1;// double(energy_keV*1e3/(3*bandgap_eV)); // this is a rule of thumb for exciton energy
	double datamax = *std::max_element(c_vholes.begin(),c_vholes.end());
	std::transform(c_vholes.begin(), c_vholes.end(), c_vholes.begin(), std::bind2nd(std::multiplies<double>(),ncarriers_final/datamax) );

	std::cerr << "\n\t\t===== ready to fit the input data for carriersvec ==========\n" << std::flush;
	for (size_t i = 0; i< 10;++i){
		std::cerr << "\n\t\t===== " << c_times[i] << "\t" << c_vholes[i] << " ==========\n" << std::flush;
	}

	//boost::math::barycentric_rational<double> interpolant(c_times.data(), c_vholes.data(), c_vholes.size());
	for (size_t i = 0 ; i < carriers.size(); i++){
		times[i] = double(i)*pulse.getdt()*Constants::fsPau<double>(); // remember, this needs to be in femtoseconds... and getdt() returns in atomic units
		carriers[i] = interpolate(c_times,c_vholes,times[i]);
	}

	std::vector<double> shortcopy(40);
	std::vector<double> shortcopy_times(40);
	sample_every(shortcopy,carriers,100);
	sample_every(shortcopy_times,times,100);
	std::cerr << shortcopy_times << "\n" << std::flush;
	std::cerr << shortcopy << "\n" << std::flush;

	carriers_set = true;
	return carriers_set;
}

bool MatResponse::fill_carriersvec(PulseFreq * pulse,double energy_keV = 9.5){return fill_carriersvec(*pulse,energy_keV);}
bool MatResponse::fill_carriersvec(PulseFreq & pulse,double energy_keV = 9.5){
	using namespace DataOps;
	double e = energy_keV*1e3;
	// ultimately, these were had fit with Nikita's 2015 derivative curves... those are too slow by 50% or so he says
	carriers.resize(pulse.getsamples(),double(0.));
	std::vector<double> times(carriers.size(),double(0.));
	std::vector<double> decay(carriers.size(),double(0.));
	//std::cerr << "pulse.getsamples() = " << pulse.getsamples() << "\t" << std::flush;
	//std::cerr << " pulse.modamp.size() = " << pulse.modamp.size() << "\n" << std::flush;
	double a = 7.e7 * std::pow(e,float(-2.5));
	double c = 128. * std::pow(e,float(-1./3));
	double p = 10. * std::pow(e,float(-0.5)); 
	double aa = 0.5 * std::exp(-std::pow(energy_keV/45.,int(2)));
	double ww = 10.*std::pow(e,e*2e-6);
	double xfall = 8.24e-05 * std::pow(e,float(4./3.));
	double wfall = 0.000373 * e;
	double y0 = 1e-4;

	carriers[0] = 1.;
	decay[0] = 0.;
	for (size_t i = 1 ; i < carriers.size()/2; i++){
		double arg = double(i)*pulse.getdt()*Constants::fsPau<double>(); // remember, this needs to be in femtoseconds... and getdt() returns in atomic units
		double f = a*std::pow(arg,int(2)) + c * std::pow(arg,p);
		double fall = 0.5*(1.-erf((arg-xfall)/wfall));
		double e = aa*std::exp(-(arg-xfall)/ww);
		double rise = 0.5 * (1.+erf((arg-xfall)/wfall));
		times[i] = arg;
		carriers[i] = carriers[i-1] + (y0 + f*fall + e*rise);
		decay[i] = decay[i-1] + 0.5*erf((arg-xfall)/wfall);
	}
	size_t ncarriers_final = 1;// size_t(energy_keV*1e3/(3*bandgap_eV)); // this is a rule of thumb for exciton energy
	double scale;
	scale = ncarriers_final / *std::max_element(carriers.begin(),carriers.end());
	std::transform(carriers.begin(), carriers.end(), carriers.begin(), std::bind2nd(std::multiplies<double>(),scale) );
	scale = ncarriers_final / *std::max_element(decay.begin(),decay.end());
	std::transform(decay.begin(), decay.end(), decay.begin(), std::bind2nd(std::multiplies<double>(),scale) );
	carriers -= decay;
	std::vector<double> shortcopy(40);
	std::vector<double> shortcopy_times(40);
	sample_every(shortcopy,carriers,100);
	sample_every(shortcopy_times,times,100);
	std::cerr << shortcopy_times << "\n" << std::flush;
	std::cerr << shortcopy << "\n" << std::flush;
	carriers_set = true;
	return carriers_set;
}

bool MatResponse::setstepvec_both_carriers(PulseFreq * pulse,double delay_in){ return setstepvec_both_carriers(*pulse,delay_in); }
bool MatResponse::setstepvec_both_carriers(PulseFreq & pulse,double delay_in)
{ 
	int arg(0);
	//std::cerr << "HERE HERE in setstepvec_both_carriers() method\n" << std::flush;
	std::fill(pulse.modamp.begin(),pulse.modamp.end(),1.);
	std::fill(pulse.modphase.begin(),pulse.modphase.end(),0.);
	return addstepvec_both_carriers(pulse,delay_in);
}
bool MatResponse::addstepvec_both_carriers(PulseFreq * pulse,double delay_in){ return addstepvec_both_carriers(*pulse,delay_in); }
bool MatResponse::addstepvec_both_carriers(PulseFreq & pulse,double delay_in)
{
	int arg(0);
	double thisamp(0);
	double thisphase(0);
	CarriersNotSet carriers_not_set;

	try {
		if (!carriers_set){throw carriers_not_set;}

	//std::cerr << "HERE  in addstepvec_both_carriers() method\n" 
		//" t0 = " << t0 << " delay_in = " << delay_in << "carriers.size() = " << carriers.size() << "\n" << "i,arg = " << std::flush;
		for (size_t i = 0;i<pulse.getsamples()/2;++i){
			thisamp = 0.;
			thisphase = 0.;
			arg = std::min(int(i) - int((t0-delay_in)/pulse.getdt()) , int(carriers.size()) - 1);
	//std::cerr << i << "," << arg << " " << std::flush;
			if (arg>0){
				thisamp = -attenuation * carriers[arg]; 
				thisphase = phase * carriers[arg];
			}
			pulse.modamp[i] *= (1.+ thisamp);
			pulse.modphase[i] += thisphase;
		}
	//std::cerr << "HERE again in addstepvec_both_carriers() method\n" << std::flush;
		for (size_t i = pulse.getsamples()/2; i< pulse.getsamples() ;++i){
			thisamp = 0.;
			thisphase = 0.;
			arg = std::min(int(i)-int(pulse.getsamples()) - int((t0-delay_in)/pulse.getdt()) , int(carriers.size()) - 1);
	//std::cerr << i << "," << arg << " " << std::flush;
			if (arg>0){
				thisamp = -attenuation * carriers[arg];
				thisphase = phase * carriers[arg];
			}
			pulse.modamp[i] *= (1.+ thisamp);
			pulse.modphase[i] += thisphase;
		}

	} catch(exception& e) {
		std::cerr << "Error in setstepvec_both_carriers() method: " << e.what() << std::endl << std::flush;
		addstepvec_amp(pulse,delay_in);
		addstepvec_phase(pulse,delay_in);
		return false;
	}
	return true;
}

void MatResponse::setstepvec_amp(PulseFreq * pulse,double delay_in){ setstepvec_amp(*pulse,delay_in); }
void MatResponse::setstepvec_amp(PulseFreq & pulse,double delay_in){
	double arg(0);
	for (unsigned i = 0 ; i < pulse.getsamples(); i++){
		pulse.modamp[i] = 0.;
		arg = (static_cast<double>(i)*pulse.getdt() - (t0-delay_in));
		if( i > pulse.getsamples()/2){
			arg -= (pulse.getdt()*pulse.getsamples());
		}
		if (arg > 0.0) {
			pulse.modamp[i] = -attenuation * (a*std::exp(-1.0*(arg*alpha)) + b*std::exp(-1.0*(arg*beta)));
			pulse.modamp[i] *= scale;
			if (arg < twidth) {
				pulse.modamp[i] *= std::pow( sin(M_PI_2*arg/twidth ) , int(2));
			}
		}
		pulse.modamp[i] += 1.0;
	}
}

void MatResponse::setstepvec_phase(PulseFreq * pulse, double delay_in){ setstepvec_phase(*pulse,delay_in); }
void MatResponse::setstepvec_phase(PulseFreq & pulse,double delay_in){
	double arg(0);
	for (unsigned i = 0 ; i < pulse.getsamples(); i++){
		pulse.modphase[i] = 0.0;
		arg = (static_cast<double>(i)*pulse.getdt() - (t0-delay_in));
		if( i > pulse.getsamples()/2){
			arg -= (pulse.getdt()*pulse.getsamples());
		}
		if (arg > 0.0) {
			pulse.modphase[i] = phase * (a*std::exp(-1.0*(arg*alpha)) + b*std::exp(-1.0*(arg*beta)));
			pulse.modphase[i] *= scale;
			if( arg < twidth){
				pulse.modphase[i] *= std::pow( sin(M_PI_2*arg/twidth ) ,int(2) ) ;
			}
		}
	}
}

void MatResponse::buffervectors(PulseFreq * pulse){ buffervectors(*pulse); }
void MatResponse::buffervectors(PulseFreq & pulse){
	// nominally setting this to 10 times the fastest twidth to minimize its effect on spectral broadening
	unsigned bufferwidth = static_cast<unsigned>(10.0*twidth/pulse.getdt());
	double bufferscale;
	double arg(0);

	pulse.modphase[pulse.getsamples()/2] *= 0.0;
	pulse.modamp[pulse.getsamples()/2] = 1.0;
	for (unsigned i=0; i < bufferwidth ; i++) {
		arg = ((double)i)/((double)bufferwidth);
		bufferscale = std::pow( sin(M_PI_2*(arg) ) , int(2)) ; // buffering the vector back to 0.0 so it can be periodic.
		pulse.modphase[pulse.getsamples()/2-i] *= bufferscale;
		pulse.modamp[pulse.getsamples()/2-i] -= 1.0;
		pulse.modamp[pulse.getsamples()/2-i] *= bufferscale;
		pulse.modamp[pulse.getsamples()/2-i] += 1.0;
	}
}


void MatResponse::addstepvec_amp(PulseFreq * pulse,double delay_in){ addstepvec_amp(*pulse,delay_in); }
void MatResponse::addstepvec_amp(PulseFreq & pulse,double delay_in){
	double arg(0);
	double thisamp(0);
	for (unsigned i = 0 ; i < pulse.getsamples(); i++){
		thisamp = 0.;
		arg = (static_cast<double>(i)*pulse.getdt() - (t0-delay_in));
		if( i > pulse.getsamples()/2){
			arg -= (pulse.getdt()*pulse.getsamples());
		}
		if (arg > 0.0) {
			thisamp = -attenuation * (a*exp(-1.0*(arg*alpha)) + b*exp(-1.0*(arg*beta)));
			thisamp *= scale;
			if (arg < twidth) {
				thisamp *= std::pow( sin(M_PI_2*arg/twidth ) , int(2) ) ;
			}
		}
		pulse.modamp[i] *= 1.0 + thisamp;
	}

}


void MatResponse::addstepvec_phase(PulseFreq * pulse,double delay_in){ addstepvec_phase(*pulse,delay_in); }
void MatResponse::addstepvec_phase(PulseFreq & pulse,double delay_in){
	double arg(0.);
	double thisphase(0.);
	for (unsigned i = 0 ; i < pulse.getsamples(); i++){
		arg = (static_cast<double>(i)*pulse.getdt() - (t0-delay_in));
		thisphase = 0.;
		if( i > pulse.getsamples()/2){
			arg -= (pulse.getdt()*pulse.getsamples());
		}
		if (arg > 0.0) {
			thisphase = phase * (a*exp(-1.0*(arg*alpha)) + b*exp(-1.0*(arg*beta)));
			thisphase *= scale;
			if( arg < twidth){
				thisphase *= std::pow( sin(M_PI_2*arg/twidth ) , int(2)) ;
			}
		}
		pulse.modphase[i] += thisphase;
	}

}

