#include <MatResponse.hpp>
#include <Constants.hpp>
#include <DataOps.hpp>
#include <exception>
#include <algorithm>

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
{
		carriers.resize(rhs.carriers.size(),0.);
		std::copy(rhs.carriers.begin(),rhs.carriers.end(),carriers.begin());
}        
hello my 

MatResponse::MatResponse(double t0_in=0.0,double width_in=10.0,double atten_in = 0.05,double phase_in = 0.03)
: t0(t0_in / fsPau<double>())
, twidth(width_in * root_pi<double>()/ fsPau<double>() / 2.0)
, attenuation(atten_in)
, phase(phase_in)
, carriers_set(false)
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

bool MatResponse::fill_carriersvec(PulseFreq * pulse,double energy_keV = 9.5){return fill_carriersvec(*pulse,energy_keV);}
bool MatResponse::fill_carriersvec(PulseFreq & pulse,double energy_keV = 9.5){
	using namespace DataOps;
	// ultimately, these were had fit with Nikita's 2015 derivative curves... those are too slow by 50% or so he says
	carriers.resize(pulse.getsamples(),double(0.));
	std::vector<double> times(carriers.size(),double(0.));
	std::vector<double> decay(carriers.size(),double(0.));
	//std::cerr << "pulse.getsamples() = " << pulse.getsamples() << "\t" << std::flush;
	//std::cerr << " pulse.modamp.size() = " << pulse.modamp.size() << "\n" << std::flush;
	size_t ncarriers_final = 1;// size_t(energy_keV*1e3/(3*bandgap_eV)); // this is a rule of thumb for exciton energy
	double wfall = 0.77863 * energy_keV; // in [fs] !!!!!!!!!!!
	double xfall = 0.0916548 * std::pow(energy_keV,int(2)) + 2.56726 * energy_keV; // in [fs] !!!!!!!!!!!
	double quad =  0.014 * std::pow(energy_keV,int(-2)); // in [fs] !!!!!!!!!!!
	double slope = 0.95 * std::pow(energy_keV,int(-2)); // in [fs] !!!!!!!!!!!
	double y0 = 0.725; // in [fs] !!!!!!!!!!!
	double arg;
	double scale;
	carriers[0] = 1.;
	decay[0] = 0.;
	for (size_t i = 1 ; i < carriers.size()/2; i++){
		arg = double(i)*pulse.getdt()*Constants::fsPau<double>(); // remember, this needs to be in femtoseconds... and getdt() returns in atomic units
		times[i] = arg;
		carriers[i] = carriers[i-1] + ((y0+slope*arg + quad * std::pow(arg,int(2))) * 0.5 * erfc((arg-xfall)/wfall));
		decay[i] = decay[i-1] + 0.5*erf((arg-xfall)/wfall);
	}
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

