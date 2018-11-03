#include "MatResponse.hpp"
#include "Constants.hpp"

MatResponse::MatResponse(MatResponse & rhs) // copy constructor
{
	t0=rhs.t0;
	twidth=rhs.twidth;
	scale = rhs.scale;
	attenuation=rhs.attenuation;
	phase=rhs.phase;
	a=rhs.a;
	b=rhs.b;
	alpha=rhs.alpha;
	beta=rhs.beta;
	etalondelay=rhs.etalondelay;
	reflectance=rhs.reflectance;
}        

MatResponse::MatResponse(double t0_in=0.0,double width_in=10.0,double atten_in = 0.95,double phase_in = 0.03) :
	t0(t0_in / fsPau<double>()),
	twidth(width_in * root_pi<double>()/ fsPau<double>() / 2.0),
	attenuation(atten_in),
	phase(phase_in)
{     
	std::cerr << "In constructor MatRespons() " << std::endl;
	a=0.75;
	b=0.25;
	alpha=5.0e-4*fsPau<double>();
	beta=0.0;
	scale=1.0;
}        
void MatResponse::setstepvec_full(PulseFreq * pulse){ setstepvec_full( *pulse); }
void MatResponse::setstepvec_full(PulseFreq & pulse){
        double arg;
        for (unsigned i = 0 ; i < pulse.getsamples(); i++){
                arg = (static_cast<double>(i)*pulse.getdt() - t0);
                if( i > pulse.getsamples()/2){
                        arg -= (pulse.getdt()*pulse.getsamples());
                }
                if( arg > twidth){
                        pulse.modphase[i] = phase;
                        pulse.modphase[i] *= a*std::exp(-1.0*(arg*alpha)) + b*exp(-1.0*(arg*beta));
                        pulse.modphase[i] *= scale;

                } else {
                        if (arg < 0.0) {
                                pulse.modphase[i] = 0.0;
                        } else {
                                pulse.modphase[i] = phase * std::pow( std::sin(half_pi<double>()*arg/twidth ), int(2) ) ;
                                pulse.modphase[i] *= a*std::exp(-1.0*(arg*alpha)) + b*std::exp(-1.0*(arg*beta));
				pulse.modphase[i] *= scale;
                        }
                }
        }
}

void MatResponse::setstepvec_amp(PulseFreq * pulse){ setstepvec_amp(*pulse); }
void MatResponse::setstepvec_amp(PulseFreq & pulse){
        double arg;
        for (unsigned i = 0 ; i < pulse.getsamples(); i++){
                arg = (static_cast<double>(i)*pulse.getdt() - t0);
                if( i > pulse.getsamples()/2){
                        arg -= (pulse.getdt()*pulse.getsamples());
                }
                if (arg > twidth) {
                        pulse.modamp[i] = -1.0*(1.0-attenuation);
                        pulse.modamp[i] *= a*std::exp(-1.0*(arg*alpha)) + b*std::exp(-1.0*(arg*beta));
                        pulse.modamp[i] += 1.0;
               
                } else {
                        if (arg < 0.0) {
                                pulse.modamp[i] = 1.0;
                        } else {
                                pulse.modamp[i] = -1.0* (1.0-attenuation) * std::pow( sin(M_PI_2*arg/twidth ) , int(2)) ;
                                pulse.modamp[i] *= a*std::exp(-1.0*(arg*alpha)) + b*std::exp(-1.0*(arg*beta));
                                pulse.modamp[i] += 1.0;
                        }
                }
        }
}

void MatResponse::setstepvec_phase(PulseFreq * pulse){ setstepvec_phase(*pulse); }
void MatResponse::setstepvec_phase(PulseFreq & pulse){
        double arg;
        for (unsigned i = 0 ; i < pulse.getsamples(); i++){
                arg = (static_cast<double>(i)*pulse.getdt() - t0);
                if( i > pulse.getsamples()/2){
                        arg -= (pulse.getdt()*pulse.getsamples());
                }
                if( arg > twidth){
                        pulse.modphase[i] = phase;
                        pulse.modphase[i] *= a*std::exp(-1.0*(arg*alpha)) + b*std::exp(-1.0*(arg*beta));
			pulse.modphase[i] *= scale;

                } else {
                        if (arg < 0.0) {
                                pulse.modphase[i] = 0.0;
                        } else {
                                pulse.modphase[i] = phase * std::pow( sin(M_PI_2*arg/twidth ) ,int(2) ) ;
                                pulse.modphase[i] *= a*std::exp(-1.0*(arg*alpha)) + b*std::exp(-1.0*(arg*beta));
				pulse.modphase[i] *= scale;
                        }
                }
        }
}

void MatResponse::buffervectors(PulseFreq * pulse){ buffervectors(*pulse); }
void MatResponse::buffervectors(PulseFreq & pulse){
        // nominally setting this to 10 times the fastest twidth to minimize its effect on spectral broadening
        unsigned bufferwidth = static_cast<unsigned>(10.0*twidth/pulse.getdt());
        double bufferscale;
        double arg;

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


void MatResponse::addstepvec_amp(PulseFreq * pulse,double delay){ addstepvec_amp(*pulse,delay); }
void MatResponse::addstepvec_amp(PulseFreq & pulse,double delay){
        double arg;
        double thisamp;
        for (unsigned i = 0 ; i < pulse.getsamples(); i++){
                arg = (static_cast<double>(i)*pulse.getdt() - (t0-delay));
                if( i > pulse.getsamples()/2){
                        arg -= (pulse.getdt()*pulse.getsamples());
                }
                if (arg > twidth) {
                        thisamp = -1.0*(1.0-attenuation);
                        thisamp *= a*exp(-1.0*(arg*alpha)) + b*exp(-1.0*(arg*beta));
			thisamp *= scale;
                        thisamp += 1.0;

                } else {
                        if (arg < 0.0) {
                                thisamp = 1.0;
                        } else {
                                thisamp = -1.0*(1.0-attenuation) * std::pow( sin(M_PI_2*arg/twidth ) , int(2) ) ;
                                thisamp *= a*exp(-1.0*(arg*alpha)) + b*exp(-1.0*(arg*beta));
				thisamp *= scale;
                                thisamp += 1.0;
                        }
                }
                pulse.modamp[i] *= thisamp;
        }

}


void MatResponse::addstepvec_phase(PulseFreq * pulse,double delay){ addstepvec_phase(*pulse,delay); }
void MatResponse::addstepvec_phase(PulseFreq & pulse,double delay){
        double arg;
        double thisphase;
        for (unsigned i = 0 ; i < pulse.getsamples(); i++){
                arg = (static_cast<double>(i)*pulse.getdt() - (t0-delay));
                if( i > pulse.getsamples()/2){
                        arg -= (pulse.getdt()*pulse.getsamples());
                }
                if( arg > twidth){
                        thisphase = phase;
                        thisphase *= a*exp(-1.0*(arg*alpha)) + b*exp(-1.0*(arg*beta));
			thisphase *= scale;

                } else {
                        if (arg < 0.0) {
                                thisphase = 0.0;
                        } else {
                                thisphase = phase * std::pow( sin(M_PI_2*arg/twidth ) , int(2)) ;
                                thisphase *= a*exp(-1.0*(arg*alpha)) + b*exp(-1.0*(arg*beta));
				thisphase *= scale;
                        }
                }
                pulse.modphase[i] += thisphase;
        }

}

