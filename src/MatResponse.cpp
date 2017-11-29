#include "MatResponse.hpp"

void MatResponse::setstepvec_full(PulseFreq & pulse){ setstepvec_full( &pulse); }
void MatResponse::setstepvec_full(PulseFreq * pulse){
        double arg;
        for (unsigned i = 0 ; i < pulse->getsamples(); i++){
                arg = (static_cast<double>(i)*pulse->getdt() - t0);
                if( i > pulse->getsamples()/2){
                        arg -= (pulse->getdt()*pulse->getsamples());
                }
                if( arg > twidth){
                        pulse->modphase->data[i] = phase;
                        pulse->modphase->data[i] *= a*std::exp(-1.0*(arg*alpha)) + b*exp(-1.0*(arg*beta));

                } else {
                        if (arg < 0.0) {
                                pulse->modphase->data[i] = 0.0;
                        } else {
                                pulse->modphase->data[i] = phase * std::pow( std::sin(M_PI_2*arg/twidth ), int(2) ) ;
                                pulse->modphase->data[i] *= a*std::exp(-1.0*(arg*alpha)) + b*std::exp(-1.0*(arg*beta));
                        }
                }
        }
}

void MatResponse::setstepvec_amp(PulseFreq & pulse){ setstepvec_amp(&pulse); }
void MatResponse::setstepvec_amp(PulseFreq * pulse){
        double arg;
        for (unsigned i = 0 ; i < pulse->getsamples(); i++){
                arg = (static_cast<double>(i)*pulse->getdt() - t0);
                if( i > pulse->getsamples()/2){
                        arg -= (pulse->getdt()*pulse->getsamples());
                }
                if (arg > twidth) {
                        pulse->modamp->data[i] = -1.0*(1.0-attenuation);
                        pulse->modamp->data[i] *= a*std::exp(-1.0*(arg*alpha)) + b*std::exp(-1.0*(arg*beta));
                        pulse->modamp->data[i] += 1.0;
               
                } else {
                        if (arg < 0.0) {
                                pulse->modamp->data[i] = 1.0;
                        } else {
                                pulse->modamp->data[i] = -1.0* (1.0-attenuation) * std::pow( sin(M_PI_2*arg/twidth ) , int(2)) ;
                                pulse->modamp->data[i] *= a*std::exp(-1.0*(arg*alpha)) + b*std::exp(-1.0*(arg*beta));
                                pulse->modamp->data[i] += 1.0;
                        }
                }
        }
}
/*
void MatResponse::setstepvec_amp(gsl_vector * modamp,double dt){
	double arg;
	for (unsigned i = 0 ; i < modamp->size; i++){
		arg = (static_cast<double>(i)*dt - t0);
		if( i > modamp->size/2){
			arg -= (dt*modamp->size);
		}
		if (arg > twidth) {
			modamp->data[i] = -1.0*(1.0-attenuation);
			modamp->data[i] *= a*std::exp(-1.0*(arg*alpha)) + b*std::exp(-1.0*(arg*beta));
			modamp->data[i] += 1.0;
			
		} else {
			if (arg < 0.0) {
				modamp->data[i] = 1.0;
			} else {
				modamp->data[i] = -1.0* (1.0-attenuation) * std::pow( std::sin(M_PI_2*arg/twidth ) , int(2)) ; 
				modamp->data[i] *= a*std::exp(-1.0*(arg*alpha)) + b*std::exp(-1.0*(arg*beta));
				modamp->data[i] += 1.0;
			}
		}
	}
}
*/
void MatResponse::setstepvec_phase(PulseFreq & pulse){ setstepvec_phase(&pulse); }
void MatResponse::setstepvec_phase(PulseFreq * pulse){
        double arg;
        for (unsigned i = 0 ; i < pulse->getsamples(); i++){
                arg = (static_cast<double>(i)*pulse->getdt() - t0);
                if( i > pulse->getsamples()/2){
                        arg -= (pulse->getdt()*pulse->getsamples());
                }
                if( arg > twidth){
                        pulse->modphase->data[i] = phase;
                        pulse->modphase->data[i] *= a*std::exp(-1.0*(arg*alpha)) + b*std::exp(-1.0*(arg*beta));

                } else {
                        if (arg < 0.0) {
                                pulse->modphase->data[i] = 0.0;
                        } else {
                                pulse->modphase->data[i] = phase * std::pow( sin(M_PI_2*arg/twidth ) ,int(2) ) ;
                                pulse->modphase->data[i] *= a*std::exp(-1.0*(arg*alpha)) + b*std::exp(-1.0*(arg*beta));
                        }
                }
        }
}

/*
void MatResponse::setstepvec_phase(gsl_vector * modphase,double dt){
	double arg;
	for (unsigned i = 0 ; i < modphase->size; i++){
		arg = (static_cast<double>(i)*dt - t0);
		if( i > modphase->size/2){
			arg -= (dt*modphase->size);
		}
		if( arg > twidth){
			modphase->data[i] = phase;
			modphase->data[i] *= a*std::exp(-1.0*(arg*alpha)) + b*std::exp(-1.0*(arg*beta));
			
		} else {
			if (arg < 0.0) {
				modphase->data[i] = 0.0;
			} else {
				modphase->data[i] = phase * std::pow( std::sin(M_PI_2*arg/twidth ) ,int(2)) ; 
				modphase->data[i] *= a*std::exp(-1.0*(arg*alpha)) + b*std::exp(-1.0*(arg*beta));
			}
		}
	}
}

void MatResponse::buffervectors(gsl_vector * modamp,gsl_vector * modphase,double dt){
	// nominally setting this to 10 times the fastest twidth to minimize its effect on spectral broadening
	unsigned bufferwidth = static_cast<unsigned>(10.0*twidth/dt);
	double bufferscale;
	double arg;
	if (modphase->size != modamp->size)
		cerr << "vectors differently sized in MaResponse::buffervectors( call" << endl;

	modphase->data[modphase->size/2] *= 0.0;
	modamp->data[modamp->size/2] = 1.0;
	for (unsigned i=0; i < bufferwidth ; i++) {
		arg = ((double)i)/((double)bufferwidth);
		bufferscale = gsl_pow_2( sin(M_PI_2*(arg) ) ) ; // buffering the vector back to 0.0 so it can be periodic.
		modphase->data[modphase->size/2-i] *= bufferscale;
		modamp->data[modamp->size/2-i] -= 1.0;
		modamp->data[modamp->size/2-i] *= bufferscale;
		modamp->data[modamp->size/2-i] += 1.0;
	}
}
*/

void MatResponse::buffervectors(PulseFreq & pulse){ buffervectors(&pulse); }
void MatResponse::buffervectors(PulseFreq * pulse){
        // nominally setting this to 10 times the fastest twidth to minimize its effect on spectral broadening
        unsigned bufferwidth = static_cast<unsigned>(10.0*twidth/pulse->getdt());
        double bufferscale;
        double arg;

        pulse->modphase->data[pulse->getsamples()/2] *= 0.0;
        pulse->modamp->data[pulse->getsamples()/2] = 1.0;
        for (unsigned i=0; i < bufferwidth ; i++) {
                arg = ((double)i)/((double)bufferwidth);
                bufferscale = std::pow( sin(M_PI_2*(arg) ) , int(2)) ; // buffering the vector back to 0.0 so it can be periodic.
                pulse->modphase->data[pulse->getsamples()/2-i] *= bufferscale;
                pulse->modamp->data[pulse->getsamples()/2-i] -= 1.0;
                pulse->modamp->data[pulse->getsamples()/2-i] *= bufferscale;
                pulse->modamp->data[pulse->getsamples()/2-i] += 1.0;
        }
}


/*
void MatResponse::addstepvec_amp(gsl_vector * modamp,double dt,double delay){
        double arg;
	double thisamp;
        for (unsigned i = 0 ; i < modamp->size; i++){
                arg = (static_cast<double>(i)*dt - (t0-delay));
                if( i > modamp->size/2){
                        arg -= (dt*modamp->size);
                }
                if (arg > twidth) {
                        thisamp = -1.0*(1.0-attenuation);
                        thisamp *= a*exp(-1.0*(arg*alpha)) + b*exp(-1.0*(arg*beta));
                        thisamp += 1.0;

                } else {
                        if (arg < 0.0) {
                                thisamp = 1.0;
                        } else {
                                thisamp = -1.0* (1.0-attenuation) * gsl_pow_2( sin(M_PI_2*arg/twidth ) ) ;
                                thisamp *= a*exp(-1.0*(arg*alpha)) + b*exp(-1.0*(arg*beta));
                                thisamp += 1.0;
                        }
                }
		modamp->data[i] *= thisamp;
        }

}
*/
void MatResponse::addstepvec_amp(PulseFreq & pulse,double delay){ addstepvec_amp(&pulse,delay); }
void MatResponse::addstepvec_amp(PulseFreq * pulse,double delay){
        double arg;
        double thisamp;
        for (unsigned i = 0 ; i < pulse->getsamples(); i++){
                arg = (static_cast<double>(i)*pulse->getdt() - (t0-delay));
                if( i > pulse->getsamples()/2){
                        arg -= (pulse->getdt()*pulse->getsamples());
                }
                if (arg > twidth) {
                        thisamp = -1.0*(1.0-attenuation);
                        thisamp *= a*exp(-1.0*(arg*alpha)) + b*exp(-1.0*(arg*beta));
                        thisamp += 1.0;

                } else {
                        if (arg < 0.0) {
                                thisamp = 1.0;
                        } else {
                                thisamp = -1.0* (1.0-attenuation) * std::pow( sin(M_PI_2*arg/twidth ) , int(2) ) ;
                                thisamp *= a*exp(-1.0*(arg*alpha)) + b*exp(-1.0*(arg*beta));
                                thisamp += 1.0;
                        }
                }
                pulse->modamp->data[i] *= thisamp;
        }

}

/*
void MatResponse::addstepvec_phase(gsl_vector * modphase,double dt,double delay){
        double arg;
	double thisphase;
        for (unsigned i = 0 ; i < modphase->size; i++){
                arg = (static_cast<double>(i)*dt - (t0-delay));
                if( i > modphase->size/2){
                        arg -= (dt*modphase->size);
                }
                if( arg > twidth){
			thisphase = phase;
			thisphase *= a*exp(-1.0*(arg*alpha)) + b*exp(-1.0*(arg*beta));
                        
                } else {
                        if (arg < 0.0) {
                                thisphase = 0.0;
                        } else {        
                                thisphase = phase * gsl_pow_2( sin(M_PI_2*arg/twidth ) ) ;
                                thisphase *= a*exp(-1.0*(arg*alpha)) + b*exp(-1.0*(arg*beta));
                        }
                }
		modphase->data[i] += thisphase;
        }

}
*/
void MatResponse::addstepvec_phase(PulseFreq & pulse,double delay){ addstepvec_phase(&pulse,delay); }
void MatResponse::addstepvec_phase(PulseFreq * pulse,double delay){
        double arg;
        double thisphase;
        for (unsigned i = 0 ; i < pulse->getsamples(); i++){
                arg = (static_cast<double>(i)*pulse->getdt() - (t0-delay));
                if( i > pulse->getsamples()/2){
                        arg -= (pulse->getdt()*pulse->getsamples());
                }
                if( arg > twidth){
                        thisphase = phase;
                        thisphase *= a*exp(-1.0*(arg*alpha)) + b*exp(-1.0*(arg*beta));

                } else {
                        if (arg < 0.0) {
                                thisphase = 0.0;
                        } else {
                                thisphase = phase * std::pow( sin(M_PI_2*arg/twidth ) , int(2)) ;
                                thisphase *= a*exp(-1.0*(arg*alpha)) + b*exp(-1.0*(arg*beta));
                        }
                }
                pulse->modphase->data[i] += thisphase;
        }

}


/*
void MatResponse::setstepvec_marco(gsl_vector * modamp,double dt){
	double arg;
	const double aa = 0.09;
	const double sig = 20.0/fsPau<float>();
	const double tau_formation = 15.0/fsPau<float>();
	const double fraction_nodecay = 0.6;
	const double tau_decay = 150.0/fsPau<float>();

        for (unsigned i = 0 ; i < modamp->size; i++){
                arg = (static_cast<double>(i)*dt - t0);
		if( i > modamp->size/2){
                        arg -= (dt*modamp->size);
                }
		modamp->data[i] = gsl_sf_erf(arg/sig/M_SQRT2);
		modamp->data[i] -= exp(-(2.0*tau_formation*arg - gsl_pow_2(sig))/2.0/gsl_pow_2(tau_formation))
					*( 1.0-gsl_sf_erf((-tau_formation*arg+gsl_pow_2(sig))/M_SQRT2/sig/tau_formation) );		
		modamp->data[i] *= aa*0.5;
		modamp->data[i] *= (1.0-fraction_nodecay)*exp(-arg/tau_decay + fraction_nodecay);
	}
}
*/
/*
def
step(t,a=0.02,t0=0,sig=10e-15,tau_formation=10e-15,fraction_nodecay=0.8,tau_decay=200e-15):
   s  =  a*0.5*(special.erf((t-t0)/sig/sqrt2)+1)
   s +=
-a*0.5*np.exp(-(2*tau_formation*(t-t0)-sig**2)/2/tau_formation**2)*(1-special.erf(
(-tau_formation*(t-t0)+sig**2)/sqrt2/tau_formation/sig))
   return s*(
(1-fraction_nodecay)*np.exp(-(t-t0)/tau_decay)+fraction_nodecay )
*/

