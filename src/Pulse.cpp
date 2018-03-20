// Pulse clas implimentation
// standard includes
#include <math.h>


// my headers
#include <Pulse.hpp>
#include <boost/math/interpolators/barycentric_rational.hpp>
#include <fftw3.h>
#include <algorithm>
#include <Constants.hpp>
#include <boost/lexical_cast.hpp>
#include <DataOps.hpp>

/*
PulseFreq & PulseFreq::operator=(const PulseFreq &rhs)
	this->omega_center(rhs.omega_center),
        this->omega_width(rhs.omega_width),
	this->omega_high(rhs.omega_high),
	this->omega_onwidth(rhs.omega_onwidth),
        this->domega(rhs.domega),
        this->intime(rhs.intime),
        this->infreq(rhs.infreq),
        this->omega(rhs.omega), // these are static vectors... i'm trying not to make copies just yet
        this->time(rhs.time), // these are static vectors... i'm trying not to make copies just yet
	this->parent(false),
	this->child(true),
	this->i_low( static_cast<unsigned>(static_cast<double>( atof( getenv("nu_low") ) )* 2.0*M_PI*fsPau/domega) ),
	this->i_high( static_cast<unsigned>(static_cast<double>( atof( getenv("nu_high") ) )* 2.0*M_PI*fsPau/domega) )
	{
	this->samples = rhs.sample
	this->startind = rhs.startind;this->stopind=rhs.stopind;this->onwidth=rhs.onwidth;this->offwidth=rhs.offwidth;
	this->tspan = rhs.tspan;
	this->lambda_center=rhs.lambda_center;this->lambda_width=rhs.lambda_width;
	this->phase_GDD=rhs.phase_GDD;this->phase_TOD=rhs.phase_TOD;this->phase_4th=rhs.phase_4th;this->phase_5th=rhs.phase_5th;

        this->dtime = rhs.dtime;this->time_center=rhs.time_center;this->time_wdith=rhs.time_wdith;

	this->cwave=rhs.cwave;this->cspace=rhs.cspace;
	this->cvec = gsl_vector_complex_alloc(samples);
	this->rhovec = gsl_vector_alloc(samples);
	this->phivec = gsl_vector_alloc(samples);
        this->nu0=rhs.nu0;
		
	// point tables to parent //
	// copy vectors explicitly //
	gsl_vector_memcpy(this->rhovec,rhs.rhovec);
	gsl_vector_memcpy(this->phivec,rhs.phivec);
	gsl_vector_complex_memcpy(this->cvec,rhs.cvec);

        return *this;
}
*/

using namespace Constants;

PulseFreq * PulseFreq::operator+=(const PulseFreq *rhs){
	gsl_vector_complex_add(cvec,rhs->cvec);
	for (unsigned i=0;i<samples;i++){
		gsl_vector_set(rhovec,i,gsl_complex_abs(gsl_vector_complex_get(cvec,i)));
		gsl_vector_set(phivec,i,gsl_complex_arg(gsl_vector_complex_get(cvec,i)));
	}
	return this;
}

PulseFreq & PulseFreq::operator+=(const PulseFreq &rhs){
	gsl_vector_complex_add(cvec,rhs.cvec);
	for (unsigned i=0;i<samples;i++){
		gsl_vector_set(rhovec,i,gsl_complex_abs(gsl_vector_complex_get(cvec,i)));
		gsl_vector_set(phivec,i,gsl_complex_arg(gsl_vector_complex_get(cvec,i)));
	}
	return *this;
}

PulseFreq & PulseFreq::operator-=(const PulseFreq &rhs){
        gsl_vector_complex_sub(cvec,rhs.cvec);
        for (unsigned i=0;i<samples;i++){
                gsl_vector_set(rhovec,i,gsl_complex_abs(gsl_vector_complex_get(cvec,i)));
                gsl_vector_set(phivec,i,gsl_complex_arg(gsl_vector_complex_get(cvec,i)));
        }
        return *this;
}

PulseFreq & PulseFreq::operator*=(const PulseFreq &rhs){
	gsl_vector_complex_mul(cvec,rhs.cvec);
	for (unsigned i=0;i<samples;i++){
		gsl_vector_set(rhovec,i,gsl_complex_abs(gsl_vector_complex_get(cvec,i)));
		gsl_vector_set(phivec,i,gsl_complex_arg(gsl_vector_complex_get(cvec,i)));
	}
	return *this;
}

PulseFreq & PulseFreq::operator*=(const double scale){
	gsl_complex z=gsl_complex_polar(scale,0.0);
	gsl_vector_complex_scale(cvec,z);
	for (unsigned i=0;i<samples;i++){
		gsl_vector_set(rhovec,i,gsl_complex_abs(gsl_vector_complex_get(cvec,i)));
		gsl_vector_set(phivec,i,gsl_complex_arg(gsl_vector_complex_get(cvec,i)));
	}
	return *this;
}

PulseFreq & PulseFreq::operator/=(const PulseFreq &rhs){
        gsl_vector_complex_div(cvec,rhs.cvec);
        for (unsigned i=0;i<samples;i++){
                gsl_vector_set(rhovec,i,gsl_complex_abs(gsl_vector_complex_get(cvec,i)));
                gsl_vector_set(phivec,i,gsl_complex_arg(gsl_vector_complex_get(cvec,i)));
        }
        return *this;
}

PulseFreq & PulseFreq::interfere(const PulseFreq &rhs){
	*this += rhs;
        return *this;
}

PulseFreq & PulseFreq::diffamps(const PulseFreq &rhs){
	gsl_complex z;
        gsl_vector_sub(rhovec,rhs.rhovec);
        for (unsigned i=0;i<samples;i++){
		z = gsl_complex_rect(gsl_vector_get(rhovec,i) , 0.0);
                gsl_vector_complex_set(cvec,i,z);
                gsl_vector_set_zero(phivec);
        }
        return *this;
}

PulseFreq & PulseFreq::normamps(const PulseFreq &rhs){
	gsl_complex z;
        gsl_vector_div(rhovec,rhs.rhovec);
        for (unsigned i=0;i<samples;i++){
		GSL_SET_COMPLEX(&z,gsl_vector_get(rhovec,i),0.0);
                gsl_vector_complex_set(cvec,i,z);
                gsl_vector_set_zero(phivec);
        }
        return *this;
}

void PulseFreq::print_amp(std::ofstream & outfile)
{
	outfile << "# amp\n";
	for (size_t i = 0; i<samples;++i){
		outfile << rhovec->data[i] << "\n";
	}
	outfile << std::endl;
}
void PulseFreq::print_phase(std::ofstream & outfile)
{
	outfile << "# phase\n";
	for (size_t i = 0; i<samples;++i){
		outfile << phivec->data[i] << "\n";
	}
	outfile << std::endl;
}
void PulseFreq::print_phase_powerspectrum(std::ofstream & outfile)
{
	double * phase = (double *) fftw_malloc(sizeof(double) * samples);
	double * phaseFT = (double *) fftw_malloc(sizeof(double) * samples);
	fftw_plan plan_r2hc = fftw_plan_r2r_1d(samples,
			phase,
			phaseFT,
			FFTW_R2HC,
			FFTW_MEASURE
			);
	fftw_execute_r2r(plan_r2hc,phase,phaseFT);

	outfile << "# power spectrum of the phaseFT\n";
	outfile << std::pow(phaseFT[0],int(2)) << "\n";
	for (size_t i = 1; i<samples/2;++i){
		outfile << std::pow(phaseFT[i],int(2)) + std::pow(phaseFT[samples-i],int(2)) << "\n";
	}
	outfile << std::pow(phaseFT[samples/2],int(2)) << std::endl;
	outfile << std::endl;
}

bool PulseFreq::addrandomphase(void)
{
	if (!infreq){
		cerr << "died here at addrandomphase()" << endl;
		return false;
	}
	size_t sz = samples*2; // doubling the vector to mirror it so that DFT hansles the phase well

	double * randphase = (double *) fftw_malloc(sizeof(double) * sz);
	double * randphaseFT = (double *) fftw_malloc(sizeof(double) * sz);

	fftw_plan plan_r2hc = fftw_plan_r2r_1d(sz,
			randphase,
			randphaseFT,
			FFTW_R2HC,
			FFTW_MEASURE
			);
	fftw_plan plan_hc2r = fftw_plan_r2r_1d(sz,
			randphaseFT,
			randphase,
			FFTW_HC2R,
			FFTW_MEASURE
			);


	std::random_device rng;
	/*
	std::normal_distribution<double> norm_distribution(
		double(atof(getenv("randphase_mean")))*Constants::pi<double>(),
		double(atof(getenv("randphase_std")))*Constants::pi<double>()
		);
	*/
	std::uniform_real_distribution<double> distribution(
		(double(atof(getenv("randphase_mean")))-double(atof(getenv("randphase_std"))))*Constants::pi<double>(),
		(double(atof(getenv("randphase_mean")))+double(atof(getenv("randphase_std"))))*Constants::pi<double>()
		);

	double phase = distribution(rng);
	randphase[0] = phase;
	randphase[samples/2] = -phase;
	for (size_t i = 1; i<samples/2;i++){
		phase = distribution(rng);
		randphase[i] = phase;
		randphase[samples-i] = -phase;
	}
	for (size_t i=sz-1;i>sz/2-1;--i){
		randphase[i] = randphase[sz-i];
	}

	size_t lowpass = boost::lexical_cast<size_t>(atoi(getenv("phaseNoiseLowpass")));
	std::cerr << "\n\n======== lowpass is " << lowpass << " =======\n\n" << std::flush;

	fftw_execute_r2r(plan_r2hc,randphase,randphaseFT);
	std::fill(randphaseFT+lowpass,randphaseFT+sz-lowpass,0.);
	for (size_t i=1;i<lowpass;++i){
		double filter = std::pow(std::cos(double(i)/(double(lowpass)) * Constants::half_pi<double>() ),int(2));
		randphaseFT[i] *= filter;
		randphaseFT[sz-i] *= filter;
	}
	randphaseFT[sz/2] = 0.;
	fftw_execute_r2r(plan_hc2r,randphaseFT,randphase);

	for (size_t i=0;i<samples;++i){
		phivec->data[i] += randphase[i]/samples;
	}

	gsl_complex z;
	z = gsl_complex_polar(rhovec->data[0],phivec->data[0]);
	gsl_vector_complex_set(cvec,0,z);
	z = gsl_complex_polar(rhovec->data[samples/2],phivec->data[samples/2]);
	gsl_vector_complex_set(cvec,samples/2,z);
	for (size_t i = 1; i<samples/2;i++){
		z = gsl_complex_polar(rhovec->data[i],phivec->data[i]);
		gsl_vector_complex_set(cvec,i,z);
		z = gsl_complex_polar(rhovec->data[samples-i],phivec->data[samples-i]);
		gsl_vector_complex_set(cvec,samples-i,z);
	}
	return true;
}


void PulseTime::setstrength(const double in)
{
  strength = in * auenergy<float>()/Eh<float>() * std::pow(aufor10PW<float>(),int(2));
}

void PulseTime::setwidth(const double in)
{
  Ctau = in * M_SQRTPI / fsPau<float>() / 2.0;
}

void PulseTime::sett0(const double in)
{
  t0 = in / fsPau<float>();
}

PulseFreq::PulseFreq(const double omcenter_in=(0.55*fsPau<float>()),const double omwidth_in=(0.15*fsPau<float>()),const double omonoff_in=(0.1*fsPau<float>()), double tspan_in=(10000.0/fsPau<float>())): // default constructor
	omega_center(omcenter_in),
	omega_width(omwidth_in ),
	omega_high( GSL_MAX_DBL(4.0*(omcenter_in + omwidth_in),10.0*omcenter_in) ),
	domega( 2.0*M_PI/tspan_in),
	omega(NULL),
	intime(false),
	infreq(true),
	time(NULL),
	parent(true),
	child(false),
	i_low( static_cast<unsigned>(static_cast<double>( atof( getenv("nu_low") ) )* 2.0*M_PI*fsPau<float>()/domega) ),
	i_high( static_cast<unsigned>(static_cast<double>( atof( getenv("nu_high") ) )* 2.0*M_PI*fsPau<float>()/domega) ),
	m_noisescale(1e-3),
	m_sampleinterval(2),
	m_saturate(4096),
	m_gain(1000000),
	m_lamsamples(2048)
{
	samples = (( static_cast<unsigned>(2.0 * omega_high / domega))/SAMPLEROUND + 1 ) *SAMPLEROUND;// dt ~ .1fs, Dt ~ 10000fs, dom = 2*pi/1e4, omhigh = 2*pi/.1fs, samples = omhigh/dom*2
	allocatetables();
	dtime = tspan_in/static_cast<double>(samples);
	omega_onwidth = omega_offwidth = omega_width/2.0; // forcing sin2 gaussian spectrum
	buildvectors();
	nu0=omcenter_in/(2.0*M_PI)*fsPau<float>();
	phase_GDD=phase_TOD=phase_4th=phase_5th=0.0;
	m_gain = (unsigned long)(atoi( getenv("gain")));
	m_noisescale = (double)( atof( getenv("noisescale") ) );
	m_sampleinterval = (unsigned)(atoi(getenv("sampleinterval")));
	m_saturate = uint16_t( atoi( getenv("saturate")));
}
PulseFreq::PulseFreq(PulseFreq &rhs): // copy constructor
	omega_center(rhs.omega_center),
	omega_width(rhs.omega_width),
	omega_high(rhs.omega_high),
	omega_onwidth(rhs.omega_onwidth),
	domega(rhs.domega),
	intime(rhs.intime),
	infreq(rhs.infreq),
	omega(rhs.omega), // these are static vectors... i'm trying not to make copies just yet
	time(rhs.time), // these are static vectors... i'm trying not to make copies just yet
	parent(false),
	child(true),
	i_low(rhs.i_low), // static_cast<unsigned>(static_cast<double>( atof( getenv("nu_low") ) )* 2.0*M_PI*fsPau/domega) ),
	i_high(rhs.i_high), /// static_cast<unsigned>(static_cast<double>( atof( getenv("nu_high") ) )* 2.0*M_PI*fsPau/domega) ),
	m_noisescale(rhs.m_noisescale),
	m_sampleinterval(rhs.m_sampleinterval),
	m_saturate(rhs.m_saturate),
	m_gain(rhs.m_gain),
	m_lamsamples(rhs.m_lamsamples)
{
	samples = rhs.samples;
	startind = rhs.startind;stopind=rhs.stopind;onwidth=rhs.onwidth;offwidth=rhs.offwidth;
	tspan = rhs.tspan;
	lambda_center=rhs.lambda_center;lambda_width=rhs.lambda_width;
	phase_GDD=rhs.phase_GDD;phase_TOD=rhs.phase_TOD;phase_4th=rhs.phase_4th;phase_5th=rhs.phase_5th;

	dtime = rhs.dtime;time_center=rhs.time_center;time_wdith=rhs.time_wdith;

	cwave=rhs.cwave;
	cspace = gsl_fft_complex_workspace_alloc(samples);
	cvec = gsl_vector_complex_alloc(samples);
	rhovec = gsl_vector_alloc(samples);
	phivec = gsl_vector_alloc(samples);
	modamp = gsl_vector_alloc(samples);
	modphase = gsl_vector_alloc(samples);

	nu0=rhs.nu0;

	// point tables to parent //
	// copy vectors explicitly //
	gsl_vector_memcpy(rhovec,rhs.rhovec);
	gsl_vector_memcpy(phivec,rhs.phivec);
	gsl_vector_complex_memcpy(cvec,rhs.cvec);
	gsl_vector_memcpy(modamp,rhs.modamp);
	gsl_vector_memcpy(modphase,rhs.modphase);

}

PulseFreq::~PulseFreq(void){
	if(child){
		killtheseonly();
	} else {
		killvectors();
	}
}

void PulseFreq::attenuate(double attenfactor){
	gsl_complex z;
	for( unsigned i=0;i<samples;i++){
		rhovec->data[i] *= attenfactor;
		z = gsl_complex_polar(rhovec->data[i],phivec->data[i]);
		gsl_vector_complex_set(cvec,i,z);
	}
}
void PulseFreq::phase(double phasein){ // expects delay in units of pi , i.e. 1.0 = pi phase flip 
	gsl_complex z;
	if(intime){
		fft_tofreq();
	}
	for (unsigned i=0;i<samples;i++){
		phivec->data[i] += phasein*M_PI;
		z = gsl_complex_polar(rhovec->data[i],phivec->data[i]);
		gsl_vector_complex_set(cvec,i,z);
	}
	if(intime){
		fft_totime();
	}
}
void PulseFreq::delay(double delayin){ // expects delay in fs
	/* 
	Somehow, this seems to casue errors, as only when etalon effect is added do we see the artifact.
	But how, it looks from this that we are adding to an existing delay... maybe mis-handeling the negative frequencies?
	*/
	gsl_complex z;
	if(intime){
		fft_tofreq();
	}
	for (unsigned i=0;i<samples;i++){
		phivec->data[i] += omega[i]*delayin/fsPau<float>();
		z = gsl_complex_polar(rhovec->data[i],phivec->data[i]);
		gsl_vector_complex_set(cvec,i,z);
	}
	if(intime){
		fft_totime();
	}
}


void PulseFreq::printfrequency(ofstream * outfile){
	double nu,lambda;
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*M_PI)/fsPau<float>());
		lambda = C_nmPfs<float>()/nu;
		(*outfile) << nu << "\t" << lambda << "\t" << rhovec->data[i] << "\t" << phivec->data[i] << "\n";
	}
}
void PulseFreq::printwavelengthbins(ofstream * outfile)
{
	std::vector<double> x(2);
	x.front() = C_nmPfs<float>()*2.0*M_PI*fsPau<float>()/omega[i_low];
	x.back() = C_nmPfs<float>()*2.0*M_PI*fsPau<float>()/omega[i_high-1];
	double dlam = (x.front()-x.back())/double(m_lamsamples);
        for (size_t i = 0;i<m_lamsamples;++i){
		(*outfile) << x.back() + i*dlam << "\t";
        }
	(*outfile) << "\n";
	return;
}
void PulseFreq::appendwavelength(ofstream * outfile)
{
	std::vector<double> x(i_high-i_low);
	std::vector<double> y(i_high-i_low);	
	for (size_t i=0;i<y.size();++i){
		x[i] = C_nmPfs<float>()*2.0*M_PI*fsPau<float>()/omega[i_low+i];
		//y[i] = std::min(rhovec->data[i_low+i] * m_gain,double(m_saturate));
		y[i] = rhovec->data[i_low+i];
		if (550<i && i<570)
			std::cerr << x[i] << "," << y[i] << "\t";

	}
	double dlam = (x.front()-x.back())/double(m_lamsamples);
	boost::math::barycentric_rational<double> interpolant(x.data(), y.data(), y.size());
	for (size_t i=0;i<m_lamsamples;++i){
		(*outfile) << uint16_t(interpolant(x.back()+i*dlam)) << "\t";
	}
	(*outfile) << std::endl;
	return;
}
void PulseFreq::appendfrequency(ofstream * outfile){
        for (unsigned i = i_low;i<i_high;i+=m_sampleinterval){ 
		uint16_t val = std::min(unsigned(rhovec->data[i] * m_gain),unsigned(m_saturate));
       		(*outfile) << val << "\t";
        }
	(*outfile) << std::endl;
}

void PulseFreq::appendnoisy(ofstream * outfile){
	time_t timeseed = std::time(NULL);
	double nu,outval;
	const gsl_rng_type * rngType = gsl_rng_taus;
	gsl_rng * rngPtr = gsl_rng_alloc(rngType);
	gsl_rng_set(rngPtr,timeseed);
        for (unsigned i = i_low;i<i_high;i+=m_sampleinterval){ 
		outval = rhovec->data[i];
		outval += gsl_ran_gaussian(rngPtr,m_noisescale);
       		(*outfile) << outval << "\t" ;
        }
	(*outfile) << std::endl;
	gsl_rng_free(rngPtr);
}

void PulseFreq::printfrequencybins(ofstream * outfile){
	double nu,lam;
        for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*M_PI)/fsPau<float>());
		lam = C_nmPfs<float>()/nu;
       		(*outfile) << nu << "\t" << lam << "\n";
        }
	(*outfile) << "\n";
}
void PulseFreq::appendfrequencybins(ofstream * outfile){
	double nu,lam;
        for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*M_PI)/fsPau<float>());
       		(*outfile) << nu << "\t";
        }
	(*outfile) << "\n";
}
void PulseFreq::printfrequencydelay(ofstream * outfile, const double *delay){
	double nu,lambda,thisdelay;
	thisdelay = (*delay)*fsPau<float>();
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\tdelay[fs]\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*M_PI)/fsPau<float>());
		lambda = C_nmPfs<float>()/nu;
		(*outfile) << nu << "\t" << lambda << "\t" << rhovec->data[i] << "\t" << phivec->data[i] << "\t" << thisdelay << "\n";
	}
	(*outfile) << "\n";
}
void PulseFreq::printfrequencydelaychirp(ofstream * outfile, const double *delay,const double *chirp){
	double nu,lambda,thisdelay,thischirp,reltime;
	thisdelay = (*delay)*fsPau<float>();
	thischirp = (*chirp)*gsl_pow_2(fsPau<float>());
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\tdelay[fs]\tchirp[fs^2]\treldelays[fs]\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*M_PI)/fsPau<float>());
		lambda = C_nmPfs<float>()/nu;
		reltime = (nu-nu0)*thischirp*10.0;
		(*outfile) << nu << "\t" << lambda << "\t" << rhovec->data[i] << "\t" << phivec->data[i] << "\t" << thisdelay << "\t" << thischirp << "\t" << reltime << "\n";
	}
	(*outfile) << "\n";
}

void PulseFreq::printwavelength(ofstream * outfile,const double *delay){
	double nu,lambda,lamlast;
	lamlast=0;
        for (unsigned i = i_high;i>i_low;i-=(unsigned)(atoi(getenv("sampleinterval")))){
                nu = (omega[i]/(2.0*M_PI)/fsPau<float>());
                lambda = (double)((int)((C_nmPfs<float>()/nu)*10.0))/10.0;
		if (lambda>lamlast & (int)(lambda*10)%5 == 0){
                	(*outfile) << lambda << "\t" << (*delay) << "\t" << rhovec->data[i] << "\n";
			lamlast = lambda;
		}
        }
	(*outfile) << "\n";
}

void PulseFreq::printtime(ofstream * outfile){
	(*outfile) << ("#time\treal\timag\n");
	/*
	for (unsigned i = 0;i<samples;i++){
		(*outfile) << (time[i]*fsPau) << "\t" << GSL_REAL(gsl_vector_complex_get(cvec,i)) << "\t" << GSL_IMAG(gsl_vector_complex_get(cvec,i)) << "\n";
	}
	*/
	for (unsigned i = samples/2;i<samples;i++){
		(*outfile) << (time[i]*fsPau<float>()) << "\t" << GSL_REAL(gsl_vector_complex_get(cvec,i)) << "\t" << GSL_IMAG(gsl_vector_complex_get(cvec,i)) << "\n";
	}
	for (unsigned i = 0;i<samples/2;i++){
		(*outfile) << (time[i]*fsPau<float>()) << "\t" << GSL_REAL(gsl_vector_complex_get(cvec,i)) << "\t" << GSL_IMAG(gsl_vector_complex_get(cvec,i)) << "\n";
	}

} 


void PulseFreq::allocatetables(void){
	cwave = gsl_fft_complex_wavetable_alloc(samples);
	while ( (cwave->factor[ (cwave->nf) -1 ] > cwave->factor[ (cwave->nf)-2 ]) ) {
		gsl_fft_complex_wavetable_free(cwave);
		samples += SAMPLEROUND;
		cwave = gsl_fft_complex_wavetable_alloc(samples);
	}
	cspace = gsl_fft_complex_workspace_alloc(samples);
}
void PulseFreq::buildvectors(void){
	static gsl_complex z;
	//factorization();

	cvec = gsl_vector_complex_alloc(samples);
	rhovec = gsl_vector_alloc(samples);
	phivec = gsl_vector_calloc(samples); // calloc since we will start with flat zero phase for a pulse centered about vector index 0 and samples
	modamp = gsl_vector_alloc(samples);
	modphase = gsl_vector_calloc(samples);
	omega = new double[samples];
	time = new double[samples];
	omega[0] = 0.0;
	time[0] = 0.0;
	rhovec->data[0] = 0.0;
	z = gsl_complex_polar(rhovec->data[0],phivec->data[0]);
	gsl_vector_complex_set(cvec,0,z);
	omega[samples/2] = -(double)(samples/2)*domega;
	time[samples/2] = -(double)(samples/2)*dtime;
	rhovec->data[samples/2] = 0.0;
	z = gsl_complex_polar(rhovec->data[samples/2],phivec->data[samples/2]);
	gsl_vector_complex_set(cvec,samples/2,z);
	startind = (unsigned)((omega_center-(omega_width/2.0))/domega);
	stopind = (unsigned)((omega_center+(omega_width/2.0))/domega);
	onwidth = (unsigned)(omega_onwidth/domega); // 2.0)/domega);//  /10.0)/domega);// sin^2 goes from0..1 in 0..pi/2
	offwidth = (unsigned)(omega_offwidth/domega); // 2.0)/domega);//  /10.0)/domega);// sin^2 goes from0..1 in 0..pi/2
	for (unsigned i = 1; i<startind;i++){
		double i_dbl = double(i);
		omega[i] = domega*i_dbl;
		time[i] = dtime*i_dbl;
		rhovec->data[i] = 0.0;
		z = gsl_complex_polar(rhovec->data[i],phivec->data[i]);
		gsl_vector_complex_set(cvec,i,z);
		omega[samples-i] = -domega*i_dbl;
		time[samples-i] = -dtime*i_dbl;
		rhovec->data[samples-i] = 0.0;
		z = gsl_complex_polar(rhovec->data[samples-i],phivec->data[samples-i]);
		gsl_vector_complex_set(cvec,samples-i,z);
	}
	for (unsigned i = startind;i<startind+onwidth; i++){
		double i_dbl = double(i);
		omega[i] = domega*i_dbl;
		time[i] = dtime*i_dbl;
		rhovec->data[i] = rising(i);
		z = gsl_complex_polar(rhovec->data[i],phivec->data[i]);
		gsl_vector_complex_set(cvec,i,z);
		omega[samples-i] = -domega*i_dbl;
		time[samples-i] = -dtime*i_dbl;
		rhovec->data[samples-i] = rising(i);
		z = gsl_complex_polar(rhovec->data[samples-i],phivec->data[samples-i]);
		gsl_vector_complex_set(cvec,samples-i,z);
	}
	for (unsigned i = startind+onwidth;i<stopind-offwidth; i++){
		double i_dbl = double(i);
		omega[i] = domega*i_dbl;
		time[i] = dtime*i_dbl;
		rhovec->data[i] = 1.0;
		z = gsl_complex_polar(rhovec->data[i],phivec->data[i]);
		gsl_vector_complex_set(cvec,i,z);
		omega[samples-i] = -domega*i_dbl;
		time[samples-i] = -dtime*i_dbl;
		rhovec->data[samples-i] = 1.0;
		z = gsl_complex_polar(rhovec->data[samples-i],phivec->data[samples-i]);
		gsl_vector_complex_set(cvec,samples-i,z);
	}
	for (unsigned i = stopind-offwidth;i<stopind; i++){
		double i_dbl = double(i);
		omega[i] = domega*i_dbl;
		time[i] = dtime*i_dbl;
		rhovec->data[i] = falling(i);
		z = gsl_complex_polar(rhovec->data[i],phivec->data[i]);
		gsl_vector_complex_set(cvec,i,z);
		omega[samples-i] = -domega*i_dbl;
		time[samples-i] = -dtime*i_dbl;
		rhovec->data[samples-i] = falling(i);
		z = gsl_complex_polar(rhovec->data[samples-i],phivec->data[samples-i]);
		gsl_vector_complex_set(cvec,samples-i,z);
	}
	for (unsigned i = stopind;i<samples/2; i++){
		double i_dbl = double(i);
		omega[i] = domega*i_dbl;
		time[i] = dtime*i_dbl;
		rhovec->data[i] = 0.0;
		z = gsl_complex_polar(rhovec->data[i],phivec->data[i]);
		gsl_vector_complex_set(cvec,i,z);
		omega[samples-i] = -domega*i_dbl;
		time[samples-i] = -dtime*i_dbl;
		rhovec->data[samples-i] = 0.0;
		z = gsl_complex_polar(rhovec->data[samples-i],phivec->data[samples-i]);
		gsl_vector_complex_set(cvec,samples-i,z);
	}
	
}
void PulseFreq::killvectors(void){
	gsl_vector_free(modamp);
	gsl_vector_free(modphase);
	gsl_vector_free(rhovec);
	gsl_vector_free(phivec);
	gsl_vector_complex_free(cvec);
	gsl_fft_complex_wavetable_free(cwave);
	gsl_fft_complex_workspace_free(cspace);
	delete omega;
	delete time;
	omega = time = NULL;
}
void PulseFreq::killtheseonly(void){
	gsl_vector_free(modamp);
	gsl_vector_free(modphase);
	gsl_vector_free(rhovec);
	gsl_vector_free(phivec);
	gsl_vector_complex_free(cvec);
	gsl_fft_complex_workspace_free(cspace);
}

void PulseFreq::factorization(void){
        /* inspecting factorization of wavetable */
        clog << samples << " samples in " << cwave->nf << " factors:\t";
        for (unsigned i = 0;i<(cwave->nf);i++){
                clog << (cwave->factor)[i] << "\t";
        }
        clog << "\n";

}
