// Pulse clas implimentation
// standard includes
#include <cmath>


// my headers
#include <Pulse.hpp>
#include <boost/math/interpolators/barycentric_rational.hpp>
#include <fftw3.h>
#include <algorithm>
#include <Constants.hpp>
#include <boost/lexical_cast.hpp>
#include <DataOps.hpp>
#include <random>

using namespace Constants;
using namespace DataOps;

PulseFreq::PulseFreq(const double omcenter_in=(0.55*fsPau<float>()),const double omwidth_in=(0.15*fsPau<float>()),const double omonoff_in=(0.1*fsPau<float>()), double tspan_in=(10000.0/fsPau<float>())): // default constructor
	omega_center(omcenter_in),
	omega_width(omwidth_in ),
	omega_high( std::max(4.0*(omcenter_in + omwidth_in),10.0*omcenter_in) ),
	domega( 2.0*pi<double>()/tspan_in),
	omega(NULL),
	intime(false),
	infreq(true),
	time(NULL),
	parent(true),
	child(false),
	i_low( (unsigned)(double( atof( getenv("nu_low") ) )* 2.0*pi<double>()*fsPau<double>()/domega) ),
	i_high( (unsigned)(double( atof( getenv("nu_high") ) )* 2.0*pi<double>()*fsPau<double>()/domega) ),
	m_noisescale(1e-3),
	m_sampleinterval(2),
	m_saturate(4096),
	m_gain(1000000),
	m_lamsamples(2048)
{
	samples = (( (unsigned)(2.0 * omega_high / domega))/SAMPLEROUND + 1 ) *SAMPLEROUND;// dt ~ .1fs, Dt ~ 10000fs, dom = 2*pi/1e4, omhigh = 2*pi/.1fs, samples = omhigh/dom*2
	dtime = tspan_in/double(samples);
	omega_onwidth = omega_offwidth = omega_width/2.0; // forcing sin2 gaussian spectrum
	buildvectors();
	nu0=omcenter_in/(2.0*pi<double>())*fsPau<float>();
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
	i_low(rhs.i_low), 
	i_high(rhs.i_high), 
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


	nu0=rhs.nu0;
	buildvectors();

	DataOps::clone(rhovec,rhs.rhovec);
	DataOps::clone(phivec,rhs.phivec);
	DataOps::clone(cvec,rhs.cvec,samples);
	DataOps::clone(modamp,rhs.modamp);
	DataOps::clone(modphase,rhs.modphase);
}
PulseFreq::~PulseFreq(void){
	if(child){
		killtheseonly();
	} else {
		killvectors();
	}
}
void PulseFreq::rhophi2cvec(void)
{
	for (size_t i=0;i<samples;i++){
		cvec[i] = std::polar(rhovec[i],phivec[i]);
	}
}
void PulseFreq::cvec2rhophi(void)
{
	for (size_t i=0;i<samples;i++){
		rhovec[i] = std::abs(cvec[i]);
		phivec[i] = std::arg(cvec[i]);
	}
}

PulseFreq * PulseFreq::operator+=(const PulseFreq *rhs){*this += *rhs;return this;}
PulseFreq * PulseFreq::operator-=(const PulseFreq *rhs){*this -= *rhs;return this;}
PulseFreq * PulseFreq::operator*=(const PulseFreq *rhs){*this *= *rhs;return this;}
PulseFreq * PulseFreq::operator/=(const PulseFreq *rhs){*this /= *rhs;return this;}

PulseFreq & PulseFreq::operator+=(const PulseFreq &rhs){
	DataOps::sum(cvec,rhs.cvec,samples);
	cvec2rhophi();
	return *this;
}

PulseFreq & PulseFreq::operator-=(const PulseFreq &rhs){
	DataOps::diff(cvec,rhs.cvec,samples);
	cvec2rhophi();
        return *this;
}

PulseFreq & PulseFreq::operator*=(const PulseFreq &rhs){
	DataOps::mul(cvec,rhs.cvec,samples);
	cvec2rhophi();
	return *this;
}

PulseFreq & PulseFreq::operator*=(const double scale){
	DataOps::mul(cvec,scale,samples);
	cvec2rhophi();
	return *this;
}

PulseFreq & PulseFreq::operator/=(const PulseFreq &rhs){
	DataOps::div(cvec,rhs.cvec,samples);
	cvec2rhophi();
        return *this;
}

PulseFreq & PulseFreq::interfere(const PulseFreq &rhs){
	*this += rhs;
        return *this;
}

PulseFreq & PulseFreq::diffamps(const PulseFreq &rhs){
	rhovec -= rhs.rhovec;
	std::fill(phivec.begin(),phivec.end(),0);
	rhophi2cvec();
        return *this;
}

PulseFreq & PulseFreq::normamps(const PulseFreq &rhs){
	rhovec /= rhs.rhovec;
	std::fill(phivec.begin(),phivec.end(),0);
	rhophi2cvec();
        return *this;
}

void PulseFreq::print_amp(std::ofstream & outfile)
{
	outfile << "# amp\n";
	outfile << rhovec << std::endl;
}
void PulseFreq::print_phase(std::ofstream & outfile)
{
	outfile << "# phase\n";
	outfile << phivec << std::endl;
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
		std::cerr << "died here at addrandomphase()" << std::endl;
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
		phivec[i] += randphase[i]/samples;
	}
	rhophi2cvec();
	return true;
}


void PulseTime::setstrength(const double in)
{
  strength = in * auenergy<double>()/Eh<double>() * std::pow(aufor10PW<double>(),int(2));
}

void PulseTime::setwidth(const double in)
{
  Ctau = in * Constants::root_pi<double>()/ Constants::fsPau<double>() / 2.0;
}

void PulseTime::sett0(const double in)
{
  t0 = in / fsPau<double>();
}

void PulseFreq::attenuate(double attenfactor){
	rhovec *= attenfactor;
	rhophi2cvec();
}
void PulseFreq::phase(double phasein){ // expects delay in units of pi , i.e. 1.0 = pi phase flip 
	if(intime){
		fft_tofreq();
	}
	phivec += phasein*Constants::pi<double>();
	rhophi2cvec();
	if(intime){
		fft_totime();
	}
}
void PulseFreq::delay(double delayin){ // expects delay in fs
	if(intime){
		fft_tofreq();
	}
	for (unsigned i=0;i<samples;i++){
		phivec[i] += omega[i]*delayin/fsPau<double>();
	}
	rhophi2cvec();
	if(intime){
		fft_totime();
	}
}


void PulseFreq::printfrequency(std::ofstream * outfile){
	double nu,lambda;
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
		lambda = C_nmPfs<double>()/nu;
		(*outfile) << nu << "\t" << lambda << "\t" << rhovec->data[i] << "\t" << phivec->data[i] << "\n";
	}
}
void PulseFreq::printwavelengthbins(std::ofstream * outfile)
{
	std::vector<double> x(2);
	x.front() = C_nmPfs<double>()*2.0*pi<double>()*fsPau<double>()/omega[i_low];
	x.back() = C_nmPfs<double>()*2.0*pi<double>()*fsPau<double>()/omega[i_high-1];
	double dlam = (x.front()-x.back())/double(m_lamsamples);
        for (size_t i = 0;i<m_lamsamples;++i){
		(*outfile) << x.back() + i*dlam << "\t";
        }
	(*outfile) << "\n";
	return;
}
void PulseFreq::appendwavelength(std::ofstream * outfile)
{
	std::vector<double> x(i_high-i_low);
	std::vector<double> y(i_high-i_low);	
	for (size_t i=0;i<y.size();++i){
		x[i] = C_nmPfs<double>()*2.0*pi<double>()*fsPau<double>()/omega[i_low+i];
		y[i] = std::min(rhovec->data[i_low+i] * m_gain,double(m_saturate));
	}
	double dlam = (x.front()-x.back())/double(m_lamsamples);
	boost::math::barycentric_rational<double> interpolant(x.data(), y.data(), y.size());
	for (size_t i=0;i<m_lamsamples;++i){
		(*outfile) << uint16_t(interpolant(x.back()+i*dlam)) << "\t";
	}
	(*outfile) << std::endl;
	return;
}
void PulseFreq::appendfrequency(std::ofstream * outfile){
        for (unsigned i = i_low;i<i_high;i+=m_sampleinterval){ 
		uint16_t val = std::min(unsigned(rhovec->data[i] * m_gain),unsigned(m_saturate));
       		(*outfile) << val << "\t";
        }
	(*outfile) << std::endl;
}

void PulseFreq::appendnoisy(std::ofstream * outfile){
	std::normal_distribution<double> norm_dist( 0.0, m_noisescale);
        for (unsigned i = i_low;i<i_high;i+=m_sampleinterval){ 
		double outval = rhovec[i] + norm_dist(rng);
       		(*outfile) << outval << "\t" ;
        }
	(*outfile) << std::endl;
}

void PulseFreq::printfrequencybins(std::ofstream * outfile){
	double nu,lam;
        for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
		lam = C_nmPfs<double>()/nu;
       		(*outfile) << nu << "\t" << lam << "\n";
        }
	(*outfile) << "\n";
}
void PulseFreq::appendfrequencybins(std::ofstream * outfile){
	double nu,lam;
        for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
       		(*outfile) << nu << "\t";
        }
	(*outfile) << "\n";
}
void PulseFreq::printfrequencydelay(std::ofstream * outfile, const double *delay){
	double nu,lambda,thisdelay;
	thisdelay = (*delay)*fsPau<double>();
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\tdelay[fs]\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
		lambda = C_nmPfs<double>()/nu;
		(*outfile) << nu << "\t" << lambda << "\t" << rhovec->data[i] << "\t" << phivec->data[i] << "\t" << thisdelay << "\n";
	}
	(*outfile) << "\n";
}
void PulseFreq::printfrequencydelaychirp(std::ofstream * outfile, const double *delay,const double *chirp){
	double nu,lambda,thisdelay,thischirp,reltime;
	thisdelay = (*delay)*fsPau<double>();
	thischirp = (*chirp)*std::pow(fsPau<double>(),int(2));
	(*outfile) << ("#omega[fs^-1]\tlambda[nm]\trho\tphi\tdelay[fs]\tchirp[fs^2]\treldelays[fs]\n");
	for (unsigned i = i_low;i<i_high;i+=(unsigned)(atoi(getenv("sampleinterval")))){
		nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
		lambda = C_nmPfs<double>()/nu;
		reltime = (nu-nu0)*thischirp*10.0;
		(*outfile) << nu << "\t" << lambda << "\t" << rhovec->data[i] << "\t" << phivec->data[i] << "\t" << thisdelay << "\t" << thischirp << "\t" << reltime << "\n";
	}
	(*outfile) << "\n";
}

void PulseFreq::printwavelength(std::ofstream * outfile,const double *delay){
	double nu,lambda,lamlast;
	lamlast=0;
        for (unsigned i = i_high;i>i_low;i-=(unsigned)(atoi(getenv("sampleinterval")))){
                nu = (omega[i]/(2.0*pi<double>())/fsPau<double>());
                lambda = (double)((int)((C_nmPfs<double>()/nu)*10.0))/10.0;
		if (lambda>lamlast & (int)(lambda*10)%5 == 0){
                	(*outfile) << lambda << "\t" << (*delay) << "\t" << rhovec->data[i] << "\n";
			lamlast = lambda;
		}
        }
	(*outfile) << "\n";
}

void PulseFreq::printtime(std::ofstream * outfile){
	(*outfile) << ("#time\treal\timag\n");
	for (unsigned i = samples/2;i<samples;i++){
		(*outfile) << (time[i]*fsPau<double>()) << "\t" << cvec[i].real() << "\t" << cvec[i].imag() << "\n";
	}
	for (unsigned i = 0;i<samples/2;i++){
		(*outfile) << (time[i]*fsPau<double>()) << "\t" << cvec[i].real() << "\t" << cvec[i].imag() << "\n";
	}

} 


void PulseFreq::buildvectors(void){
		
	if (cvec != NULL){
		delete cvec;
		cvec = NULL;
	}
	cvec = (std::complex<double> *) fftw_malloc(sizeof(std::complex<double>) * samples);
	std::fill(cvec,samples,std::complex<double>(0));

	FTplan_forward = fftw_plan_dft_1d(samples, cvec, cvec, FFTW_FORWARD, FFTW_ESTIMATE);
	FTplan_backward = fftw_plan_dft_1d(samples, cvec, cvec, FFTW_BACKWARD, FFTW_ESTIMATE);

	static std::complex<double> z;
	//factorization();

	rhovec.resize(samples,0.0);
	phivec.resize(samples,0.0);
	modamp.resize(samples,1.0);
	modphase.resize(samples,0.0);
	omega.resize(samples);
	time.resize(samples);
	//omega = new double[samples];
	//time = new double[samples];
	omega[0] = 0.0;
	time[0] = 0.0;
	omega[samples/2] = -(double)(samples/2)*domega;
	time[samples/2] = -(double)(samples/2)*dtime;

	startind = (unsigned)((omega_center-(omega_width/2.0))/domega);
	stopind = (unsigned)((omega_center+(omega_width/2.0))/domega);
	onwidth = (unsigned)(omega_onwidth/domega); // 2.0)/domega);//  /10.0)/domega);// sin^2 goes from0..1 in 0..pi/2
	offwidth = (unsigned)(omega_offwidth/domega); // 2.0)/domega);//  /10.0)/domega);// sin^2 goes from0..1 in 0..pi/2
	for (unsigned i = 1; i<startind;i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		omega[samples-i] = -domega*i;
		time[samples-i] = -dtime*i;
	}
	for (unsigned i = startind;i<startind+onwidth; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		rhovec[i] = rising(i);
		cvec[i] = std::polar(rhovec[i],phivec[i]);
		omega[samples-i] = -domega*i_dbl;
		time[samples-i] = -dtime*i_dbl;
		rhovec[samples-i] = rising(i);
		cvec[samples-i] = std::polar(rhovec[samples-i],phivec[samples-i]);
	}
	for (unsigned i = startind+onwidth;i<stopind-offwidth; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		rhovec->data[i] = 1.0;
		cvec[i] = std::polar(rhovec[i],phivec[i]);
		omega[samples-i] = -domega*i_dbl;
		time[samples-i] = -dtime*i_dbl;
		rhovec[samples-i] = 1.0;
		cvec[samples-i] = std::polar(rhovec[samples-i],phivec[samples-i]);
	}
	for (unsigned i = stopind-offwidth;i<stopind; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		rhovec[i] = falling(i);
		cvec[i] = std::polar(rhovec[i],phivec[i]);
		omega[samples-i] = -domega*i;
		time[samples-i] = -dtime*i;
		rhovec->data[samples-i] = falling(i);
		cvec[samples - i] = std::polar(rhovec[samples - i],phivec[samples - i]);
	}
	for (unsigned i = stopind;i<samples/2; i++){
		omega[i] = domega*i;
		time[i] = dtime*i;
		omega[samples-i] = -domega*i;
		time[samples-i] = -dtime*i;
	}
	
}
void PulseFreq::killvectors(void){
	fftw_destroy_plan(FTplan);
	fftw_free(cvec);
	delete omega;
	delete time;
	omega = time = NULL;
}
void PulseFreq::killtheseonly(void){
	fftw_destroy_plan(FTplan);
	fftw_free(cvec);
}

