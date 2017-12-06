#include "scan_material.hpp"
#include <cstdlib>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <complex>

#include <vector>
#include <random>
#include <chrono>

#include <DebugOps.hpp>

using namespace Constants;


/* Here's main */
int main(int argc, char* argv[])
{
	/* start reading in from the command line, make and write fake data */

	size_t nimages = size_t(atoi(getenv("nimages")));

	/*** Spawn a parallel region explicitly scoping all variables ***/
	unsigned nthreads = (unsigned)atoi( getenv("nthreads") );

	//std::default_random_engine rng;
	std::random_device rng;
	std::normal_distribution<double> distribution(
			double(atof(getenv("delays_mean"))),
			double(atof(getenv("delays_std"))));
	std::uniform_real_distribution<double> uni_distribution(
			double(atof(getenv("delays_mean")))-double(atof(getenv("delays_std"))),
			double(atof(getenv("delays_mean")))+double(atof(getenv("delays_std"))));

	std::normal_distribution<float> xray_pos_distribution(
			float(atof(getenv("xray_pos_distribution_normal_mean"))),
			float(atof(getenv("xray_pos_distribution_normal_std"))));
	std::normal_distribution<float> laser_pos_distribution(
			float(atof(getenv("laser_pos_distribution_normal_mean"))),
			float(atof(getenv("laser_pos_distribution_normal_std"))));
	std::lognormal_distribution<float> xray_distribution(
			float(atof(getenv("xray_distribution_lognormal_mean"))),
			float(atof(getenv("xray_distribution_lognormal_std"))));
	std::lognormal_distribution<float> laser_distribution(
			float(atof(getenv("laser_distribution_lognormal_mean"))),
			float(atof(getenv("laser_distribution_lognormal_std"))));


	std::string filename;
	std::string filebase = std::string(getenv("filebase"));



	float dalpha(float(atof(getenv("drifting_alpha")))*pi<double>()/nimages);


	bool doublepulse = ((int)atoi(getenv("doublepulse")) > 0 ? true : false);
	double doublepulsedelay(0.);
	if (doublepulse){
		doublepulsedelay = (double)atof( getenv("doublepulsedelay") ) ; // this one gets used directly in atomic units
	}
	double lambda0 = (double)( atof( getenv("lambda0") ) );
	double lambda_width = (double)( atof( getenv("lambda_width") ) );
	double lambda_onoff = (double)( atof( getenv("lambda_onoff") ) );
	double tspan = (double)( atof( getenv("tspan") ) )/fsPau<float>();
	double omega_low = twopi<double>()/(lambda0+lambda_width/2.0)*C_nmPfs<double>()*fsPau<double>();
	double omega_high = twopi<double>()/(lambda0-lambda_width/2.0)*C_nmPfs<double>()*fsPau<double>();
	double omega_width = omega_high-omega_low;
	double omega_onoff = twopi<double>()/(lambda0+lambda_width/2.0-lambda_onoff)*C_nmPfs<double>()*fsPau<double>() - omega_low;

	double omega0 = (omega_high+omega_low)/2.0;

	unsigned ngroupsteps = (unsigned)( atoi( getenv("ngroupsteps") ) );
	double groupdelay = (double)(atof(getenv("groupdelay"))); //this one gets used in fs
	double backdelay = (double)(atof(getenv("backdelay")));
	double groupstep = groupdelay/ngroupsteps ;
	double backstep = backdelay/ngroupsteps;
	const unsigned netalon = (unsigned)(atoi(getenv("netalon")));

	PulseFreq masterpulse(omega0,omega_width,omega_onoff,tspan);

	if (atoi(getenv("addrandomphase"))>0){
		masterpulse.addrandomphase();
		std::string filename = filebase + "spectralphaseFTpower.dat";
		std::ofstream outfile(filename,std::ios::out);
		masterpulse.print_phase_powerspectrum(outfile);
		outfile.close();
		filename = filebase + "spectralphase.dat";
		outfile.open(filename);
		masterpulse.print_phase(outfile);
		outfile.close();
		filename = filebase + "spectralamp.dat";
		outfile.open(filename);
		masterpulse.print_amp(outfile);
		outfile.close();
	}

	std::vector<double> chirpvec(4,0.);
	std::vector<double> chirpnoisevec(4,0.);
	chirpvec[0] = (double)( atof( getenv("chirp") ) ) / std::pow(fsPau<float>(),int(2));// the difference in slopes at omega_low versus omega_high must equal tspan
	chirpvec[1] = (double)( atof( getenv("TOD") ) ) / std::pow(fsPau<float>(),int(3));
	chirpvec[2] = (double)( atof( getenv("FOD") ) ) / std::pow(fsPau<float>(),int(4));
	chirpvec[3] = (double)( atof( getenv("fifthOD") ) ) / std::pow(fsPau<float>(),int(5));

	masterpulse.addchirp(chirpvec);							// chirp that ref pulse


	FiberBundle masterbundle(boost::lexical_cast<size_t>(atoi(getenv("nfibers"))));
	masterbundle.fiberdiameter(boost::lexical_cast<float>(atof(getenv("fiberdiam"))));
	masterbundle.laserdiameter(boost::lexical_cast<float>(atof(getenv("laserdiam"))));
	masterbundle.xraydiameter(boost::lexical_cast<float>(atof(getenv("xraydiam"))));
	masterbundle.set_fsPmm(boost::lexical_cast<float>(atof(getenv("bundle_fsPmm"))));
	masterbundle.scalePolarCoords();
	masterbundle.shuffle_output();
	masterbundle.Ixray(float(1.));
	masterbundle.Ilaser(float(1.));
	filename = filebase + "fibermap.out";
	std::cerr << filename << std::endl;
	std::ofstream mapfile(filename.c_str(),std::ios::out);
	masterbundle.print_mapping(mapfile);
	mapfile.close();
	MatResponse masterresponse(
			0,															// stepdelay
			boost::lexical_cast<double>( atof( getenv("stepwidth") ) ),								// stepwidth
			(boost::lexical_cast<double>( atof( getenv("attenuation") ) ) - 1.0) * masterbundle.Ixray() / ngroupsteps + 1.0,	// attenuation
			boost::lexical_cast<double>( atof( getenv("phase") ) ) * masterbundle.Ixray() / ngroupsteps				// phase
			);
	masterresponse.aalphabbeta(
			(double)( atof( getenv("a") ) ),		// a
			(double)( atof( getenv("alpha" ) ) ),		// alpha
			(double)( atof( getenv("b") ) ),		// b
			(double)( atof( getenv("beta") ) )		// beta
			);
	masterresponse.setreflectance((double)(atof(getenv("etalon"))));
	masterresponse.setetalondelay((double)(atof(getenv("etalondelay"))));



#pragma omp parallel num_threads(nthreads) shared(masterbundle,masterpulse)
	{ // begin parallel region
		size_t tid = omp_get_thread_num();
		FiberBundle parabundle(masterbundle);
		std::vector< PulseFreq* > pulsearray(parabundle.get_nfibers());
		std::vector< PulseFreq* > crosspulsearray(parabundle.get_nfibers());

//#pragma omp for // ordered 
		for (size_t n=0;n<nimages;++n){ // outermost loop for nimages to produce //
			double t0 = uni_distribution(rng);
			double startdelay;
			std::cout << "Running for t0 = " << t0 << std::endl;
			DebugOps::pushout(std::string("in threaded for loop, thread" + std::to_string(tid)));

			parabundle.Ixray(float(xray_distribution(rng)));
			parabundle.Ilaser(float(laser_distribution(rng)));
			parabundle.delay_angle(dalpha*float(n));
			parabundle.center_Ixray(xray_pos_distribution(rng),xray_pos_distribution(rng));
			parabundle.center_Ilaser(laser_pos_distribution(rng),laser_pos_distribution(rng));

			PulseFreq * pulsearray_i = new PulseFreq(masterpulse);
			PulseFreq * crosspulsearray_i = new PulseFreq(masterpulse);
			MatResponse pararesponse(masterresponse);

			for(size_t i = 0; i< parabundle.get_nfibers(); i++){ // begin fibers loop

				*pulsearray_i = masterpulse;
				*crosspulsearray_i = masterpulse;


				std::normal_distribution<double> chirpnoiseDist( 
						(double)( atof( getenv("chirp") ) ) / std::pow(fsPau<float>(),(int)2), 
						(double)( atof( getenv("chirpnoise") ) ) / std::pow(fsPau<float>(),(int)2) );
				std::normal_distribution<double> TODnoiseDist( 
						(double)( atof( getenv("TOD") ) ) / std::pow(fsPau<float>(),(int)3),
						(double)( atof( getenv("TODnoise") ) ) / std::pow(fsPau<float>(),(int)3) );
				std::normal_distribution<double> FODnoiseDist(
						(double)( atof( getenv("FOD") ) ) / std::pow(fsPau<float>(),(int)4),
						(double)( atof( getenv("FODnoise") ) ) / std::pow(fsPau<float>(),(int)4) );
				std::normal_distribution<double> fifthODnoiseDist(
						(double)( atof( getenv("fifthOD") ) ) / std::pow(fsPau<float>(),(int)5),
						(double)( atof( getenv("fifthODnoise") ) ) / std::pow(fsPau<float>(),(int)5) );
				double interferedelay = (double)atof(getenv("interferedelay"));

				startdelay = parabundle.delay(i) + t0;
				pararesponse.set_delay(startdelay);
				pararesponse.set_attenuation( (double( atof( getenv("attenuation") ) ) - 1.0) * parabundle.Ixray(size_t(i)) / ngroupsteps + 1.0 ); 
				pararesponse.set_phase(double( atof( getenv("phase") ) ) * parabundle.Ixray(size_t(i)) / ngroupsteps);

				std::vector<double> chirpvec(4,0.);
				chirpvec[0] = chirpnoiseDist(rng);
				chirpvec[1] = TODnoiseDist(rng);
				chirpvec[2] = FODnoiseDist(rng);
				chirpvec[3] = fifthODnoiseDist(rng);
				pulsearray_i->addchirp(chirpvec); 
				crosspulsearray_i->addchirp(chirpvec); 
				crosspulsearray_i->delay(interferedelay); // delay in the frequency domain
				pulsearray_i->fft_totime();
				crosspulsearray_i->fft_totime();

				DebugOps::pushout("\nhere before oneresponsePtr setting in tid ",tid);


				for(size_t g=0;g<ngroupsteps;g++){ // begin groupsteps loop
					pararesponse.setdelay(startdelay - g*groupstep); // forward propagating, x-rays advance on the optical
					pararesponse.setstepvec_amp(pulsearray_i);
					pararesponse.setstepvec_phase(pulsearray_i);
					pararesponse.setstepvec_amp(crosspulsearray_i);
					pararesponse.setstepvec_phase(crosspulsearray_i);
					if (doublepulse){
						pararesponse.addstepvec_amp(pulsearray_i,doublepulsedelay);
						pararesponse.addstepvec_phase(pulsearray_i,doublepulsedelay);
						pararesponse.addstepvec_amp(crosspulsearray_i,doublepulsedelay);
						pararesponse.addstepvec_phase(crosspulsearray_i,doublepulsedelay);
					}
					// this pulls down the tail of the response so vector is periodic on nsamples	
					pararesponse.buffervectors(pulsearray_i); 
					pararesponse.buffervectors(crosspulsearray_i); 
					pulsearray_i->modulateamp_time();
					pulsearray_i->modulatephase_time();
					crosspulsearray_i->modulateamp_time();
					crosspulsearray_i->modulatephase_time();
				}// end groupsteps loop

				for (size_t e=0;e<netalon;e++){ // begin etalon loop
					// back propagation step //
					double etalondelay = startdelay - double(e+1) * (pararesponse.getetalondelay()); // at front surface, x-rays see counter-propagating light from one full etalon delay
					PulseFreq etalonpulse(*(pulsearray_i));
					PulseFreq crossetalonpulse(*(crosspulsearray_i));

					for(size_t g=0;g<ngroupsteps;g++){
						pararesponse.setdelay(etalondelay + g*backstep); // counterpropagating, x-rays work backwards through the optical
						pararesponse.setstepvec_amp(&etalonpulse);
						pararesponse.setstepvec_phase(&etalonpulse);
						pararesponse.setstepvec_amp(&crossetalonpulse);
						pararesponse.setstepvec_phase(&crossetalonpulse);
						if (doublepulse){
							pararesponse.addstepvec_amp(&etalonpulse,doublepulsedelay);
							pararesponse.addstepvec_phase(&etalonpulse,doublepulsedelay);
							pararesponse.addstepvec_amp(&crossetalonpulse,doublepulsedelay);
							pararesponse.addstepvec_phase(&crossetalonpulse,doublepulsedelay);
						}
						pararesponse.buffervectors(&etalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						pararesponse.buffervectors(&crossetalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						etalonpulse.modulateamp_time();
						etalonpulse.modulatephase_time();
						crossetalonpulse.modulateamp_time();
						crossetalonpulse.modulatephase_time();
					}
					// forward propagation //
					for(size_t g=0;g<ngroupsteps;g++){
						pararesponse.setdelay(startdelay - g*groupstep); // forward propagating, x-rays advance on the optical
						pararesponse.setstepvec_amp(&etalonpulse);
						pararesponse.setstepvec_phase(&etalonpulse);
						pararesponse.setstepvec_amp(&crossetalonpulse);
						pararesponse.setstepvec_phase(&crossetalonpulse);
						if (doublepulse){
							pararesponse.addstepvec_amp(&etalonpulse,doublepulsedelay);
							pararesponse.addstepvec_phase(&etalonpulse,doublepulsedelay);
							pararesponse.addstepvec_amp(&crossetalonpulse,doublepulsedelay);
							pararesponse.addstepvec_phase(&crossetalonpulse,doublepulsedelay);

						}
						pararesponse.buffervectors(&etalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						pararesponse.buffervectors(&crossetalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						etalonpulse.modulateamp_time();
						etalonpulse.modulatephase_time();
						crossetalonpulse.modulateamp_time();
						crossetalonpulse.modulatephase_time();
					}
					etalonpulse.fft_tofreq();
					crossetalonpulse.fft_tofreq();
					etalonpulse.delay(pararesponse.getetalondelay()); // delay and attenuate in frequency domain
					etalonpulse.attenuate(pow(pararesponse.getreflectance(),(int)2));
					crossetalonpulse.delay(pararesponse.getetalondelay()); // delay and attenuate in frequency domain
					crossetalonpulse.attenuate(pow(pararesponse.getreflectance(),(int)2));
					etalonpulse.fft_totime();
					crossetalonpulse.fft_totime();
					*(pulsearray_i) += etalonpulse;
					*(crosspulsearray_i) += crossetalonpulse;

				} // end etalon loop

				pulsearray_i->fft_tofreq();
				crosspulsearray_i->fft_tofreq();

				crosspulsearray_i->delay(-interferedelay); // expects this in fs // time this back up to the crosspulse

				pulsearray[i] = new PulseFreq(*pulsearray_i);
				crosspulsearray[i] = new PulseFreq(*crosspulsearray_i);

			} // end nfibers loop




			filename = filebase + "interference.out.thread" + std::to_string(tid) + "." + std::to_string(n);
			ofstream interferestream(filename.c_str(),ios::out); // use app to append delays to same file.
			std::cout << "interfere filename out = " << filename << std::endl;
			std::complex<float> z_laser = parabundle.center_Ilaser();
			std::complex<float> z_xray = parabundle.center_Ixray();
			interferestream << "#delay for image = \t" << t0 
				<< "\n#Ilaser = \t" << parabundle.Ilaser()
				<< "\n#Ixray = \t" << parabundle.Ixray()
				<< "\n#center laser = \t" << z_laser.real() << "\t" << z_laser.imag() 
				<< "\n#center xray = \t" << z_xray.real() << "\t" << z_xray.imag()
				<< "\n#alpha = \t" << parabundle.delay_angle() 
				<< std::endl;
			interferestream << "#";
			pulsearray[0]->printwavelengthbins(&interferestream);

			for (size_t i=0;i<parabundle.get_nfibers();i++){
				*pulsearray[i] -= *crosspulsearray[i];
				*pulsearray[i] *= parabundle.Ilaser(size_t(i)); 
				pulsearray[i]->appendwavelength(&interferestream);
			}
			interferestream.close();

			DebugOps::pushout(std::string("At end of images loop in thread" + std::to_string(tid) + "\t"),n);

		} // outermost loop for nimages to produce //


	// file for delay bins
	filename = filebase + "delaybins.out.thread" + std::to_string(tid);
	ofstream outbins(filename.c_str(),ios::out); 
	for (size_t i=0;i<parabundle.get_nfibers();++i){
		outbins << parabundle.delay(i) << "\n";
	}
	outbins.close();
		for (size_t i=0;i<pulsearray.size();++i){
			delete pulsearray[i];
			delete crosspulsearray[i];
		}
	} // end parallel region

	return 0;
}


