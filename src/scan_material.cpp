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
	PulseFreq * masterpulsePtr = &masterpulse;

	if (atoi(getenv("addrandomphase"))>0){
		masterpulsePtr->addrandomphase();
		std::string filename = filebase + "spectralphaseFTpower.dat";
		std::ofstream outfile(filename,std::ios::out);
		masterpulsePtr->print_phase_powerspectrum(outfile);
		outfile.close();
		filename = filebase + "spectralphase.dat";
		outfile.open(filename);
		masterpulsePtr->print_phase(outfile);
		outfile.close();
		filename = filebase + "spectralamp.dat";
		outfile.open(filename);
		masterpulsePtr->print_amp(outfile);
		outfile.close();
	}

	FiberBundle bundle(boost::lexical_cast<size_t>(atoi(getenv("nfibers"))));
	bundle.fiberdiameter(boost::lexical_cast<float>(atof(getenv("fiberdiam"))));
	bundle.laserdiameter(boost::lexical_cast<float>(atof(getenv("laserdiam"))));
	bundle.xraydiameter(boost::lexical_cast<float>(atof(getenv("xraydiam"))));
	bundle.set_fsPmm(boost::lexical_cast<float>(atof(getenv("bundle_fsPmm"))));
	bundle.scalePolarCoords();
	bundle.shuffle_output();
	bundle.Ixray(float(1.));
	bundle.Ilaser(float(1.));
	filename = filebase + "fibermap.out";
	std::cerr << filename << std::endl;
	std::ofstream mapfile(filename.c_str(),std::ios::out);
	bundle.print_mapping(mapfile);
	mapfile.close();

	FiberBundle * bundlePtr = new FiberBundle(bundle);



	std::vector<double> chirpvec(4,0.);
	std::vector<double> chirpnoisevec(4,0.);
	chirpvec[0] = (double)( atof( getenv("chirp") ) ) / std::pow(fsPau<float>(),int(2));// the difference in slopes at omega_low versus omega_high must equal tspan
	chirpvec[1] = (double)( atof( getenv("TOD") ) ) / std::pow(fsPau<float>(),int(3));
	chirpvec[2] = (double)( atof( getenv("FOD") ) ) / std::pow(fsPau<float>(),int(4));
	chirpvec[3] = (double)( atof( getenv("fifthOD") ) ) / std::pow(fsPau<float>(),int(5));

	masterpulse.addchirp(chirpvec);							// chirp that ref pulse

#pragma omp parallel  num_threads(nthreads)
	{ // begin parallel region


	// =========== Definately, make the pulses outside of the for loop, but one per thread //
	std::vector< PulseFreq* > pulsearray(bundlePtr->get_nfibers());
	std::vector< PulseFreq* > crosspulsearray(bundlePtr->get_nfibers());

	for (size_t i=0;i<pulsearray.size();++i){
		pulsearray[i] = new PulseFreq(*masterpulsePtr);
		crosspulsearray[i] = new PulseFreq(*masterpulsePtr);
	}


		size_t tid = omp_get_thread_num();

#pragma omp for // ordered 
		for (size_t n=0;n<nimages;++n){ // outermost loop for nimages to produce //
			double t0 = uni_distribution(rng);
			double startdelay;
			std::cout << "Running for t0 = " << t0 << std::endl;

			bundle.Ixray(float(xray_distribution(rng)));
			bundle.Ilaser(float(laser_distribution(rng)));
			bundle.delay_angle(dalpha*float(n));
			bundle.center_Ixray(xray_pos_distribution(rng),xray_pos_distribution(rng));
			bundle.center_Ilaser(laser_pos_distribution(rng),laser_pos_distribution(rng));

			for(size_t i = 0; i< bundle.get_nfibers(); i++){ // begin fibers loop
				DebugOps::pushout(std::string("in threaded for loop, thread" + std::to_string(tid)));
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

				startdelay = bundle.delay(i) + t0;

				std::vector<double> chirpvec(4,0.);
				chirpvec[0] = chirpnoiseDist(rng);
				chirpvec[1] = TODnoiseDist(rng);
				chirpvec[2] = FODnoiseDist(rng);
				chirpvec[3] = fifthODnoiseDist(rng);
				pulsearray[i]->addchirp(chirpvec); 
				crosspulsearray[i]->addchirp(chirpvec); 
				crosspulsearray[i]->delay(interferedelay); // delay in the frequency domain
				pulsearray[i]->fft_totime();
				crosspulsearray[i]->fft_totime();
				MatResponse * oneresponsePtr = new MatResponse(
						startdelay,											// stepdelay
						boost::lexical_cast<double>( atof( getenv("stepwidth") ) ),					// stepwidth
						(boost::lexical_cast<double>( atof( getenv("attenuation") ) ) - 1.0) * bundle.Ixray(size_t(i)) / ngroupsteps + 1.0,// attenuation
						boost::lexical_cast<double>( atof( getenv("phase") ) ) * bundle.Ixray(size_t(i)) / ngroupsteps				// phase
						);
				oneresponsePtr->aalphabbeta(
						(double)( atof( getenv("a") ) ),		// a
						(double)( atof( getenv("alpha" ) ) ),		// alpha
						(double)( atof( getenv("b") ) ),		// b
						(double)( atof( getenv("beta") ) )		// beta
						);
				oneresponsePtr->setreflectance((double)(atof(getenv("etalon"))));
				oneresponsePtr->setetalondelay((double)(atof(getenv("etalondelay"))));

				for(size_t g=0;g<ngroupsteps;g++){ // begin groupsteps loop
					oneresponsePtr->setdelay(startdelay - g*groupstep); // forward propagating, x-rays advance on the optical
					oneresponsePtr->setstepvec_amp(pulsearray[i]);
					oneresponsePtr->setstepvec_phase(pulsearray[i]);
					oneresponsePtr->setstepvec_amp(crosspulsearray[i]);
					oneresponsePtr->setstepvec_phase(crosspulsearray[i]);
					if (doublepulse){
						oneresponsePtr->addstepvec_amp(pulsearray[i],doublepulsedelay);
						oneresponsePtr->addstepvec_phase(pulsearray[i],doublepulsedelay);
						oneresponsePtr->addstepvec_amp(crosspulsearray[i],doublepulsedelay);
						oneresponsePtr->addstepvec_phase(crosspulsearray[i],doublepulsedelay);
					}
					// this pulls down the tail of the response so vector is periodic on nsamples	
					oneresponsePtr->buffervectors(pulsearray[i]); 
					oneresponsePtr->buffervectors(crosspulsearray[i]); 
					pulsearray[i]->modulateamp_time();
					pulsearray[i]->modulatephase_time();
					crosspulsearray[i]->modulateamp_time();
					crosspulsearray[i]->modulatephase_time();
				}// end groupsteps loop

				for (size_t e=0;e<netalon;e++){ // begin etalon loop
					// back propagation step //
					double etalondelay = startdelay - double(e+1) * (oneresponsePtr->getetalondelay()); // at front surface, x-rays see counter-propagating light from one full etalon delay
					PulseFreq etalonpulse(*(pulsearray[i]));
					PulseFreq crossetalonpulse(*(crosspulsearray[i]));

					for(size_t g=0;g<ngroupsteps;g++){
						oneresponsePtr->setdelay(etalondelay + g*backstep); // counterpropagating, x-rays work backwards through the optical
						oneresponsePtr->setstepvec_amp(&etalonpulse);
						oneresponsePtr->setstepvec_phase(&etalonpulse);
						oneresponsePtr->setstepvec_amp(&crossetalonpulse);
						oneresponsePtr->setstepvec_phase(&crossetalonpulse);
						if (doublepulse){
							oneresponsePtr->addstepvec_amp(&etalonpulse,doublepulsedelay);
							oneresponsePtr->addstepvec_phase(&etalonpulse,doublepulsedelay);
							oneresponsePtr->addstepvec_amp(&crossetalonpulse,doublepulsedelay);
							oneresponsePtr->addstepvec_phase(&crossetalonpulse,doublepulsedelay);
						}
						oneresponsePtr->buffervectors(&etalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						oneresponsePtr->buffervectors(&crossetalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						etalonpulse.modulateamp_time();
						etalonpulse.modulatephase_time();
						crossetalonpulse.modulateamp_time();
						crossetalonpulse.modulatephase_time();
					}
					// forward propagation //
					for(size_t g=0;g<ngroupsteps;g++){
						oneresponsePtr->setdelay(startdelay - g*groupstep); // forward propagating, x-rays advance on the optical
						oneresponsePtr->setstepvec_amp(&etalonpulse);
						oneresponsePtr->setstepvec_phase(&etalonpulse);
						oneresponsePtr->setstepvec_amp(&crossetalonpulse);
						oneresponsePtr->setstepvec_phase(&crossetalonpulse);
						if (doublepulse){
							oneresponsePtr->addstepvec_amp(&etalonpulse,doublepulsedelay);
							oneresponsePtr->addstepvec_phase(&etalonpulse,doublepulsedelay);
							oneresponsePtr->addstepvec_amp(&crossetalonpulse,doublepulsedelay);
							oneresponsePtr->addstepvec_phase(&crossetalonpulse,doublepulsedelay);

						}
						oneresponsePtr->buffervectors(&etalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						oneresponsePtr->buffervectors(&crossetalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						etalonpulse.modulateamp_time();
						etalonpulse.modulatephase_time();
						crossetalonpulse.modulateamp_time();
						crossetalonpulse.modulatephase_time();
					}
					etalonpulse.fft_tofreq();
					crossetalonpulse.fft_tofreq();
					etalonpulse.delay(oneresponsePtr->getetalondelay()); // delay and attenuate in frequency domain
					etalonpulse.attenuate(pow(oneresponsePtr->getreflectance(),(int)2));
					crossetalonpulse.delay(oneresponsePtr->getetalondelay()); // delay and attenuate in frequency domain
					crossetalonpulse.attenuate(pow(oneresponsePtr->getreflectance(),(int)2));
					etalonpulse.fft_totime();
					crossetalonpulse.fft_totime();
					*(pulsearray[i]) += etalonpulse;
					*(crosspulsearray[i]) += crossetalonpulse;

				} // end etalon loop

				pulsearray[i]->fft_tofreq();
				crosspulsearray[i]->fft_tofreq();

				crosspulsearray[i]->delay(-interferedelay); // expects this in fs // time this back up to the crosspulse
				DebugOps::pushout(std::string("Before delete oneresponsePtr loop in thread" + std::to_string(tid) + "\t"),i);

				delete oneresponsePtr;

				DebugOps::pushout(std::string("At end of nfibers loop in thread" + std::to_string(tid) + "\t"),i);

			} // end nfibers loop




			filename = filebase + "interference.out.thread" + std::to_string(tid) + "." + std::to_string(n);
			ofstream interferestream(filename.c_str(),ios::out); // use app to append delays to same file.
			std::cout << "interfere filename out = " << filename << std::endl;
			std::complex<float> z_laser = bundle.center_Ilaser();
			std::complex<float> z_xray = bundle.center_Ixray();
			interferestream << "#delay for image = \t" << t0 
				<< "\n#Ilaser = \t" << bundle.Ilaser()
				<< "\n#Ixray = \t" << bundle.Ixray()
				<< "\n#center laser = \t" << z_laser.real() << "\t" << z_laser.imag() 
				<< "\n#center xray = \t" << z_xray.real() << "\t" << z_xray.imag()
				<< "\n#alpha = \t" << bundle.delay_angle() 
				<< std::endl;
			interferestream << "#";
			pulsearray[0]->printwavelengthbins(&interferestream);

			for (size_t i=0;i<bundle.get_nfibers();i++){
				*pulsearray[i] -= *crosspulsearray[i];
				*pulsearray[i] *= bundle.Ilaser(size_t(i)); 
				pulsearray[i]->appendwavelength(&interferestream);
			}
			interferestream.close();

			DebugOps::pushout(std::string("At end of images loop in thread" + std::to_string(tid) + "\t"),n);


		} // outermost loop for nimages to produce //

		for (size_t i=0;i<pulsearray.size();++i){
			delete pulsearray[i];
			delete crosspulsearray[i];
		}
		delete masterpulsePtr;
		delete bundlePtr;
	} // end parallel region
	// file for delay bins
	filename = filebase + "delaybins.out";
	ofstream outbins(filename.c_str(),ios::out); 
	for (size_t i=0;i<bundle.get_nfibers();++i){
		outbins << bundle.delay(i) << "\n";
	}
	outbins.close();

	return 0;
}


