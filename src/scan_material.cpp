#include "scan_material.hpp"
#include <cstdlib>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <complex>

#include <vector>
#include <random>
#include <chrono>

#include <DebugOps.hpp>
#include <ScanParams.hpp>

using namespace Constants;

/* Here's main */
int main(int argc, char* argv[])
{
	unsigned nthreads = (unsigned)atoi( getenv("nthreads") );

	ScanParams scanparams;
	scanparams.nimages(size_t(atoi(getenv("nimages"))));
	scanparams.filebase(std::string(getenv("filebase")));

	scanparams.dalpha((atof(getenv("drifting_alpha")))*pi<double>()/scanparams.nimages());

	if (scanparams.doublepulse(atoi(getenv("doublepulse")) > 0 )){
		scanparams.doublepulsedelay(atof( getenv("doublepulsedelay") ) ) ; // this one gets used directly in atomic units
	}
	scanparams.lambda_0(atof( getenv("lambda0") ));
	scanparams.lambda_width( atof( getenv("lambda_width") ));
	scanparams.lambda_onoff( atof( getenv("lambda_onoff") ));
	scanparams.tspan((atof( getenv("tspan") ) )/fsPau<double>());

	scanparams.ngroupsteps(atoi( getenv("ngroupsteps") ));
	scanparams.groupdelay(atof(getenv("groupdelay")));
	scanparams.backdelay(atof(getenv("backdelay")));
	scanparams.netalon(atoi(getenv("netalon")));


	scanparams.etalonreflectance(atof(getenv("etalon")));
	scanparams.etalondelay(atof(getenv("etalondelay")));
	scanparams.interferedelay((double)atof(getenv("interferedelay")));

	scanparams.chirp(
			( atof( getenv("chirp") ) ) / std::pow(fsPau<float>(),int(2)), // the difference in slopes at omega_low versus omega_high must equal tspan
			( atof( getenv("TOD") ) ) / std::pow(fsPau<float>(),int(3)),
			( atof( getenv("FOD") ) ) / std::pow(fsPau<float>(),int(4)),
			( atof( getenv("fifthOD") ) ) / std::pow(fsPau<float>(),int(5))
			);


	if (scanparams.addchirpnoise(atoi(getenv("usechirpnoise"))>0)){
		scanparams.initchirpnoise( 
				( atof( getenv("chirpnoise") ) ) / std::pow(fsPau<float>(),int(2)), 
				( atof( getenv("TODnoise") ) ) / std::pow(fsPau<float>(),int(3)),
				( atof( getenv("FODnoise") ) ) / std::pow(fsPau<float>(),int(4)),
				( atof( getenv("fifthODnoise") ) ) / std::pow(fsPau<float>(),int(5))
				);
	}



	FiberBundle masterbundle(boost::lexical_cast<size_t>(atoi(getenv("nfibers"))));
	masterbundle.fiberdiameter(boost::lexical_cast<float>(atof(getenv("fiberdiam"))));
	masterbundle.laserdiameter(boost::lexical_cast<float>(atof(getenv("laserdiam"))));
	masterbundle.xraydiameter(boost::lexical_cast<float>(atof(getenv("xraydiam"))));
	masterbundle.set_fsPmm(boost::lexical_cast<float>(atof(getenv("bundle_fsPmm"))));
	masterbundle.scalePolarCoords();
	masterbundle.shuffle_output();
	masterbundle.Ixray(float(1.));
	masterbundle.Ilaser(float(1.));
	std::string filename = scanparams.filebase() + "fibermap.out";
	std::cerr << filename << std::endl;
	std::ofstream mapfile(filename.c_str(),std::ios::out);
	masterbundle.print_mapping(mapfile);
	mapfile.close();




	std::cerr << "just now setting masterresponse()" << std::endl;
	
	// file for delay bins
	ofstream outbins(std::string(scanparams.filebase() + "delaybins.out").c_str(),ios::out); 
	for (size_t i=0;i<masterbundle.get_nfibers();++i){
		outbins << masterbundle.delay(i) << "\n";
	}
	outbins.close();

	std::cerr << "Setting original masterresponse" << std::endl;
	MatResponse masterresponse(
			0,															// stepdelay
			(double)( atof( getenv("stepwidth") ) ),								// stepwidth
			((double)( atof( getenv("attenuation") ) ) - 1.0) * masterbundle.Ixray() / scanparams.ngroupsteps() + 1.0,	// attenuation
			(double)( atof( getenv("phase") ) ) * masterbundle.Ixray() / scanparams.ngroupsteps()				// phase
			);
	masterresponse.aalphabbeta(
			(double)( atof( getenv("a") ) ),		// a
			(double)( atof( getenv("alpha" ) ) ),		// alpha
			(double)( atof( getenv("b") ) ),		// b
			(double)( atof( getenv("beta") ) )		// beta
			);

	masterresponse.setreflectance(scanparams.etalonreflectance());
	masterresponse.setetalondelay(scanparams.etalondelay());

	std::cerr << "Just finished modifying chirp for original masterresponse" << std::endl;

#pragma omp parallel num_threads(nthreads) shared(masterbundle,masterresponse,scanparams)
//#pragma omp parallel num_threads(nthreads) 
	{ // begin parallel region
	size_t tid = omp_get_thread_num();
	if (tid==0) std::cerr << "inside parallel region for tid = " << tid << std::endl;

	std::cerr << "made it here" << std::endl;
	PulseFreq masterpulse(scanparams.omega0(),scanparams.omega_width(),scanparams.omega_onoff(),scanparams.tspan());
	std::cerr << tid << " = " << atoi(getenv("addrandomphase")) << std::endl;
	masterpulse.addchirp(scanparams.getchirp());							// chirp that ref pulse

	if (scanparams.addrandomphase(atoi(getenv("addrandomphase"))>0) ){
		masterpulse.addrandomphase();
		if (tid ==0){
			std::string filename = scanparams.filebase() + "spectralphaseFTpower.dat";
			std::ofstream outfile(filename.c_str(),std::ios::out);
			masterpulse.print_phase_powerspectrum(outfile);
			outfile.close();
			filename = scanparams.filebase() + "spectralphase.dat";
			outfile.open(filename.c_str(),std::ios::out);
			masterpulse.print_phase(outfile);
			outfile.close();
			filename = scanparams.filebase() + "spectralamp.dat";
			outfile.open(filename.c_str(),std::ios::out);
			masterpulse.print_amp(outfile);
			outfile.close();
		}
	}

#pragma omp for //schedule(static) // ordered 
	for (size_t n=0;n<scanparams.nimages();++n){ // outermost loop for nimages to produce //
		if (tid==0) std::cerr << "threading on image " << n << " for tid = " << tid << std::endl;

		if (tid=0) std::cerr << " HERE HERE HERE HERE Come back to here, probably go back to vector and use push_back(), likely fftw plan alloc" << std::endl;

		PulseFreq * pulsearray = (PulseFreq*) std::malloc(masterbundle.get_nfibers() * sizeof(masterpulse));
		PulseFreq * crosspulsearray = (PulseFreq*) std::malloc(masterbundle.get_nfibers() * sizeof(masterpulse));
		//std::vector< PulseFreq * > pulsearray(masterbundle.get_nfibers(),NULL);
		//std::vector< PulseFreq * > crosspulsearray(masterbundle.get_nfibers(),NULL);
		if (pulsearray == NULL) 
			std::cerr << "pulsearray is NULL !" << std::endl;
		if (tid==0) std::cerr << "just used masterbundle.get_nfibers() = " << masterbundle.get_nfibers() << std::endl;
		for (size_t i=0;i<masterbundle.get_nfibers();++i){
			std::cerr << "(f" << i << ",t" << tid << ")" << std::endl;
			pulsearray[i] = masterpulse;
			crosspulsearray[i] = masterpulse;
		}
		std::cerr << "sizeof(masterpulse) " << sizeof(masterpulse) << " sizeof(pulsearray[2]) " << sizeof(pulsearray[2]) << std::endl;

		if (tid==0) std::cerr << "pulsearray & crosspulsearray populated by masterpulse for tid = " << tid << std::endl;


		double t0 = scanparams.delays_uniform();
		double startdelay(0);

		FiberBundle parabundle(masterbundle);
		parabundle.Ixray(scanparams.xray_inten_rand());
		parabundle.Ilaser(scanparams.laser_inten_rand());
		parabundle.delay_angle(scanparams.dalpha()*double(n));
		parabundle.center_Ixray(scanparams.xray_pos_rand(),scanparams.xray_pos_rand());
		parabundle.center_Ilaser(scanparams.laser_pos_rand(),scanparams.laser_pos_rand());

		MatResponse pararesponse(masterresponse);


			if (tid==0){
				std::cout << "Running for t0 = " << t0 << std::endl;
				std::cerr << "\n\n\t\t===== scanparams.delays_uniform() = " << scanparams.delays_uniform() << std::endl;
				DebugOps::pushout(std::string("in threaded for loop, thread" + std::to_string(tid)));
			}


			for(size_t i = 0; i< parabundle.get_nfibers(); i++)
			{ // begin fibers loop
				startdelay = t0 + parabundle.delay(i);

				if (scanparams.addchirpnoise()){
					std::vector<double> noise(scanparams.getchirpnoise());
					pulsearray[i].addchirp(noise); 
					crosspulsearray[i].addchirp(noise); 
				}

				crosspulsearray[i].delay(scanparams.interferedelay()); // delay in the frequency domain
				pulsearray[i].fft_totime();
				crosspulsearray[i].fft_totime();


				for(size_t g=0;g<scanparams.ngroupsteps();g++){ // begin groupsteps loop
					pararesponse.setdelay(startdelay - g*scanparams.groupstep()); // forward propagating, x-rays advance on the optical
					pararesponse.setstepvec_amp(pulsearray[i]);
					pararesponse.setstepvec_phase(pulsearray[i]);
					pararesponse.setstepvec_amp(crosspulsearray[i]);
					pararesponse.setstepvec_phase(crosspulsearray[i]);
					if (scanparams.doublepulse()){
						pararesponse.addstepvec_amp(pulsearray[i],scanparams.doublepulsedelay());
						pararesponse.addstepvec_phase(pulsearray[i],scanparams.doublepulsedelay());
						pararesponse.addstepvec_amp(crosspulsearray[i],scanparams.doublepulsedelay());
						pararesponse.addstepvec_phase(crosspulsearray[i],scanparams.doublepulsedelay());
					}
					// this pulls down the tail of the response so vector is periodic on nsamples	
					pararesponse.buffervectors(pulsearray[i]); 
					pararesponse.buffervectors(crosspulsearray[i]); 
					pulsearray[i].modulateamp_time();
					pulsearray[i].modulatephase_time();
					crosspulsearray[i].modulateamp_time();
					crosspulsearray[i].modulatephase_time();
				}// end groupsteps loop


				for (size_t e=0;e<scanparams.netalon();e++){ // begin etalon loop
					// back propagation step //
					double etalondelay = startdelay - double(e+1) * (pararesponse.getetalondelay()); // at front surface, x-rays see counter-propagating light from one full etalon delay
					PulseFreq etalonpulse(pulsearray[i]);
					PulseFreq crossetalonpulse(crosspulsearray[i]);

					for(size_t g=0;g<scanparams.ngroupsteps();g++){
						pararesponse.setdelay(etalondelay + g*scanparams.backstep()); // counterpropagating, x-rays work backwards through the optical
						pararesponse.setstepvec_amp(etalonpulse);
						pararesponse.setstepvec_phase(etalonpulse);
						pararesponse.setstepvec_amp(crossetalonpulse);
						pararesponse.setstepvec_phase(crossetalonpulse);
						if (scanparams.doublepulse()){
							pararesponse.addstepvec_amp(etalonpulse,scanparams.doublepulsedelay());
							pararesponse.addstepvec_phase(etalonpulse,scanparams.doublepulsedelay());
							pararesponse.addstepvec_amp(crossetalonpulse,scanparams.doublepulsedelay());
							pararesponse.addstepvec_phase(crossetalonpulse,scanparams.doublepulsedelay());
						}
						pararesponse.buffervectors(etalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						pararesponse.buffervectors(crossetalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						etalonpulse.modulateamp_time();
						etalonpulse.modulatephase_time();
						crossetalonpulse.modulateamp_time();
						crossetalonpulse.modulatephase_time();
					}
					// forward propagation //
					for(size_t g=0;g<scanparams.ngroupsteps();g++){
						pararesponse.setdelay(startdelay - g*scanparams.groupstep()); // forward propagating, x-rays advance on the optical
						pararesponse.setstepvec_amp(etalonpulse);
						pararesponse.setstepvec_phase(etalonpulse);
						pararesponse.setstepvec_amp(crossetalonpulse);
						pararesponse.setstepvec_phase(crossetalonpulse);
						if (scanparams.doublepulse()){
							pararesponse.addstepvec_amp(etalonpulse,scanparams.doublepulsedelay());
							pararesponse.addstepvec_phase(etalonpulse,scanparams.doublepulsedelay());
							pararesponse.addstepvec_amp(crossetalonpulse,scanparams.doublepulsedelay());
							pararesponse.addstepvec_phase(crossetalonpulse,scanparams.doublepulsedelay());

						}
						pararesponse.buffervectors(etalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						pararesponse.buffervectors(crossetalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
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
					pulsearray[i] += etalonpulse;
					crosspulsearray[i] += crossetalonpulse;

				} // end etalon loop


				pulsearray[i].fft_tofreq();
				crosspulsearray[i].fft_tofreq();

				pulsearray[i].delay(scanparams.interferedelay()); // expects this in fs // time this back up to the crosspulse

			} // end nfibers loop


		std::string filename = scanparams.filebase() + "interference.out." + std::to_string(n);
		ofstream interferestream(filename.c_str(),ios::out); // use app to append delays to same file.

		std::cout << "interfere filename out = " << filename << std::endl;
		std::complex<double> z_laser = parabundle.center_Ilaser();
		std::complex<double> z_xray = parabundle.center_Ixray();
		interferestream << "#delay for image = \t" << t0 
			<< "\n#Ilaser = \t" << parabundle.Ilaser()
			<< "\n#Ixray = \t" << parabundle.Ixray()
			<< "\n#center laser = \t" << z_laser.real() << "\t" << z_laser.imag() 
			<< "\n#center xray = \t" << z_xray.real() << "\t" << z_xray.imag()
			<< "\n#alpha = \t" << parabundle.delay_angle() 
			<< std::endl;
		interferestream << "#";
		pulsearray[0].printwavelengthbins(&interferestream);
		for (size_t i=0;i<parabundle.get_nfibers();i++){
			pulsearray[i] -= crosspulsearray[i];
			pulsearray[i] *= parabundle.Ilaser(size_t(i)); 
			pulsearray[i].appendwavelength(&interferestream);
		}
		interferestream.close();


/*
	for (size_t i=0;i<masterbundle.get_nfibers();i++){
		delete pulsearray[i];
		delete crosspulsearray[i];
	}
*/
	std::free(pulsearray);
	std::free(crosspulsearray);

	} // outermost loop for nimages to produce //

	} // end parallel region

	return 0;
}


