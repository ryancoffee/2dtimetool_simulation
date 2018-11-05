#include "scan_material.hpp"
#include <cstdlib>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <complex>
#include <memory>
#include <fftw3.h>

#include <vector>
#include <random>
#include <chrono>

#include <DebugOps.hpp>
#include <ScanParams.hpp>

#include <cstddef> // for nullptr

using namespace Constants;

/* Here's main */
int main(int argc, char* argv[])
{
	std::time_t tstart = std::time(nullptr);
	std::cout << "\t\t======================================\n"
		<< "\t\t======= scan_material started ========\n"
		<< "\t\t===== " << std::asctime(std::localtime(&tstart)) 
		<< "\t\t======================================\n" << std::flush;

	unsigned nthreads = (unsigned)atoi( getenv("nthreads") );
	std::cerr << "Scaling fibers =\t";
	if (getenv("scalefibers")){
		std::cerr << "yes\n";
	} else {
		std::cerr << "no\n";
	}
	std::cerr << std::flush;

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

	std::cerr << "\t\tshuffle fibers?\t";
	if (getenv("shuffle_fibers"))
	{
		std::cerr << "yes\n";masterbundle.shuffle_output();
	} else {
		std::cerr << "no\n";
	}

	masterbundle.Ixray(float(1.));
	masterbundle.Ilaser(float(1.));
	std::string filename = scanparams.filebase() + "fibermap.out";
	std::cout << "fibermap file = " << filename << std::endl << std::flush;
	std::ofstream mapfile(filename.c_str(),std::ios::out);
	//masterbundle.print_mapping(mapfile,double(0.0));
	mapfile.close();

	// file for delay bins
	ofstream outbins(std::string(scanparams.filebase() + "delaybins.out").c_str(),ios::out); 
	for (size_t f=0;f<masterbundle.get_nfibers();++f){
		outbins << masterbundle.delay(f) << "\n";
	}
	outbins.close();

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

	std::cout << "initializing masterpulse and masterplans" << std::endl << std::flush;
	PulseFreq masterpulse(scanparams.omega0(),scanparams.omega_width(),scanparams.omega_onoff(),scanparams.tspan());
	fftw_plan forward;
	fftw_plan backward;
	fftw_plan plan_r2hc;
	fftw_plan plan_hc2r;
	fftw_plan plan_r2hc_2x;
	fftw_plan plan_hc2r_2x;
	masterpulse.setmasterplans(&forward,&backward);
	masterpulse.setancillaryplans(& plan_r2hc,& plan_hc2r,& plan_r2hc_2x,& plan_hc2r_2x);
	masterpulse.addchirp(scanparams.getchirp());							// chirp that ref pulse


	// Setup the shared pulse arrays

	CalibMat calibration(boost::lexical_cast<size_t>(atoi(getenv("ncalibdelays")))
			, boost::lexical_cast<double>(atof(getenv("fsWindow"))));
	if (!getenv("skipcalibration"))
	{
		calibration.set_center(boost::lexical_cast<double>(atof(getenv("delays_mean"))));
		std::cout << "\n\n\t\t ====== delays =======\n";
		for (size_t i = 0 ; i< calibration.get_ndelays(); ++i){
			std::cout << calibration.get_delay(i) << " ";
		}
		std::cout << std::endl << std::flush;
	}

	std::vector< PulseFreq* > calpulsearray(calibration.get_ndelays(),NULL);

#pragma omp parallel num_threads(nthreads) shared(calpulsearray) 
	//#pragma omp parallel num_threads(nthreads) shared(calpulsearray,masterresponse,masterpulse,scanparams) 
	//#pragma omp parallel num_threads(nthreads) 
	{ // begin parallel region
		size_t tid = omp_get_thread_num();

		MatResponse calibresponse(masterresponse);
		std::vector< PulseFreq* > pulsearray(masterbundle.get_nfibers(),NULL);
		std::vector< PulseFreq* > crosspulsearray(masterbundle.get_nfibers(),NULL);
		for (size_t f=0;f<pulsearray.size();++f){
			pulsearray[f] = new PulseFreq(masterpulse);
			crosspulsearray[f] = new PulseFreq(masterpulse);
		}
		if (!getenv("skipcalibration"))
		{
#pragma omp for schedule(static) ordered 
			for (size_t d=0;d<calpulsearray.size();++d)
			{ // outermost loop for calibration.get_ndelays() to produce //
				std::cerr << "\tinside parallel region for actual loop\td = " << d << "\ttid = " << tid << std::endl << std::flush;
				//initialize with masterpulse
				PulseFreq calpulse(masterpulse);
				PulseFreq* calpulsePtr = &calpulse;
				PulseFreq calcrosspulse(masterpulse);
				PulseFreq* calcrosspulsePtr = &calcrosspulse;

				double startdelay(calibration.get_delay(d));

				calcrosspulsePtr->delay(scanparams.interferedelay()); // delay in the frequency domain
				calpulsePtr->fft_totime();
				calcrosspulsePtr->fft_totime();

				for(size_t g=0;g<scanparams.ngroupsteps();g++){ // begin groupsteps loop
					calibresponse.setdelay(startdelay - g*scanparams.groupstep()); // forward propagating, x-rays advance on the optical
					calibresponse.setstepvec_amp(calpulsePtr);
					calibresponse.setstepvec_phase(calpulsePtr);
					calibresponse.setstepvec_amp(calcrosspulsePtr);
					calibresponse.setstepvec_phase(calcrosspulsePtr);
					if (scanparams.doublepulse()){
						calibresponse.addstepvec_amp(calpulsePtr,scanparams.doublepulsedelay());
						calibresponse.addstepvec_phase(calpulsePtr,scanparams.doublepulsedelay());
						calibresponse.addstepvec_amp(calcrosspulsePtr,scanparams.doublepulsedelay());
						calibresponse.addstepvec_phase(calcrosspulsePtr,scanparams.doublepulsedelay());
					}
					// this pulls down the tail of the response so vector is periodic on nsamples	
					calibresponse.buffervectors(calpulsePtr); 
					calibresponse.buffervectors(calcrosspulsePtr); 
					calpulsePtr->modulateamp_time();
					calpulsePtr->modulatephase_time();
					calcrosspulsePtr->modulateamp_time();
					calcrosspulsePtr->modulatephase_time();
				}// end groupsteps loop



				for (size_t e=0;e<scanparams.netalon();e++){ // begin etalon loop
					// back propagation step //
					double etalondelay = startdelay - double(e+1) * (calibresponse.getetalondelay()); 
					// at front surface, x-rays see counter-propagating light from one full etalon delay

					PulseFreq etalonpulse=*(calpulsePtr);
					PulseFreq crossetalonpulse=*(calcrosspulsePtr);

					for(size_t g=0;g<scanparams.ngroupsteps();g++){
						calibresponse.setdelay(etalondelay + g*scanparams.backstep()); 
						// counterpropagating, x-rays work backwards through the optical

						calibresponse.setstepvec_amp(etalonpulse);
						calibresponse.setstepvec_phase(etalonpulse);
						calibresponse.setstepvec_amp(crossetalonpulse);
						calibresponse.setstepvec_phase(crossetalonpulse);
						if (scanparams.doublepulse()){
							calibresponse.addstepvec_amp(etalonpulse,scanparams.doublepulsedelay());
							calibresponse.addstepvec_phase(etalonpulse,scanparams.doublepulsedelay());
							calibresponse.addstepvec_amp(crossetalonpulse,scanparams.doublepulsedelay());
							calibresponse.addstepvec_phase(crossetalonpulse,scanparams.doublepulsedelay());
						}
						calibresponse.buffervectors(etalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						calibresponse.buffervectors(crossetalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						etalonpulse.modulateamp_time();
						etalonpulse.modulatephase_time();
						crossetalonpulse.modulateamp_time();
						crossetalonpulse.modulatephase_time();
					}
					// forward propagation //
					for(size_t g=0;g<scanparams.ngroupsteps();g++){
						calibresponse.setdelay(startdelay - g*scanparams.groupstep()); // forward propagating, x-rays advance on the optical
						calibresponse.setstepvec_amp(etalonpulse);
						calibresponse.setstepvec_phase(etalonpulse);
						calibresponse.setstepvec_amp(crossetalonpulse);
						calibresponse.setstepvec_phase(crossetalonpulse);
						if (scanparams.doublepulse()){
							calibresponse.addstepvec_amp(etalonpulse,scanparams.doublepulsedelay());
							calibresponse.addstepvec_phase(etalonpulse,scanparams.doublepulsedelay());
							calibresponse.addstepvec_amp(crossetalonpulse,scanparams.doublepulsedelay());
							calibresponse.addstepvec_phase(crossetalonpulse,scanparams.doublepulsedelay());

						}
						calibresponse.buffervectors(etalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						calibresponse.buffervectors(crossetalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						etalonpulse.modulateamp_time();
						etalonpulse.modulatephase_time();
						crossetalonpulse.modulateamp_time();
						crossetalonpulse.modulatephase_time();
					}
					etalonpulse.fft_tofreq();
					crossetalonpulse.fft_tofreq();
					etalonpulse.delay(calibresponse.getetalondelay()); // delay and attenuate in frequency domain
					etalonpulse.attenuate(pow(calibresponse.getreflectance(),(int)2));
					crossetalonpulse.delay(calibresponse.getetalondelay()); // delay and attenuate in frequency domain
					crossetalonpulse.attenuate(pow(calibresponse.getreflectance(),(int)2));
					etalonpulse.fft_totime();
					crossetalonpulse.fft_totime();
					*(calpulsePtr) += etalonpulse;
					*(calcrosspulsePtr) += crossetalonpulse;
				} // end etalon loop


				calpulsePtr->fft_tofreq();
				calcrosspulsePtr->fft_tofreq();
				calpulsePtr->delay(scanparams.interferedelay()); // expects this in fs // time this back up to the crosspulse

				calpulse -= calcrosspulse;
				// reversing order for sake of chirp calib matrix
				calpulsearray[calpulsearray.size()-d-1] = new PulseFreq(calpulse);
			} // end of loop calibration.get_ndelays() to produce //

#pragma omp flush 

#pragma omp master
			{
				// print out the calibration as ascii for now //
				// print rows in order, eventually in tf_record or matrix or so. //
				std::string calfilename = scanparams.filebase() + "interference.calibration";
				ofstream calibrationstream(calfilename.c_str(),ios::out); 
				std::cout << "\tcalibration filename out = " << calfilename << std::endl;
				calibrationstream << "#";
				calpulsearray[0]->printwavelengthbins(&calibrationstream);

				for (size_t n=0;n<calpulsearray.size();++n){
					calpulsearray[n]->appendwavelength(&calibrationstream);
					delete calpulsearray[n];
					calpulsearray[n] = NULL;
				}
				calibrationstream.close();
				std::cout << "Finished with the calibration image/matrix" << std::endl << std::flush;

				if (scanparams.addrandomphase(atoi(getenv("addrandomphase"))>0))
				{
					masterpulse.addrandomphase();
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
		}

		if (!getenv("skipimages"))
		{

#pragma omp for schedule(static) ordered 
			for (size_t n=0;n<scanparams.nimages();++n)
			{ // outermost loop for nimages to produce //
				if (n==0) {
					std::cout << "\n\t ==== http://www.fftw.org/fftw3_doc/Advanced-Complex-DFTs.html ===="
						<< "\n\t ==== use this for defining multiple fibers as ===="
						<< "\n\t ==== contiguous blocks for row-wise FFT as 2D ====\n" << std::flush;
				}

				for (size_t f=0;f<pulsearray.size();++f){
					*(pulsearray[f]) = masterpulse;
					*(crosspulsearray[f]) = masterpulse;
				}

				double t0 = scanparams.delays_uniform();
				double startdelay(0);

				FiberBundle parabundle(masterbundle);


				parabundle.Ixray(scanparams.xray_inten_rand());
				parabundle.Ilaser(scanparams.laser_inten_rand());
				parabundle.delay_angle(scanparams.dalpha()*double(n));
				parabundle.center_Ixray(scanparams.xray_pos_rand(),scanparams.xray_pos_rand());
				parabundle.center_Ilaser(scanparams.laser_pos_rand(),scanparams.laser_pos_rand());



				if (tid==0){
					DebugOps::pushout(std::string("Running for t0 = " + std::to_string(t0) + "in threaded for loop, thread " + std::to_string(tid)));
				}
				std::string mapfilename = scanparams.filebase() + "fibermap.out." + std::to_string(n);
				//std::cout << "fibermap file = " << mapfilename << std::endl << std::flush;
				std::ofstream mapfile(mapfilename.c_str(),std::ios::out);
				parabundle.print_mapping(mapfile,t0);
				mapfile.close();


				for(size_t f = 0; f < parabundle.get_nfibers(); f++)
				{ // begin fibers loop
					startdelay = t0 + parabundle.delay(f);
					pulsearray[f]->scale(parabundle.Ilaser(f));
					crosspulsearray[f]->scale(parabundle.Ilaser(f));

					MatResponse pararesponse(masterresponse);

					if (getenv("scalefibers")){
						pararesponse.setscale(parabundle.Ixray(f));
						if (tid==0){
							std::cerr << "parabundle.Ixray(" << f << ") = " << parabundle.Ixray(f) << std::endl << std::flush;
						}
					}

					if (scanparams.addchirpnoise()){
						std::vector<double> noise(scanparams.getchirpnoise());
						pulsearray[f]->addchirp(noise); 
						crosspulsearray[f]->addchirp(noise); 
					}

					crosspulsearray[f]->delay(scanparams.interferedelay()); // delay in the frequency domain
					pulsearray[f]->fft_totime();
					crosspulsearray[f]->fft_totime();

					for(size_t g=0;g<scanparams.ngroupsteps();g++){ // begin groupsteps loop
						pararesponse.setdelay(startdelay - g*scanparams.groupstep()); // forward propagating, x-rays advance on the optical
						pararesponse.setstepvec_amp(pulsearray[f]);
						pararesponse.setstepvec_phase(pulsearray[f]);
						pararesponse.setstepvec_amp(crosspulsearray[f]);
						pararesponse.setstepvec_phase(crosspulsearray[f]);
						if (scanparams.doublepulse()){
							pararesponse.addstepvec_amp(pulsearray[f],scanparams.doublepulsedelay());
							pararesponse.addstepvec_phase(pulsearray[f],scanparams.doublepulsedelay());
							pararesponse.addstepvec_amp(crosspulsearray[f],scanparams.doublepulsedelay());
							pararesponse.addstepvec_phase(crosspulsearray[f],scanparams.doublepulsedelay());
						}
						// this pulls down the tail of the response so vector is periodic on nsamples	
						pararesponse.buffervectors(pulsearray[f]); 
						pararesponse.buffervectors(crosspulsearray[f]); 
						pulsearray[f]->modulateamp_time();
						pulsearray[f]->modulatephase_time();
						crosspulsearray[f]->modulateamp_time();
						crosspulsearray[f]->modulatephase_time();
					}// end groupsteps loop


					for (size_t e=0;e<scanparams.netalon();e++){ // begin etalon loop
						//std::cerr << "\n\t\t ---- starting etalon at " << e << " ----" << std::endl;
						// back propagation step //
						double etalondelay = startdelay - double(e+1) * (pararesponse.getetalondelay()); 
						// at front surface, x-rays see counter-propagating light from one full etalon delay

						PulseFreq etalonpulse=*(pulsearray[f]);
						PulseFreq crossetalonpulse=*(crosspulsearray[f]);

						for(size_t g=0;g<scanparams.ngroupsteps();g++){
							pararesponse.setdelay(etalondelay + g*scanparams.backstep()); 
							// counterpropagating, x-rays work backwards through the optical

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
						*(pulsearray[f]) += etalonpulse;
						*(crosspulsearray[f]) += crossetalonpulse;
					} // end etalon loop


					pulsearray[f]->fft_tofreq();
					crosspulsearray[f]->fft_tofreq();
					pulsearray[f]->delay(scanparams.interferedelay()); // expects this in fs // time this back up to the crosspulse
				} // end nfibers loop


				std::string filename = scanparams.filebase() + "interference.out." + std::to_string(n);// + ".tid" + std::to_string(tid);
				ofstream interferestream(filename.c_str(),ios::out); // use app to append delays to same file.

				//std::cout << "tid = " << tid << ": interfere filename out = " << filename << std::endl;
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
				pulsearray[0]->printwavelengthbins(&interferestream);
				for (size_t f=0;f<parabundle.get_nfibers();f++){
					*(pulsearray[f]) -= *(crosspulsearray[f]);
					pulsearray[f]->scale(parabundle.Ilaser(f)); 
					if (tid == 0){
						int max = boost::lexical_cast<double>(getenv("gain")) * std::pow(pulsearray[f]->maxsignal(),int(2));
						for (size_t i=0;i<std::log(max);++i){
							std::cout << '.';
						}
						std::cout << "\n" << std::flush;
					}
					pulsearray[f]->appendwavelength(&interferestream);
				}
				interferestream.close();

				// std::cerr << "\n\t ---- moving to next image in nimages loop ----" << std::endl;
			} // outermost loop for nimages to produce //
			for (size_t f=0;f<pulsearray.size();++f){
				delete pulsearray[f];
				pulsearray[f] = NULL;
				delete crosspulsearray[f];
				crosspulsearray[f] = NULL;
			}

		}

		//std::cout << "\n\t... trying to leave parallel region" << std::endl;
	} // end parallel region
	//std::cout << "\n ---- just left parallel region ----" << std::endl;
	std::cout << "masterresponse reflectance: " << masterresponse.getreflectance() << std::endl;
	std::cout << "masterbundle fiberdiameter: " << masterbundle.fiberdiameter() << std::endl;
	std::cout << "scanparams lambda_0: " << scanparams.lambda_0() << std::endl;
	fftw_destroy_plan(forward);
	fftw_destroy_plan(backward);

	std::time_t tstop = std::time(nullptr);
	std::cout << "\t\t======================================\n"
		<< "\t\t======== scan_material stopped =======\n"
		<< "\t\t===== " << std::asctime(std::localtime(&tstop)) 
		<< "\t\t===== in " << int(tstop-tstart) << " s \n"
		<< "\t\t======================================\n" << std::flush;

	return 0;
}


