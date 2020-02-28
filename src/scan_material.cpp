#include "scan_material.hpp"
#include <cstdlib>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <complex>
#include <memory>
#include <fftw3.h>

#include <cstdint> // for sake of defining uint16_t for the OpenCV mat to be filled.
// OpenCV includes
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/imgcodecs.hpp>

#include <vector>
#include <random>
#include <chrono>

#include <DebugOps.hpp>
#include <ScanParams.hpp>

#include <cstddef> // for nullptr
#include <cstdint> // for writing as an int32_t and int16_t

using namespace Constants;

/* Here's main */
int main(int argc, char* argv[])
{
	std::time_t tstart = std::time(nullptr);
	std::cout << "\t\t======================================\n"
		<< "\t\t======= scan_material started ========\n"
		<< "\t\t===== " << std::asctime(std::localtime(&tstart)) 
		<< "\t\t===== on host " << getenv("HOSTNAME") << "\n"
		<< "\t\t===== for " << getenv("nimages") << " images\n"
		<< "\t\t===== with " << getenv("nfibers") << " fibers\n"
		<< "\t\t======================================\n" << std::flush;

	unsigned nthreads = (unsigned)atoi( getenv("nthreads") );
	std::cout << "Scaling fibers =\t";
	if (getenv("scale_fibers")){
		std::cout << "yes\n";
	} else {
		std::cout << "no\n";
	}
	std::cout << std::flush;


	ScanParams scanparams;
	scanparams.nimages(size_t(atoi(getenv("nimages"))));
	scanparams.filebase(std::string(getenv("filebase")));
	scanparams.calfilebase(std::string(getenv("calfilebase")));

	std::vector<float> imagetimes(scanparams.nimages(),0); // for benchmarking the processors

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
	scanparams.etalondelay(double(atof(getenv("etalondelay"))));
	scanparams.interferedelay((double)atof(getenv("interferedelay")));
	scanparams.interferephase((double)atof(getenv("interferephase")));

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
	masterbundle.thermaldiameter(boost::lexical_cast<float>(atof(getenv("thermaldiam"))));
	masterbundle.center_thermal(boost::lexical_cast<float>(atof(getenv("thermalcenter_x"))), boost::lexical_cast<float>(atof(getenv("thermalcenter_y"))));

	masterbundle.setTmax_Tbase(boost::lexical_cast<float>(atof(getenv("TmaxK"))),boost::lexical_cast<float>(atof(getenv("TbaseK"))));
	masterbundle.set_fsPmm(boost::lexical_cast<float>(atof(getenv("bundle_fsPmm"))));
	masterbundle.scalePolarCoords();

	std::vector<uint16_t>keyinds(masterbundle.get_nfibers());
	DataOps::ramp(keyinds);
	std::cout << "\t\tshuffle fibers?\t";
	if (getenv("shuffle_fibers"))
	{
		std::cout << "yes\n";
		std::random_device rng;
		std::seed_seq seed{rng(), rng(), rng(), rng(), rng(), rng(), rng(), rng()};
		std::mt19937 e(seed);
		std::shuffle(keyinds.begin(),keyinds.end(),e);
		/*
		if (!masterbundle.shuffle_inds())
			std::cerr << "The rerturn from masterbundle.shuffle_output() was false\t-- something failed\n" << std::flush;
		*/
	} else {
		std::cout << "no\n";
	}
	ofstream outkey(std::string(scanparams.filebase() + "fiberkey.out").c_str(),ios::out);
	outkey << std::flush;
	outkey.close();
	masterbundle.set_inds(keyinds);

	masterbundle.Ixray(float(1.));
	masterbundle.Ilaser(float(1.));
	std::string filename = scanparams.filebase() + "fibermap.out";
	std::cout << "fibermap file = " << filename << std::endl << std::flush;
	std::ofstream mapfile(filename.c_str(),std::ios::out);
	masterbundle.print_mapping(mapfile,double(0.0));
	mapfile.close();

	// file for delay bins
	if (!getenv("skipcalibration")){
		ofstream outbins(std::string(scanparams.filebase() + "delaybins.out").c_str(),ios::out); 
		for (size_t f=0;f<masterbundle.get_nfibers();++f){
			outbins << masterbundle.delay(f) << "\n";
			outbins.close();
		}
	}

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
	masterresponse.set_thickness(double(atof(getenv("l_thickness"))));
	masterresponse.set_n_refractive(double(atof(getenv("n_0"))),
				double(atof(getenv("n_a"))),
				double(atof(getenv("n_b")))
				);


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

	double xrayphoton_energy = double(atof(getenv("xrayphoton_energy")));
	std::cout << " xrayphoton_energy = " << xrayphoton_energy << " keV\n" << std::flush;
	masterresponse.bandgap(double(atof(getenv("bandgap_eV")))); //
	if (getenv("usediamond")){
		std::cerr << "using xrayphoton_energy to compute carriers\n" << std::flush;
		masterresponse.fill_carriersvec(masterpulse,xrayphoton_energy);
	} else {
		std::string carriersfilename = getenv("carriersfile");
		std::cerr << "carriersfilename = " << carriersfilename << "\n" << std::flush;
		std::ifstream Nikita_file(carriersfilename.c_str(),std::ios::in);
		masterresponse.fill_carriersvec(masterpulse,Nikita_file);
	}


	std::time_t tstop = std::time(nullptr);
	std::cout << "\tIt has taken " << (tstop-tstart) << " s so far for initializing masterpulse and building fftw plans\n" << std::flush;

	CalibMat calibration(boost::lexical_cast<size_t>(atoi(getenv("ncalibdelays")))
			, boost::lexical_cast<double>(atof(getenv("fsWindow"))));
	if (!getenv("skipcalibration"))
	{
		std::cout << "\t\t############ entering calibration ###########\n" << std::flush;
		calibration.set_center(boost::lexical_cast<double>(atof(getenv("delays_mean"))));
		std::cout << "\t\t====== delays =======\n";
		for (size_t i = 0 ; i< calibration.get_ndelays(); ++i){
			std::cout << calibration.get_delay(i) << " ";
		}
		std::cout << std::endl << std::flush;
	}

	if (!getenv("skipcalibration"))
	{
		// Setup the shared pulse arrays
		std::vector< PulseFreq > calpulsearray(calibration.get_ndelays(),masterpulse);

#pragma omp parallel num_threads(nthreads) default(shared) shared(masterpulse)
		{ // begin parallel region 1
			size_t tid = omp_get_thread_num();

			// all non-shared objects must be created inside the parallel section for default is shared if defined outside
			// http://pages.tacc.utexas.edu/~eijkhout/pcse/html/omp-data.html
			PulseFreq etalonpulse(masterpulse);
			PulseFreq crossetalonpulse(masterpulse);
			PulseFreq calpulse(masterpulse);
			PulseFreq calcrosspulse(masterpulse);

			// initialize with masterpulse/masterresponse
			MatResponse calibresponse(masterresponse);



#pragma omp for schedule(dynamic) 
			for (size_t d=0;d<calpulsearray.size();++d)
			{ // outermost loop for calibration.get_ndelays() to produce //
				//std::cout << "\tinside parallel region for actual loop d = " 
				//<< d << "\twith tid = " << tid << "\n" << std::flush;
				std::cout << '+' << std::flush;
				//before each delay, reset to masterpulse
				//std::cerr << "\tmasterpulse.domain() = " << masterpulse.domain() << "\n" << std::flush;
				calpulse = masterpulse;
				calcrosspulse = masterpulse;
				//std::cerr << "\tcalpulse.domain() = " << calpulse.domain() << "\n" << std::flush;
				//std::cerr << "\tcalcrosspulse.domain() = " << calcrosspulse.domain() << "\n" << std::flush;
				//std::cerr << "\n\n\t\t###### Made it here HERE HERE HERE debugging seg fault ######\n\n" << std::flush;

				double startdelay(calibration.get_delay(d));

				calcrosspulse.delay(scanparams.interferedelay()); // delay in the frequency domain

				calpulse.fft_totime();
				calcrosspulse.fft_totime();

				//std::cerr << "\tafter fft_totime() calpulse.domain() = " << calpulse.domain() << "\n" << std::flush;
				//std::cerr << "\tafter fft_totime() calcrosspulse.domain() = " << calcrosspulse.domain() << "\n" << std::flush;

				for(size_t g=0;g<scanparams.ngroupsteps();g++){ // begin groupsteps loop
					calibresponse.setdelay(startdelay - g*scanparams.groupstep()); // forward propagating, x-rays advance on the optical
					calibresponse.setstepvec_both_carriers(calpulse);
					calibresponse.setstepvec_both_carriers(calcrosspulse);
					if (scanparams.doublepulse()){
						calibresponse.addstepvec_both_carriers(calpulse,scanparams.doublepulsedelay());
						calibresponse.addstepvec_both_carriers(calcrosspulse,scanparams.doublepulsedelay());
					}
					// this pulls down the tail of the response so vector is periodic on nsamples	
					calibresponse.buffervectors(calpulse); 
					calibresponse.buffervectors(calcrosspulse); 
					calpulse.modulateamp_time();
					calpulse.modulatephase_time();
					calcrosspulse.modulateamp_time();
					calcrosspulse.modulatephase_time();
				}// end groupsteps loop



				for (size_t e=0;e<scanparams.netalon();e++){ // begin etalon loop
					// back propagation step //
					double etalondelay = startdelay - double(e+1) * (calibresponse.getetalondelay()*(1. + float(d)/float(calpulsearray.size()) ) ); 
					// at front surface, x-rays see counter-propagating light from one full etalon delay

					// reset back to calpulse for each round
					etalonpulse = calpulse;
					crossetalonpulse = calcrosspulse;
					//std::cerr << "\t\t\t -- inside etalon: etalonpulse.domain() = " << etalonpulse.domain() << "\n" << std::flush;
					//std::cerr << "\t\t\t -- inside etalon: crossetalonpulse.domain() = " << crossetalonpulse.domain() << "\n" << std::flush;

					for(size_t g=0;g<scanparams.ngroupsteps();g++){
						calibresponse.setdelay(etalondelay + g*scanparams.backstep()); 
						// counterpropagating, x-rays work backwards through the optical

						calibresponse.setstepvec_both_carriers(etalonpulse);
						calibresponse.setstepvec_both_carriers(crossetalonpulse);
						/*
						calibresponse.setstepvec_amp(etalonpulse);
						calibresponse.setstepvec_phase(etalonpulse);
						calibresponse.setstepvec_amp(crossetalonpulse);
						calibresponse.setstepvec_phase(crossetalonpulse);
						*/
						if (scanparams.doublepulse()){
							calibresponse.addstepvec_both_carriers(etalonpulse,scanparams.doublepulsedelay());
							calibresponse.addstepvec_both_carriers(crossetalonpulse,scanparams.doublepulsedelay());
						}
						calibresponse.buffervectors(etalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						calibresponse.buffervectors(crossetalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						//std::cerr << "\t\t\t -- inside etalon groupsteps: before: etalonpulse.domain() = " << etalonpulse.domain() << "\n" << std::flush;
						//std::cerr << "\t\t\t -- inside etalon groupsteps: before: crossetalonpulse.domain() = " << crossetalonpulse.domain() << "\n" << std::flush;
						etalonpulse.modulateamp_time();
						etalonpulse.modulatephase_time();
						crossetalonpulse.modulateamp_time();
						crossetalonpulse.modulatephase_time();
						//std::cerr << "\t\t\t -- inside etalon groupsteps: after: etalonpulse.domain() = " << etalonpulse.domain() << "\n" << std::flush;
						//std::cerr << "\t\t\t -- inside etalon groupsteps: after: crossetalonpulse.domain() = " << crossetalonpulse.domain() << "\n" << std::flush;
					}
					// forward propagation //
					for(size_t g=0;g<scanparams.ngroupsteps();g++){
						calibresponse.setdelay(startdelay - g*scanparams.groupstep()); // forward propagating, x-rays advance on the optical
						calibresponse.setstepvec_both_carriers(etalonpulse);
						calibresponse.setstepvec_both_carriers(crossetalonpulse);
						if (scanparams.doublepulse()){
							calibresponse.addstepvec_both_carriers(etalonpulse,scanparams.doublepulsedelay());
							calibresponse.addstepvec_both_carriers(crossetalonpulse,scanparams.doublepulsedelay());
						}
						calibresponse.buffervectors(etalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						calibresponse.buffervectors(crossetalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
						//std::cerr << "\t\t\t -- inside etalon groupsteps forward prop: before: etalonpulse.domain() = " << etalonpulse.domain() << "\n" << std::flush;
						//std::cerr << "\t\t\t -- inside etalon groupsteps forward prop: before: crossetalonpulse.domain() = " << crossetalonpulse.domain() << "\n" << std::flush;
						etalonpulse.modulateamp_time();
						etalonpulse.modulatephase_time();
						crossetalonpulse.modulateamp_time();
						crossetalonpulse.modulatephase_time();
						//std::cerr << "\t\t\t -- inside etalon groupsteps forward prop: after: etalonpulse.domain() = " << etalonpulse.domain() << "\n" << std::flush;
						//std::cerr << "\t\t\t -- inside etalon groupsteps forward prop: after: crossetalonpulse.domain() = " << crossetalonpulse.domain() << "\n" << std::flush;
					}
					//std::cerr << "\t\t\t -- inside etalon finished groupsteps: etalonpulse.domain() = " << etalonpulse.domain() << "\ttid = "<< tid << "\n" << std::flush;
					//std::cerr << "\t\t\t -- inside etalon finished groupsteps: crossetalonpulse.domain() = " << crossetalonpulse.domain() << "\n" << std::flush;
					etalonpulse.fft_tofreq();
					crossetalonpulse.fft_tofreq();
					etalonpulse.delay(calibresponse.getetalondelay() * (1. + float(d)/float(calpulsearray.size())) ); // delay and attenuate in frequency domain
					etalonpulse.attenuate(pow(calibresponse.getreflectance(),(int)2));
					crossetalonpulse.delay(calibresponse.getetalondelay() * (1. + float(d)/float(calpulsearray.size())) ); // delay and attenuate in frequency domain
					crossetalonpulse.attenuate(pow(calibresponse.getreflectance(),(int)2));
					etalonpulse.fft_totime();
					crossetalonpulse.fft_totime();
					//std::cerr << "\t\t\t -- end etalon: etalonpulse.domain() = " << etalonpulse.domain() << "\n" << std::flush;
					//std::cerr << "\t\t\t -- end etalon: crossetalonpulse.domain() = " << crossetalonpulse.domain() << "\n" << std::flush;
					//std::cerr << "\t\t\t -- end etalon: calpulse.domain() = " << calpulse.domain() << "\n" << std::flush;
					//std::cerr << "\t\t\t -- end etalon: calcrosspulse.domain() = " << calcrosspulse.domain() << "\n" << std::flush;
					calpulse += etalonpulse;
					calcrosspulse += crossetalonpulse;
				} // end etalon loop


				//std::cerr << "\t\tbefore fft_tofreq(): calpulse.domain() = " << calpulse.domain() << "\n" << std::flush;
				//std::cerr << "\t\tbefore fft_tofreq(): calcrosspulse.domain() = " << calcrosspulse.domain() << "\n" << std::flush;
				calpulse.fft_tofreq();
				calcrosspulse.fft_tofreq();
				calpulse.delay(scanparams.interferedelay()); // expects this in fs // time this back up to the crosspulse

				calpulse -= calcrosspulse;
				// reversing order for sake of chirp calib matrix
				calpulsearray[calpulsearray.size()-d-1] = calpulse;

				/*
				   std::cerr << "\t\tcalpulse.domain() = " << calpulse.domain() << "\n" << std::flush;
				   std::cerr << "\t\tcalcrosspulse.domain() = " << calcrosspulse.domain() << "\n" << std::flush;
				   std::cerr << "\t\tcalpulsearray[ " << (calpulsearray.size()-d-1) << " ].domain() = " 
				   << calpulsearray[calpulsearray.size()-d-1].domain() << "\n" << std::flush;
				   */
			} // end of loop calibration.get_ndelays() to produce //


#pragma omp barrier

#pragma omp master
			{
				std::cout << "|\t done with calibration delays\n" << std::flush;
				//std::cerr << "\t\t###### Made it here too #####\n\t\t##### should call only once in master #######\n" << std::flush;
				// print out the calibration as ascii for now //
				// print rows in order, eventually in tf_record or matrix or so. //
				std::string calfilename = scanparams.calfilebase() + "interference.calibration";
				std::string derivfilename = scanparams.calfilebase() + "interference.calibration.derivative";
				std::string calfilename_delays = scanparams.calfilebase() + "interference.calibration.delays";
				std::string calfilename_wavelengths = scanparams.calfilebase() + "interference.calibration.wavelengths";
				ofstream calibrationstream(calfilename.c_str(),ios::out); 
				ofstream derivstream(derivfilename.c_str(),ios::out); 
				ofstream calibrationstream_delays(calfilename_delays.c_str(),ios::out); 
				ofstream calibrationstream_wavelengths(calfilename_wavelengths.c_str(),ios::out); 
				/*
				std::string bin_calfilename = scanparams.filebase() + "interference.calibration.bin";
				ofstream bin_calibrationstream(bin_calfilename.c_str(),ios::out | ios::binary); 
				std::cout << "\tcalibration filename out = " << calfilename << "\n\t and \t" << bin_calfilename << std::endl;
				*/
				calibrationstream << "# wavelengths\n#";
				derivstream << "# wavelengths\n#";
				calpulsearray[0].printwavelengthbins(&calibrationstream);
				calpulsearray[0].printwavelengthbins(&derivstream);
				calpulsearray[0].printwavelengthbins(&calibrationstream_wavelengths);
				calibrationstream << "# delays\n#";
				derivstream << "# delays\n#";
				calibrationstream_delays << "# delays\n";
				for (size_t i = 0 ; i< calibration.get_ndelays(); ++i){
					calibrationstream << calibration.get_delay(i) << "\t";
					derivstream << calibration.get_delay(i) << "\t";
					calibrationstream_delays << calibration.get_delay(i) << "\t";
				}
				calibrationstream << "\n";
				derivstream << "\n";
				calibrationstream_delays << "\n";

				for (size_t n=0;n<calpulsearray.size();++n){
					calpulsearray[n].appendwavelength(&calibrationstream);
					calpulsearray[n].appendwavelength_deriv(&derivstream);
					//	calpulsearray[n].appendwavelength_bin(&bin_calibrationstream);
				}

				calibrationstream.close();
				derivstream.close();
				calibrationstream_delays.close();
				calibrationstream_wavelengths.close();
				//bin_calibrationstream.close();
				std::cout << "Finished with the calibration image/matrix\n" << std::flush;


			}

#pragma omp master
			{
				std::cout << "\t\t############ ending parallel region 1 ###########\n" << std::flush;
			}
		} // end parallel region 1

	} // end if (!getenv("skipcalibration"))


	//############## Images section ##############

#pragma omp parallel num_threads(nthreads) default(shared) shared(masterpulse,masterbundle,scanparams)
		{
	if (!getenv("skipimages"))
	{
		std::random_device rd{};
		std::mt19937 rng{rd()};
		std::normal_distribution<> xrayshadow_x{double(atof(getenv("xrayshadowcorner_x"))),double(atof(getenv("xrayshadowcorner_xjitter")))};
		std::normal_distribution<> xrayshadow_y{double(atof(getenv("xrayshadowcorner_y"))),double(atof(getenv("xrayshadowcorner_yjitter")))};

		std::cout << "\t\t############ entering parallel/images ###########\n" << std::flush;
			size_t tid = omp_get_thread_num();
			size_t nfibers = masterbundle.get_nfibers();


			FiberBundle parabundle(masterbundle);
			MatResponse pararesponse(masterresponse);

			PulseFreq pulse(masterpulse);
			PulseFreq crosspulse(masterpulse);
			PulseFreq etalonpulse(masterpulse);
			PulseFreq crossetalonpulse(masterpulse);
			std::vector< PulseFreq > pulsearray(nfibers,PulseFreq(masterpulse));
#pragma omp barrier

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


/*
	HERE HERE HERE HERE
	Here we need to figure out how to run a rolling TDI style image superposition so that each thread 
	produces a time-ordered series of images who are comprised of the most recent n_overlay images that
	are superimposed with a 1x row-shift
*/
			//std::cerr << "Entering parallel for loop with nimages = "  << scanparams.nimages() << "\n" << std::flush;
#pragma omp for schedule(dynamic)
			for (size_t n=0;n<scanparams.nimages();++n)
			{ // outermost loop for nimages to produce //
				//std::cerr << "\tinside the parallel region 2 for images loop n = " << n << " in thread " << tid << "\n" << std::flush;
				if (n<nthreads & tid==0) {
					std::cout << "========================================================================="
						<<   "\n\t\t ==== http://www.fftw.org/fftw3_doc/Advanced-Complex-DFTs.html ===="
						<<   "\n\t\t ====         use this for defining multiple fibers as         ===="
						<<   "\n\t\t ====         contiguous blocks for row-wise FFT as 2D         ===="
						<<   "\n\t\t ==================================================================\n" << std::flush;
				}

				std::time_t imgstart = std::time(nullptr);

				double t0 = scanparams.delays_uniform();
				double startdelay(0);

				parabundle = masterbundle;


				parabundle.Ixray(scanparams.xray_inten_rand());
				parabundle.Ilaser(scanparams.laser_inten_rand());
				parabundle.delay_angle(scanparams.dalpha()*double(n));
				parabundle.center_Ixray(scanparams.xray_pos_rand(),scanparams.xray_pos_rand());
				parabundle.center_Ilaser(scanparams.laser_pos_rand(),scanparams.laser_pos_rand());

				parabundle.fillIxray();
				parabundle.shadow_xrays(xrayshadow_x(rng) , xrayshadow_y(rng));


				//DebugOps::pushout(std::string("Running image " + std::to_string(n) + " for t0 = " + std::to_string(t0) + " in threaded for loop, thread " + std::to_string(tid)));
				std::string mapfilename = scanparams.filebase() + "fibermap.out." + std::to_string(n);
				//std::cout << "fibermap file = " << mapfilename << std::endl << std::flush;
				std::ofstream mapfile(mapfilename.c_str(),std::ios::out);
				if (!parabundle.print_mapping(mapfile,t0))
					std::cerr << "Something failed in printing this fibermapping out\n" << std::flush;
				mapfile.close();


				for(size_t f = 0; f < parabundle.get_nfibers(); f++)
				{ // begin fibers loop
					pulse = masterpulse;
					crosspulse = masterpulse;
					startdelay = t0 + parabundle.delay(f);
					pulse.scale(parabundle.Ixray(f));
					crosspulse.scale(parabundle.Ixray(f));
					pararesponse = masterresponse;

					if (getenv("scale_fibers")){
						pararesponse.setscale(parabundle.Ixray(f));
						//std::cerr << "parabundle.Ixray(" << f << ") = " << parabundle.Ixray(f) << "\n" << std::flush;
					}

					if (scanparams.addchirpnoise()){
						std::vector<double> noise(scanparams.getchirpnoise());
						pulse.addchirp(noise); 
						crosspulse.addchirp(noise); 
					}

					crosspulse.delay(scanparams.interferedelay()); // delay in the frequency domain
					pulse.fft_totime();
					crosspulse.fft_totime();

					for(size_t g=0;g<scanparams.ngroupsteps();g++){ // begin groupsteps loop
						pararesponse.setdelay(startdelay - g*scanparams.groupstep()); // forward propagating, x-rays advance on the optical
						pararesponse.setstepvec_both_carriers(pulse,0.,parabundle.Ixray(f));
						pararesponse.setstepvec_both_carriers(crosspulse,0.,parabundle.Ixray(f));
						if (scanparams.doublepulse()){
							pararesponse.addstepvec_both_carriers(pulse,scanparams.doublepulsedelay(),parabundle.Ixray(f));
							pararesponse.addstepvec_both_carriers(crosspulse,scanparams.doublepulsedelay(),parabundle.Ixray(f));
						}
						// this pulls down the tail of the response so vector is periodic on nsamples	
						pararesponse.buffervectors(pulse); 
						pararesponse.buffervectors(crosspulse); 
						pulse.modulateamp_time();
						pulse.modulatephase_time();
						crosspulse.modulateamp_time();
						crosspulse.modulatephase_time();
					}// end groupsteps loop
					//std::cerr << "tid = " << tid << "\tpulse/crosspulse.domain() = " << pulse.domain() << "/" << crosspulse.domain() << "\n" << std::flush;


					for (size_t e=0;e<scanparams.netalon();e++){ // begin etalon loop
						//std::cerr << "\n\t\t ---- starting etalon at " << e << " ----\n" << std::flush;
						//std::cerr << "parabundle.TinK( " << f << " ) is " << parabundle.TinK(f) << "\n" << std::flush;
						// back propagation step //
						double etalondelay = startdelay - double(e+1) * pararesponse.thermaletalondelay(parabundle.TinK(f)); 
						// at front surface, x-rays see counter-propagating light from one full etalon delay

						etalonpulse = pulse;
						crossetalonpulse = crosspulse;
						//std::cerr << "etalonpulse/crossetalonpulse.domain() = " << etalonpulse.domain() << "/" << crossetalonpulse.domain() << "\n" << std::flush;

						for(size_t g=0;g<scanparams.ngroupsteps();g++){
							pararesponse.setdelay(etalondelay + g*scanparams.backstep()); 
							// counterpropagating, x-rays work backwards through the optical
							pararesponse.setstepvec_both_carriers(etalonpulse,0.,parabundle.Ixray(f));
							pararesponse.setstepvec_both_carriers(crossetalonpulse,0.,parabundle.Ixray(f));
							if (scanparams.doublepulse()){
								pararesponse.addstepvec_both_carriers(etalonpulse,scanparams.doublepulsedelay(),parabundle.Ixray(f));
								pararesponse.addstepvec_both_carriers(crossetalonpulse,scanparams.doublepulsedelay(),parabundle.Ixray(f));
							}
							pararesponse.buffervectors(etalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
							pararesponse.buffervectors(crossetalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
							etalonpulse.modulateamp_time();
							etalonpulse.modulatephase_time();
							crossetalonpulse.modulateamp_time();
							crossetalonpulse.modulatephase_time();
						}
						// forward propagation //
						//std::cerr << "\t\t\t ########### // forward propagation // #############\n" << std::flush;
						for(size_t g=0;g<scanparams.ngroupsteps();g++){
							pararesponse.setdelay(startdelay - g*scanparams.groupstep()); // forward propagating, x-rays advance on the optical
							pararesponse.setstepvec_both_carriers(etalonpulse,0.,parabundle.Ixray(f));
							pararesponse.setstepvec_both_carriers(crossetalonpulse,0.,parabundle.Ixray(f));
							if (scanparams.doublepulse()){
								pararesponse.addstepvec_both_carriers(etalonpulse,scanparams.doublepulsedelay(),parabundle.Ixray(f));
								pararesponse.addstepvec_both_carriers(crossetalonpulse,scanparams.doublepulsedelay(),parabundle.Ixray(f));
							}
							pararesponse.buffervectors(etalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
							pararesponse.buffervectors(crossetalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
							etalonpulse.modulateamp_time();
							etalonpulse.modulatephase_time();
							crossetalonpulse.modulateamp_time();
							crossetalonpulse.modulatephase_time();
						}
						//std::cerr << "etalonpulse/crossetalonpulse.domain() = " << etalonpulse.domain() << "/" << crossetalonpulse.domain() << "\n" << std::flush;
						etalonpulse.fft_tofreq();
						crossetalonpulse.fft_tofreq();
						etalonpulse.delay(pararesponse.thermaletalondelay(parabundle.TinK(f))) ; // delay and attenuate in frequency domain
						etalonpulse.attenuate(pow(pararesponse.getreflectance(),(int)2));
						crossetalonpulse.delay(pararesponse.thermaletalondelay(parabundle.TinK(f)));
						crossetalonpulse.attenuate(pow(pararesponse.getreflectance(),(int)2));
						etalonpulse.fft_totime();
						crossetalonpulse.fft_totime();
						pulse += etalonpulse;
						crosspulse += crossetalonpulse;
					} // end etalon loop


					pulse.fft_tofreq();
					crosspulse.fft_tofreq();
					pulse.delay(scanparams.interferedelay()); // expects this in fs // time this back up to the crosspulse
					//crosspulse.scale(0.9);
					pulse -= crosspulse;
					//pulse.interfere(crosspulse,scanparams.interferephase());
					//std::cerr << "\n\n\t\t\t\t============== testing... just before the push_back() ==============\n\n" << std::flush;
					pulsearray[f] = pulse;
				} // end nfibers loop
				//

				size_t img_nsamples(1024);
				size_t img_stride(10);
				uint16_t * imdata = (uint16_t*)std::calloc(pulsearray.size() * img_stride * img_nsamples , sizeof(uint16_t));

				for (size_t f=0;f<pulsearray.size();++f){
					pulsearray[f].scale(parabundle.Ilaser(f)); 
					pulsearray[f].fillrow_uint16(imdata + parabundle.get_key(f) * img_stride * img_nsamples,img_nsamples);
				}
				std::complex<double> z_laser = parabundle.center_Ilaser();
				std::complex<double> z_xray = parabundle.center_Ixray();

				// direct image

				ofstream interferestream;
				std::string filename;

				int kr(2*7 + 1);
				int kc(2*10 + 1);

				cv::Mat imageMat_in(pulsearray.size()*img_stride, img_nsamples, CV_16UC1, imdata );	// imageMat_in is 16bit unsigned data
				cv::Mat imageMat(pulsearray.size()*img_stride, img_nsamples, CV_32FC1);
				imageMat_in.convertTo(imageMat,CV_32FC1);	//imageMat is 32 bit float data

				cv::Mat kernel_raw(cv::Mat::zeros(kr,kc,CV_32FC1));

				// vertical bluring //
				std::vector<float> kblur(kr);
				DataOps::sinsqr(kblur);
				cv::Mat cblur(kr,1,CV_32F,kblur.data());
				std::vector<float> zeros(kc,0.);
				zeros[zeros.size()/2] = 1.;
				cv::flip(cblur*cv::Mat(1,kc,CV_32F,zeros.data()) , kernel_raw , -1);

				// constructing raw output, but blurred in the vertical //
				cv::Mat imageMat_raw(cv::Mat(pulsearray.size()*img_stride, img_nsamples, CV_32FC1));
				cv::Mat rawMatout(cv::Mat(pulsearray.size()*img_stride, img_nsamples, CV_16UC1));
				cv::filter2D(imageMat,imageMat_raw,-1,kernel_raw);

				std::vector<int> compression_params;
				compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
				compression_params.push_back(0);

				if (bool(getenv("printASCIIimages"))){

					filename = scanparams.filebase() + "interference.out." + std::to_string(n);
					interferestream.open(filename.c_str(),ios::out); // use app to append delays to same file.
					interferestream << "#delay for image = \t" << t0 
						<< "\n#Ilaser = \t" << parabundle.Ilaser()
						<< "\n#Ixray = \t" << parabundle.Ixray()
						<< "\n#center laser = \t" << z_laser.real() << "\t" << z_laser.imag() 
						<< "\n#center xray = \t" << z_xray.real() << "\t" << z_xray.imag()
						<< "\n#alpha = \t" << parabundle.delay_angle() 
						<< std::endl;
					interferestream << "#";
					pulsearray[0].printwavelengthbins(&interferestream);
					for (size_t f=0;f<pulsearray.size();f++){
						pulsearray[f].scale(parabundle.Ilaser(f));
						pulsearray[f].appendwavelength(&interferestream);
					}
					interferestream.close();
				} else {
					filename = scanparams.filebase() + "interference.params." + std::to_string(n);
					interferestream.open(filename.c_str(),ios::out); // use app to append delays to same file.
					interferestream << "#delay for image = \t" << t0 
						<< "\n#Ilaser = \t" << parabundle.Ilaser()
						<< "\n#Ixray = \t" << parabundle.Ixray()
						<< "\n#center laser = \t" << z_laser.real() << "\t" << z_laser.imag() 
						<< "\n#center xray = \t" << z_xray.real() << "\t" << z_xray.imag()
						<< "\n#alpha = \t" << parabundle.delay_angle() 
						<< "\n#img_stride = \t" << img_stride 
						<< std::endl;
					interferestream << "#";
					pulsearray[0].printwavelengthbins(&interferestream);
					parabundle.print_mapping(interferestream,t0);
					interferestream.close();
					double min,max,scale,offset;
					cv::minMaxLoc(imageMat_raw,&min,&max);
					offset = -min; //*scale;
					imageMat_raw.convertTo(rawMatout,CV_16UC1,1.0,offset);
					cv::minMaxLoc(imageMat_raw,&min,&max);
					scale = float(std::pow(int(2),int(16))-2)/(max-min); // HERE HERE HERE HERE trying to get saturation to be protected
					imageMat_raw.convertTo(rawMatout,CV_16UC1,scale,0.0);
					cv::flip(rawMatout,rawMatout,0);
					filename = scanparams.filebase() + "interference.image." + std::to_string(n) + ".png";
					cv::imwrite(filename.c_str(),rawMatout,compression_params);
				}
				
				// HERE HERE HERE HERE // 

				// initialize kernels vector //
				const unsigned nkernels = 6; // hard coding for now since the storage will be in 2x3channel bgra png + the k0 as greyscale image
				std::vector<cv::Mat> kernels;
				for (unsigned k = 0; k< nkernels; ++k){
					kernels.push_back(cv::Mat::zeros(kr,kc,CV_32FC1));
				}

				std::vector< float > leg(kc,0.);
				for (unsigned k = 0 ; k<nkernels; ++k){		// filling kernels
					DataOps::legendre( leg, k);
					cv::flip(cblur*cv::Mat(1,kc,CV_32F,leg.data()) , kernels[k] , 1);
				}

				if (tid==0 and n<nthreads){ // print the kernels, but only once
					/*
					 * OK, we should use Grahm-Schmidt to come up with orthogonal set, start with sin defined from 0..pi, then cos, then sin*cos, then sin**2, 
					 * then cos**2 but maybe just neg of sin**2, then sin**2*cos... and so forth
					 */
					std::string kfilename;
					ofstream kernelstream;
					for (unsigned k = 0 ; k < nkernels; ++k){
						kfilename = scanparams.filebase() + "kernel" + std::to_string(k);
						kernelstream.open(kfilename.c_str(),ios::out); 
						for (size_t r=0;r<kernels[k].rows;++r){
							for (size_t c=0;c<kernels[k].cols;++c){
								kernelstream << kernels[k].at<float>(r,c) << "\t";
							}
							kernelstream << "\n";
						}
						kernelstream.close(); 
					}
				} // end printing kernels



				std::vector< cv::Mat > imageMat_vec;
				for (unsigned k = 0; k < nkernels; ++k){	// setting up imageMat_vec
					imageMat_vec.push_back(cv::Mat(pulsearray.size()*img_stride, img_nsamples, CV_32FC1));
				}

				const unsigned nchannels = 3; // Stop using alpha channel, that is just awkward
				std::vector< cv::Mat > imageMatout_batch;
				for (unsigned b = 0 ; b < nkernels/nchannels; ++b){	// setting up imageMatout_batch
					imageMatout_batch.push_back(cv::Mat(imageMat_vec[0].rows/8, imageMat_vec[0].cols/8, CV_16UC3));
				}

				for (unsigned i=0;i<nkernels;++i){ // filling imageMat_vec
					cv::filter2D(imageMat, imageMat_vec[i], -1, kernels[i]);
				}
				std::vector<cv::Mat> imageMatout_vec;

				for (unsigned k = 0 ; k < nkernels ; ++k ){
					double min,max,scale,offset;
					cv::minMaxLoc(imageMat_vec[k],&min,&max);
					scale = float(std::pow(int(2),int(16))-1)/(max-min);
					offset = -min*scale;
					imageMatout_vec.push_back(cv::Mat(imageMat_vec[k].rows,imageMat_vec[k].cols,CV_16UC1));
					imageMat_vec[k].convertTo(imageMatout_vec[k],CV_16UC1,scale,offset);
				}

				if ( !(getenv("skipdisplayframes")) and tid==0 ) {
					char FrameStr[15];
					sprintf(FrameStr,"Frame_%i",int(tid));
					cv::namedWindow(FrameStr,cv::WINDOW_NORMAL);
					cv::resizeWindow(FrameStr,img_nsamples*5,pulsearray.size()*5);
					cv::imshow(FrameStr, imageMatout_vec[0]);
					cv::waitKey(0);
					cv::destroyAllWindows();
				}

				std::string pngfilename;
				for (unsigned b = 0; b< imageMatout_batch.size(); ++b){
					std::vector<cv::Mat> imageMat_4vec(nchannels);
					for (unsigned i=0; i<imageMatout_batch[b].channels(); ++i){
						imageMat_4vec[i] = imageMatout_vec[i + b * imageMatout_batch[b].channels()];
						for (unsigned j=0; j<3; ++j) // three times cut both dimensions in half
							cv::pyrDown(imageMat_4vec[i],imageMat_4vec[i],cv::Size(imageMat_4vec[i].cols/2,imageMat_4vec[i].rows/2));
					}
					cv::merge(imageMat_4vec,imageMatout_batch[b]);
					cv::flip(imageMatout_batch[b],imageMatout_batch[b],0);
					pngfilename = scanparams.filebase() + "interference.out.batch" + std::to_string(b) + "." + std::to_string(n) + ".png";
					cv::imwrite(pngfilename.c_str(),imageMatout_batch[b],compression_params);
				}


				/*
				for (unsigned i=0; i<imageMat_4vec.size(); ++i){
					imageMat_4vec[i] = imageMatout_vec[i+1+4];
				}
				cv::merge(imageMat_4vec,imageMat_4chan);
				cv::flip(imageMat_4chan,imageMat_4chan,0);
				pngfilename = scanparams.filebase() + "interference.out.k5678." + std::to_string(n) + ".png";
				cv::imwrite(pngfilename.c_str(),imageMat_4chan,compression_params);
				


				// kernel0
				filename = scanparams.filebase() + "interference.out.K0." + std::to_string(n);
				interferestream.open(filename.c_str(),ios::out); // use app to append delays to same file.
				interferestream << "#delay for image = \t" << t0 
					<< "\n#Ilaser = \t" << parabundle.Ilaser()
					<< "\n#Ixray = \t" << parabundle.Ixray()
					<< "\n#center laser = \t" << z_laser.real() << "\t" << z_laser.imag() 
					<< "\n#center xray = \t" << z_xray.real() << "\t" << z_xray.imag()
					<< "\n#alpha = \t" << parabundle.delay_angle() 
					<< std::endl;
				interferestream << "#";
				pulsearray[0].printwavelengthbins(&interferestream);
				for (size_t r=0;r<imageMatK0.rows;++r){
					for (size_t c=0; c<imageMatK0.cols;++c)
						interferestream << imageMatK0.at<float>(r,c) << "\t";	
					interferestream << "\n";
				}
				interferestream.close();
				// kernel1
				filename = scanparams.filebase() + "interference.out.K1." + std::to_string(n);
				interferestream.open(filename.c_str(),ios::out); 
				interferestream << "#delay for image = \t" << t0 
					<< "\n#Ilaser = \t" << parabundle.Ilaser()
					<< "\n#Ixray = \t" << parabundle.Ixray()
					<< "\n#center laser = \t" << z_laser.real() << "\t" << z_laser.imag() 
					<< "\n#center xray = \t" << z_xray.real() << "\t" << z_xray.imag()
					<< "\n#alpha = \t" << parabundle.delay_angle() 
					<< std::endl;
				interferestream << "#";
				pulsearray[0].printwavelengthbins(&interferestream);
				for (size_t r=0;r<imageMatK1.rows;++r){
					for (size_t c=0; c<imageMatK1.cols;++c)
						interferestream << imageMatK1.at<float>(r,c) << "\t";	
					interferestream << "\n";
				}
				interferestream.close();
				// kernel2
				filename = scanparams.filebase() + "interference.out.K2." + std::to_string(n);
				interferestream.open(filename.c_str(),ios::out); // use app to append delays to same file.
				interferestream << "#delay for image = \t" << t0 
					<< "\n#Ilaser = \t" << parabundle.Ilaser()
					<< "\n#Ixray = \t" << parabundle.Ixray()
					<< "\n#center laser = \t" << z_laser.real() << "\t" << z_laser.imag() 
					<< "\n#center xray = \t" << z_xray.real() << "\t" << z_xray.imag()
					<< "\n#alpha = \t" << parabundle.delay_angle() 
					<< std::endl;
				interferestream << "#";
				pulsearray[0].printwavelengthbins(&interferestream);
				for (size_t r=0;r<imageMatK2.rows;++r){
					for (size_t c=0; c<imageMatK2.cols;++c)
						interferestream << imageMatK2.at<float>(r,c) << "\t";	
					interferestream << "\n";
				}
				interferestream.close();
				// kernel3
				filename = scanparams.filebase() + "interference.out.K3." + std::to_string(n);
				interferestream.open(filename.c_str(),ios::out); // use app to append delays to same file.
				interferestream << "#delay for image = \t" << t0 
					<< "\n#Ilaser = \t" << parabundle.Ilaser()
					<< "\n#Ixray = \t" << parabundle.Ixray()
					<< "\n#center laser = \t" << z_laser.real() << "\t" << z_laser.imag() 
					<< "\n#center xray = \t" << z_xray.real() << "\t" << z_xray.imag()
					<< "\n#alpha = \t" << parabundle.delay_angle() 
					<< std::endl;
				interferestream << "#";
				pulsearray[0].printwavelengthbins(&interferestream);
				for (size_t r=0;r<imageMatK2.rows;++r){
					for (size_t c=0; c<imageMatK2.cols;++c)
						interferestream << imageMatK3.at<float>(r,c) << "\t";	
					interferestream << "\n";
				}
				interferestream.close();
				*/

				std::free(imdata); // this may be able to free right after making hte cv::Mat for this.

				if (tid % 10 < 2){
					for (size_t f=0;f<parabundle.get_nfibers();f++){
						int max = boost::lexical_cast<double>(getenv("gain")) * pulsearray[f].maxsignal();
						for (size_t i=0;i<std::log(max);++i){
							std::cout << '.';
						}
						std::cout << "|";
					}
					std::cout << "\timg = " << n << " in tid = " << tid << "\n" << std::flush;
				}

				std::time_t imgstop = std::time(nullptr);
				imagetimes[n] = float(imgstop - imgstart);

			} // outermost loop for nimages to produce //

#pragma omp master
			{
				std::cout << "\t\t############ ending parallel region 2 ###########\n" << std::flush;
			}

#pragma omp barrier

			//std::cerr << "\n\t... trying to leave parallel region 2" << std::endl;
	} // end if (!getenv("skipimages")
		} // end parallel region

	//std::cout << "\n ---- just left parallel region ----" << std::endl;
	std::cout << "masterresponse reflectance: " << masterresponse.getreflectance() << std::endl;
	std::cout << "masterbundle fiberdiameter: " << masterbundle.fiberdiameter() << std::endl;
	std::cout << "scanparams lambda_0: " << scanparams.lambda_0() << std::endl;
	fftw_destroy_plan(forward);
	fftw_destroy_plan(backward);

	tstop = std::time(nullptr);
	tstop -= tstart;
	std::cout << "\t\t======================================\n"
		<< "\t\t======== scan_material stopped =======\n"
		<< "\t\t===== " << std::asctime(std::localtime(&tstop)) 
		<< "\t\t===== in " << tstop << " s \n"
		<< "\t\t======================================\n" << std::flush;
	std::string timesfilename = scanparams.filebase() + "runtimes.log";
	std::ofstream timesout(timesfilename.c_str(),std::ios::app);
	timesout << "#########################################\n" << std::flush;
	timesout << "# HOSTNAME:\t" << getenv("HOSTNAME") << "\n" << std::flush;
	timesout << "# total time (seconds):\t" << tstop << "\n" << std::flush;
	timesout << "# nfibers :\t" << masterbundle.get_nfibers() << "\n" << std::flush;
	timesout << "# nthreads :\t" << nthreads << "\n" << std::flush;
	timesout << "# mean time / image:\t" << DataOps::mean(imagetimes) 
		<< "\t/ fiber\t" << DataOps::mean(imagetimes)/float(masterbundle.get_nfibers()) << "\n" << std::flush;
	for (size_t i=0 ;i< imagetimes.size();++i){
		timesout << imagetimes[i] << "\t";
	}
	timesout << "\n" << std::flush;
	timesout.close();

	return 0;
}


