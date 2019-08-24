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
	std::cerr << "\n\t\tHERE HERE HERE HERE\n\n" << std::flush;
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
	std::cerr << "\n\t\tHERE HERE before FiberBundle HERE HERE\n\n" << std::flush;
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
	masterbundle.thermaldiameter(boost::lexical_cast<float>(atof(getenv("thermaldiam"))));
	masterbundle.center_thermal(boost::lexical_cast<float>(atof(getenv("thermalcenter_x"))), boost::lexical_cast<float>(atof(getenv("thermalcenter_y"))));

	masterbundle.setTmax_Tbase(boost::lexical_cast<float>(atof(getenv("TmaxK"))),boost::lexical_cast<float>(atof(getenv("TbaseK"))));
	masterbundle.set_fsPmm(boost::lexical_cast<float>(atof(getenv("bundle_fsPmm"))));
	masterbundle.scalePolarCoords();

	std::cout << "\t\tshuffle fibers?\t";
	if (getenv("shuffle_fibers"))
	{
		std::cout << "yes\n";masterbundle.shuffle_output();
	} else {
		std::cout << "no\n";
	}

	std::cerr << "HERE HERE HERE HERE\n\t\tWHy can't we shuffle fibers on the home machines? (pavoni)\n" << std::flush;
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
	masterresponse.bandgap(double(atof(getenv("bandgap_eV")))); //
	if (getenv("usediamond")){
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
					/*
					calibresponse.setstepvec_amp(calpulse);
					calibresponse.setstepvec_phase(calpulse);
					calibresponse.setstepvec_amp(calcrosspulse);
					calibresponse.setstepvec_phase(calcrosspulse);
					*/
					if (scanparams.doublepulse()){
						calibresponse.addstepvec_both_carriers(calpulse,scanparams.doublepulsedelay());
						calibresponse.addstepvec_both_carriers(calcrosspulse,scanparams.doublepulsedelay());
						/*
						calibresponse.addstepvec_amp(calpulse,scanparams.doublepulsedelay());
						calibresponse.addstepvec_phase(calpulse,scanparams.doublepulsedelay());
						calibresponse.addstepvec_amp(calcrosspulse,scanparams.doublepulsedelay());
						calibresponse.addstepvec_phase(calcrosspulse,scanparams.doublepulsedelay());
						*/
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
							/*
							calibresponse.addstepvec_amp(etalonpulse,scanparams.doublepulsedelay());
							calibresponse.addstepvec_phase(etalonpulse,scanparams.doublepulsedelay());
							calibresponse.addstepvec_amp(crossetalonpulse,scanparams.doublepulsedelay());
							calibresponse.addstepvec_phase(crossetalonpulse,scanparams.doublepulsedelay());
							*/
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
						/*
						calibresponse.setstepvec_amp(etalonpulse);
						calibresponse.setstepvec_phase(etalonpulse);
						calibresponse.setstepvec_amp(crossetalonpulse);
						calibresponse.setstepvec_phase(crossetalonpulse);
						*/
						if (scanparams.doublepulse()){
							calibresponse.addstepvec_both_carriers(etalonpulse,scanparams.doublepulsedelay());
							calibresponse.addstepvec_both_carriers(crossetalonpulse,scanparams.doublepulsedelay());
							/*
							calibresponse.addstepvec_amp(etalonpulse,scanparams.doublepulsedelay());
							calibresponse.addstepvec_phase(etalonpulse,scanparams.doublepulsedelay());
							calibresponse.addstepvec_amp(crossetalonpulse,scanparams.doublepulsedelay());
							calibresponse.addstepvec_phase(crossetalonpulse,scanparams.doublepulsedelay());
							*/

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

#pragma omp parallel num_threads(nthreads) default(shared) shared(masterpulse)
		{
	if (!getenv("skipimages"))
	{
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
#pragma omp for schedule(dynamic)
			for (size_t n=0;n<scanparams.nimages();++n)
			{ // outermost loop for nimages to produce //
				//std::cerr << "\tinside the parallel region 2 for images loop n = " << n << " in thread " << tid << "\n" << std::flush;
				if (n==0 & tid==0) {
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



				//DebugOps::pushout(std::string("Running image " + std::to_string(n) + " for t0 = " + std::to_string(t0) + " in threaded for loop, thread " + std::to_string(tid)));
				std::string mapfilename = scanparams.filebase() + "fibermap.out." + std::to_string(n);
				//std::cout << "fibermap file = " << mapfilename << std::endl << std::flush;
				std::ofstream mapfile(mapfilename.c_str(),std::ios::out);
				parabundle.print_mapping(mapfile,t0);
				mapfile.close();


				for(size_t f = 0; f < parabundle.get_nfibers(); f++)
				{ // begin fibers loop
					pulse = masterpulse;
					crosspulse = masterpulse;
					startdelay = t0 + parabundle.delay(f);
					pulse.scale(parabundle.Ilaser(f));
					crosspulse.scale(parabundle.Ilaser(f));

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
						pararesponse.setstepvec_amp(pulse);
						pararesponse.setstepvec_phase(pulse);
						pararesponse.setstepvec_amp(crosspulse);
						pararesponse.setstepvec_phase(crosspulse);
						if (scanparams.doublepulse()){
							pararesponse.addstepvec_amp(pulse,scanparams.doublepulsedelay());
							pararesponse.addstepvec_phase(pulse,scanparams.doublepulsedelay());
							pararesponse.addstepvec_amp(crosspulse,scanparams.doublepulsedelay());
							pararesponse.addstepvec_phase(crosspulse,scanparams.doublepulsedelay());
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
						//std::cerr << "tid = " << tid << "\tpulse/crosspulse.domain() = " << pulse.domain() << "/" << crosspulse.domain() << "\n" << std::flush;
						//std::cerr << "\n\t\t ---- starting etalon at " << e << " ----\n" << std::flush;
						// back propagation step //
						double etalondelay = startdelay - double(e+1) * pararesponse.thermaletalondelay(parabundle.TinK(f)); 
						// at front surface, x-rays see counter-propagating light from one full etalon delay

						etalonpulse = pulse;
						crossetalonpulse = crosspulse;
						//std::cerr << "etalonpulse/crossetalonpulse.domain() = " << etalonpulse.domain() << "/" << crossetalonpulse.domain() << "\n" << std::flush;

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
						//std::cerr << "\t\t\t ########### // forward propagation // #############\n" << std::flush;
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
						//std::cerr << "etalonpulse/crossetalonpulse.domain() = " << etalonpulse.domain() << "/" << crossetalonpulse.domain() << "\n" << std::flush;
						etalonpulse.fft_tofreq();
						crossetalonpulse.fft_tofreq();
						//etalonpulse.delay(pararesponse.getetalondelay() * ( 1. + parabundle.ThermalEtalonDelta(f)) ); // delay and attenuate in frequency domain
						etalonpulse.delay(pararesponse.thermaletalondelay(parabundle.TinK(f))) ; // delay and attenuate in frequency domain
						etalonpulse.attenuate(pow(pararesponse.getreflectance(),(int)2));
						crossetalonpulse.delay(pararesponse.thermaletalondelay(parabundle.TinK(f)));
								//pararesponse.getetalondelay() * (1. + parabundle.ThermalEtalonDelta(f)) ); // delay and attenuate in frequency domain
						crossetalonpulse.attenuate(pow(pararesponse.getreflectance(),(int)2));
						etalonpulse.fft_totime();
						crossetalonpulse.fft_totime();
						pulse += etalonpulse;
						crosspulse += crossetalonpulse;
					} // end etalon loop


					pulse.fft_tofreq();
					crosspulse.fft_tofreq();
					pulse.delay(scanparams.interferedelay()); // expects this in fs // time this back up to the crosspulse
					pulse -= crosspulse;
					// std::cerr << "\n\n\t\t\t\t============== testing... just before the push_back() ==============\n\n" << std::flush;
					pulsearray[f] = pulse;
				} // end nfibers loop


				std::string filename = scanparams.filebase() + "interference.out." + std::to_string(n);// + ".tid" + std::to_string(tid);
				//std::cerr << "testing: filename = " << filename << "\n" << std::flush;
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
				pulsearray[0].printwavelengthbins(&interferestream);
				for (size_t f=0;f<pulsearray.size();f++){
					pulsearray[f].scale(parabundle.Ilaser(f)); 
					pulsearray[f].appendwavelength(&interferestream);
				}
/*
	HERE HERE HERE HERE
	Here we need to figure out how to run a rolling TDI style image superposition so that each thread 
	produces a time-ordered series of images who are comprised of the most recent n_overlay images that
	are superimposed with a 1x row-shift
	Looks like maybe we need to write another Pulse.accumwavelength(&interferestream,nshift)
*/
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
				interferestream.close();

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


