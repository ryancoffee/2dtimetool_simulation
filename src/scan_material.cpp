#include "scan_material.hpp"
#include <DebugOps.hpp>
#include <ScanParams.hpp>

using namespace Constants;


/***************/
/* Here's main */
/***************/

int main(int argc, char* argv[])
{
	std::time_t tstart = std::time(nullptr);
	std::cout << "\t\t======================================\n"
		<< "\t\t======= scan_material started ========\n"
		<< "\t\t===== " << std::asctime(std::localtime(&tstart)) << "====\n"
		<< "\t\t===== on host " << getenv("HOSTNAME") << "====\n"
		<< "\t\t===== for " << getenv("nimages") << " images ====\n"
		<< "\t\t===== with " << getenv("nfibers") << " fibers ====\n"
		<< "\t\t======================================\n" << std::endl << std::flush;

	unsigned nthreads = (unsigned)atoi( getenv("nthreads") );
	std::cout << "Scaling fibers =\t";
	if (getenv("scale_fibers")){
		std::cout << "yes\n";
	} else {
		std::cout << "no\n";
	}
	std::cout << std::endl << std::flush;


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

	FiberBundle masterbundle((size_t)(atoi(getenv("nfibers"))));
	masterbundle.fiberdiameter((float)(atof(getenv("fiberdiam"))));
	masterbundle.laserdiameter((float)(atof(getenv("laserdiam"))));
	masterbundle.xraydiameter((float)(atof(getenv("xraydiam"))));
	masterbundle.thermaldiameter((float)(atof(getenv("thermaldiam"))));
	masterbundle.center_thermal((float)(atof(getenv("thermalcenter_x"))), (float)(atof(getenv("thermalcenter_y"))));

	masterbundle.setTmax_Tbase((float)(atof(getenv("TmaxK"))),(float)(atof(getenv("TbaseK"))));
	masterbundle.set_fsPmm((float)(atof(getenv("bundle_fsPmm"))));
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
			(double)( -1.0*atof( getenv("attenuation") ) ) * masterbundle.Ixray() / scanparams.ngroupsteps() ,	// attenuation
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
		std::cout << "using xrayphoton_energy to compute carriers\n" << std::flush;
		masterresponse.fill_carriersvec(masterpulse,xrayphoton_energy);
	} else {
		std::string carriersfilename = getenv("carriersfile");
		std::cout << "carriersfilename = " << carriersfilename << "\n" << std::flush;
		std::ifstream Nikita_file(carriersfilename.c_str(),std::ios::in);
		masterresponse.fill_carriersvec(masterpulse,Nikita_file);
	}


	std::time_t tstop = std::time(nullptr);
	std::cout << "\tIt has taken " << (tstop-tstart) << " s so far for initializing masterpulse and building fftw plans\n" << std::flush;

	CalibMat calibration((size_t)(atoi(getenv("ncalibdelays")))
			, (double)(atof(getenv("fsWindow"))));
	if (!getenv("skipcalibration"))
	{
		std::cout << "\t\t############ entering calibration ###########\n" << std::flush;
		calibration.set_center((double)(atof(getenv("delays_mean"))));
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
				// derivstream << "# wavelengths\n#";
				calpulsearray[0].printwavelengthbins(&calibrationstream);
				calpulsearray[0].printwavelengthbins(&derivstream);
				calpulsearray[0].printwavelengthbins(&calibrationstream_wavelengths);
				calibrationstream << "# delays\n#";
				// derivstream << "# delays\n#";
				calibrationstream_delays << "# delays\n";
				for (size_t i = 0 ; i< calibration.get_ndelays(); ++i){
					calibrationstream << calibration.get_delay(i) << "\t";
					// derivstream << calibration.get_delay(i) << "\t";
					calibrationstream_delays << calibration.get_delay(i) << "\t";
				}
				calibrationstream << "\n";
				// derivstream << "\n";
				calibrationstream_delays << "\n";

				for (size_t n=0;n<calpulsearray.size();++n){
					calpulsearray[n].appendwavelength(&calibrationstream);
					// calpulsearray[n].appendwavelength_deriv(&derivstream);
					// calpulsearray[n].appendwavelength_bin(&bin_calibrationstream);
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


	//############################################//
	//############## Images section ##############//
	//############################################//

/* Commenting out HDF5 for sake of flight, need to build hdf5 build environment for laptop
    */


	auto localnow = std::chrono::system_clock::now();
	std::time_t ttime = std::chrono::system_clock::to_time_t(localnow);
	std::cout << "chrono localtime = " << std::ctime(& ttime) << std::endl;
	std::tm * local_time = std::localtime(& ttime);


	size_t flatimgsize = masterbundle.get_nfibers()*masterpulse.get_freqsamples();
	std::vector< std::vector< uint16_t > > datablock;
	std::vector< double > delays(scanparams.nimages()*nthreads);
	std::vector< double > alphas(scanparams.nimages()*nthreads);
	std::vector< double > ilaser(scanparams.nimages()*nthreads);
	std::vector< double > ixray(scanparams.nimages()*nthreads);
	std::vector< std::array <double,2> > laserposition(scanparams.nimages()*nthreads);
	std::vector< std::array <double,2> > xrayposition(scanparams.nimages()*nthreads);

	for (size_t i=0;i<scanparams.nimages()*nthreads;i++){
		datablock.push_back(std::vector<uint16_t>(masterbundle.get_nfibers()*masterpulse.get_freqsamples()));
	}

#pragma omp parallel num_threads(nthreads) default(shared) shared(masterpulse,masterbundle,scanparams,datablock,ilaser,ixray,laserposition,xrayposition,delays)
	{
        size_t tid = omp_get_thread_num();
		if (!getenv("skipimages"))
		{
			std::random_device rd{};
			std::mt19937 rng{rd()};
			std::normal_distribution<> xrayshadow_x{double(atof(getenv("xrayshadowcorner_x"))),double(atof(getenv("xrayshadowcorner_xjitter")))};
			std::normal_distribution<> xrayshadow_y{double(atof(getenv("xrayshadowcorner_y"))),double(atof(getenv("xrayshadowcorner_yjitter")))};

			size_t nfibers = masterbundle.get_nfibers();
			std::cout << "\t\t############ entering parallel image: tid = " << int(tid) << " ###########\n" << std::flush;


			FiberBundle parabundle(masterbundle);
			MatResponse pararesponse(masterresponse);
			MatResponse pararesponse1(masterresponse);

			PulseFreq pulse(masterpulse);
			PulseFreq crosspulse(masterpulse);
			PulseFreq etalonpulse(masterpulse);
			PulseFreq crossetalonpulse(masterpulse);
			std::vector< PulseFreq > pulsearray(nfibers,PulseFreq(masterpulse));



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



			std::cout << "Entering parallel for loop with nimages = "  << scanparams.nimages() << "\n" << std::flush;
			std::cout << "\t\t#################Entering parallel region 2 #######################\n" << std::flush;



			for (size_t n=0;n<scanparams.nimages();++n)
			{ // outermost loop for nimages to produce //
			  	//std::cerr << "\tinside the parallel region 2 for images loop n = " << n << " in thread " << tid << "\n" << std::flush;

#pragma omp master
				{
					std::cout << "========================================================================="
						<<   "\n\t\t ==== http://www.fftw.org/fftw3_doc/Advanced-Complex-DFTs.html ===="
						<<   "\n\t\t ====         use this for defining multiple fibers as         ===="
						<<   "\n\t\t ====         contiguous blocks for row-wise FFT as 2D         ===="
						<<   "\n\t\t ==================================================================\n" << std::flush;
				}

				std::time_t imgstart = std::time(nullptr);

				double t0 = scanparams.delays_uniform();
				double t1 = t0 + scanparams.rel_delays_normal();
				double startdelay(0);
				double startdelay1(0);

				parabundle = masterbundle;


				parabundle.Ixray(scanparams.xray_inten_rand());
				parabundle.Ilaser(scanparams.laser_inten_rand());
				parabundle.delay_angle(scanparams.dalpha()*double(n));
				parabundle.center_Ixray(scanparams.xray_pos_rand(),scanparams.xray_pos_rand());
				parabundle.center_Ilaser(scanparams.laser_pos_rand(),scanparams.laser_pos_rand());

				parabundle.fillIxray();
				parabundle.shadow_xrays(xrayshadow_x(rng) , xrayshadow_y(rng));


				DebugOps::pushout(std::string("Running image " + std::to_string(n) + " for t0 = " + std::to_string(t0) + " in threaded for loop, thread " + std::to_string(tid)));
				std::string mapfilename = scanparams.filebase() + "fibermap.out." + std::to_string(n);
				//std::cout << "fibermap file = " << mapfilename << std::endl << std::flush;
				std::ofstream mapfile(mapfilename.c_str(),std::ios::out);
				if (!parabundle.print_mapping(mapfile,t0))
					std::cerr << "Something failed in printing this fibermapping out\n" << std::flush;
				mapfile.close();


                std::vector<double> x,y;
                std::vector< uint16_t > imdata(parabundle.get_nfibers()*masterpulse.get_freqsamples());
                /* HERE HERE HERE HERE 
                Fis that you fill from f* nlamsamples insie the fibers loop.
                */
				for(size_t f = 0; f < parabundle.get_nfibers(); f++)
				{ // begin fibers loop
					pulse = masterpulse;
					crosspulse = masterpulse;
					startdelay = t0 + parabundle.delay(f);
					startdelay1 = t1 + parabundle.delay(f);
					pulse.scale(parabundle.Ixray(f));
					crosspulse.scale(parabundle.Ixray(f));
					pararesponse = masterresponse;
					pararesponse1 = masterresponse;

					if (getenv("scale_fibers")){
						pararesponse.setscale(parabundle.Ixray(f));
						pararesponse1.setscale(parabundle.Ixray(f));
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
						pararesponse1.setdelay(startdelay1 - g*scanparams.groupstep()); // forward propagating, x-rays advance on the optical
						pararesponse.setstepvec_both_carriers(pulse,0.,parabundle.Ixray(f));
						pararesponse1.setstepvec_both_carriers(crosspulse,0.,parabundle.Ixray(f));
						if (scanparams.doublepulse()){
							pararesponse.addstepvec_both_carriers(pulse,scanparams.doublepulsedelay(),parabundle.Ixray(f));
							pararesponse1.addstepvec_both_carriers(crosspulse,scanparams.doublepulsedelay(),parabundle.Ixray(f));
						}
						// this pulls down the tail of the response so vector is periodic on nsamples	
						pararesponse.buffervectors(pulse); 
						pararesponse1.buffervectors(crosspulse); 
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
						double etalondelay = startdelay; //- double(e+1) * pararesponse.thermaletalondelay(parabundle.TinK(f)); 
						double etalondelay1 = startdelay1; //- double(e+1) * pararesponse1.thermaletalondelay(parabundle.TinK(f)); 
						// at front surface, x-rays see counter-propagating light from one full etalon delay

						etalonpulse = pulse;
						crossetalonpulse = crosspulse;
						//std::cerr << "etalonpulse/crossetalonpulse.domain() = " << etalonpulse.domain() << "/" << crossetalonpulse.domain() << "\n" << std::flush;

						for(size_t g=0;g<scanparams.ngroupsteps();g++){
							pararesponse.setdelay(etalondelay + g*scanparams.backstep()); 
							pararesponse1.setdelay(etalondelay1 + g*scanparams.backstep()); 
							// counterpropagating, x-rays work backwards through the optical
							pararesponse.setstepvec_both_carriers(etalonpulse,0.,parabundle.Ixray(f));
							pararesponse1.setstepvec_both_carriers(crossetalonpulse,0.,parabundle.Ixray(f));
							if (scanparams.doublepulse()){
								pararesponse.addstepvec_both_carriers(etalonpulse,scanparams.doublepulsedelay(),parabundle.Ixray(f));
								pararesponse1.addstepvec_both_carriers(crossetalonpulse,scanparams.doublepulsedelay(),parabundle.Ixray(f));
							}
							pararesponse.buffervectors(etalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
							pararesponse1.buffervectors(crossetalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
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
								pararesponse1.addstepvec_both_carriers(crossetalonpulse,scanparams.doublepulsedelay(),parabundle.Ixray(f));
							}
							pararesponse.buffervectors(etalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
							pararesponse1.buffervectors(crossetalonpulse); // this pulls down the tail of the response so vector is periodic on nsamples
							etalonpulse.modulateamp_time();
							etalonpulse.modulatephase_time();
							crossetalonpulse.modulateamp_time();
							crossetalonpulse.modulatephase_time();
						}
						//std::cerr << "etalonpulse/crossetalonpulse.domain() = " << etalonpulse.domain() << "/" << crossetalonpulse.domain() << "\n" << std::flush;
						etalonpulse.fft_tofreq();
						crossetalonpulse.fft_tofreq();
						//etalonpulse.delay(pararesponse.thermaletalondelay(parabundle.TinK(f))) ; // delay and attenuate in frequency domain
						etalonpulse.attenuate(pow(pararesponse.getreflectance(),(int)2));
						//crossetalonpulse.delay(pararesponse1.thermaletalondelay(parabundle.TinK(f)));
						crossetalonpulse.attenuate(pow(pararesponse1.getreflectance(),(int)2));
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
					pulsearray[f] = pulse.scale(parabundle.Ilaser(f));
					//std::cerr << double(*(pulse.minmaxvals()).second) << " ";
					//pulsearray[f].scale(parabundle.Ilaser(f)); 
                    if (f==0)
                        x = pulsearray[f].getLamVec(x);
                    pulsearray[f].fillfrequency_bytes(x,
                                                        pulsearray[f].getSigVec(y),
                                                        datablock[scanparams.nimages()*tid + n],
                                                        f);
		    /*
                    pulsearray[f].fillwavelength_bytes(x,
                                                        pulsearray[f].getSigVec(y),
                                                        datablock[scanparams.nimages()*tid + n],
                                                        f);
							*/
				} // end nfibers loop

                /* ########### HERE HERE HERE HERE ############
                    pulsearray is an array of fibers now.
                    Give the flattened image scaled to one byte values and stored flattened with dimensions in a dims vector added to the parameter vector.
                    Do the same for the parameters used, like z_laser and such.
                    */


				std::complex<double> z_laser = parabundle.center_Ilaser();
				std::complex<double> z_xray = parabundle.center_Ixray();

                laserposition[scanparams.nimages()*tid + n]=std::array<double,2>({z_laser.real(),z_laser.imag()});
                xrayposition[scanparams.nimages()*tid + n]=std::array<double,2>({z_xray.real(),z_xray.imag()});
                delay[scanparams.nimages()*tid + n] = double(t0);
                alpha[scanparams.nimages()*tid + n] = double(parabundle.delay_angle();
                ilaser[scanparams.nimages()*tid + n] = double(parabundle.Ilaser();
                ixray[scanparams.nimages()*tid + n] = double(parabundle.Ixray();



				ofstream interferestream;
				if (bool(getenv("printASCIIimages"))){
					std::cout << "\n\n \t\t ###########\tPrinting ASCII images for frame # " << int(n) << "\t#############\n\n" << std::flush;

					filename = scanparams.filebase() + "interference.out." + std::to_string(n);
					interferestream.open(filename.c_str(),ios::out); // use app to append delays to same file.
					interferestream << "#delay for image = \t" << t0 << "and\t" << t1
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

				} 

				std::time_t imgstop = std::time(nullptr);
				imagetimes[n] = float(imgstop - imgstart);

			} // outermost loop for nimages to produce //



            /*##################################################*/
            /*#################### HDF5 storage ################*/
            /*##################################################*/


			//std::cerr << "\n\t... trying to leave parallel region 2" << std::endl;

		} // end if (!getenv("skipimages")

		std::cout << "\t\t############ ending parallel region 2 for tid " << tid << "###########\n" << std::flush;

	} // ends parallel region 2

	std::cout << "\n ---- just left parallel region 2 ----" << std::endl;

	/*
	 *
	 * Attributes of the individual images to be used as truth
	std::vector< double > delays(scanparams.nimages()*nthreads);
	std::vector< double > alphas(scanparams.nimages()*nthreads);
	std::vector< double > ilaser(scanparams.nimages()*nthreads);
	std::vector< double > ixray(scanparams.nimages()*nthreads);
	std::vector< std::array <double,2> > laserposition(scanparams.nimages()*nthreads);
	std::vector< std::array <double,2> > xrayposition(scanparams.nimages()*nthreads);
	 *
	 * The params should be the things we've been so far encoding in the filename
	 */
	if (bool(getenv("H5OUTPUT"))){
		H5::H5File * h5filePtr;
		H5::Group * paramgrp;
		H5::Group * dsetgrp;

		std::string imfilename(scanparams.filebase());
		std::stringstream filetail;
		filetail << "interference.h5";
		imfilename += filetail.str();
		std::cout << "Opening h5 file for intermittent saving:\t" << imfilename << std::endl << std::flush;
		h5filePtr = new H5::H5File ( imfilename , H5F_ACC_TRUNC );
    		std::string name = "/params";
    		paramgrp = new H5::Group( h5filePtr->createGroup( name ) );
    		std::string dsetname = "/datasets";
    		dsetgrp = new H5::Group( h5filePtr->createGroup( dsetname ) );

		const int rank(1);
		size_t dims[1] = {masterbundle.get_nfibers()*masterpulse.get_freqsamples()};
		size_t posdims[1] = {2};
		//size_t dims[1] = {masterbundle.get_nfibers()*masterpulse.get_lamsamples()};
		H5::DataSpace * dataspace = new H5::DataSpace( rank , dims ); 	
		//std::cerr << "Made it here in H5 output\n" << std::flush;
		for (size_t im=0;im<datablock.size();im++){
			H5::DataSet * datasetPtr;
			std::string imname = "/datasets/im_" + std::to_string((int)im);
			datasetPtr = new H5::DataSet( dsetgrp->createDataSet( imname, H5::PredType::NATIVE_USHORT, *dataspace ) );
			datasetPtr->write( datablock[im].data(), H5::PredType::NATIVE_USHORT);
            H5::Attribute * attrPtr;
            ilaserPtr = new H5::Attribute( datasetPtr->createAtribute( "Ilaser", H5::PredType::NATIVE_DOUBLE, ilaser[im]) );
            xrayPtr = new H5::Attribute( datasetPtr->createAtribute( "Ixray", H5::PredType::NATIVE_DOUBLE, ixray[im]) );
            poslaserPtr = new H5::Attribute( datasetPtr->createAtribute( "positionlaser", H5::PredType::NATIVE_DOUBLE, laserposition[im]) );
            posxrayPtr = new H5::Attribute( datasetPtr->createAtribute( "positionxray", H5::PredType::NATIVE_DOUBLE, xrayposition[im]) );
            delayPtr = new H5::Attribute( datasetPtr->createAtribute( "delay", H5::PredType::NATIVE_DOUBLE, delays[im]) );
            anglePtr = new H5::Attribute( datasetPtr->createAtribute( "phaseangle", H5::PredType::NATIVE_DOUBLE, alphas[im]) );
			delete datasetPtr;
            delete ilaserPtr;
            delete xrayPtr;
            delete poslaserPtr;
            delete posxrayPtr;
            delete delayPtr;
            delete anglePtr;
		}
		//std::cerr << "Made it past H5 image fill" << std::endl << std::flush;
		delete dataspace;
		delete paramgrp;
		delete dsetgrp;
		delete h5filePtr;
	} 


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
	//std::cerr << "Hmmmm, are we done yet?\n" << std::flush;


	return 0;
}


