#include "scan_material.hpp"
#include <cstdlib>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <complex>

#include <vector>
#include <random>
#include <chrono>

using namespace std;
using namespace Constants;


/* Here's main */
int main(int argc, char* argv[])
{
	/* start reading in from the command line, make and write fake data */

	size_t nimages = size_t(atoi(getenv("nimages")));

	//std::default_random_engine rng;
	std::random_device rng;
	std::normal_distribution<double> distribution(
			double(atof(getenv("delays_mean"))),
			double(atof(getenv("delays_std"))));
	std::uniform_real_distribution<double> uni_distribution(
			double(atof(getenv("delays_mean")))-double(atof(getenv("delays_std"))),
			double(atof(getenv("delays_mean")))+double(atof(getenv("delays_std"))));
	//std::chrono::system_clock::time_point tp = std::chrono::system_clock::now();
	//std::chrono::system_clock::duration dtm = tp.time_since_epoch();
	//rng.seed(dtm.count());


	/*
	HERE HERE HERE HERE
		Add a method to include slow drift for x-ray center (x,y) and delay center, and also the angle of the delay relative to the input coordinates.
		Include jitter in x-ray center (x,y) from shot to shot as well

		Also simulate various x-ray diameters.
	*/
	FiberBundle bundle(boost::lexical_cast<size_t>(atoi(getenv("nfibers"))));
	bundle.fiberdiameter(boost::lexical_cast<float>(atof(getenv("fiberdiam"))));
	bundle.laserdiameter(boost::lexical_cast<float>(atof(getenv("laserdiam"))));
	bundle.xraydiameter(boost::lexical_cast<float>(atof(getenv("xraydiam"))));
	bundle.set_fsPmm(boost::lexical_cast<float>(atof(getenv("bundle_fsPmm"))));
	bundle.scalePolarCoords();
	bundle.shuffle_output();
	bundle.Ixray(float(1.));
	bundle.Ilaser(float(1.));
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

	filename = filebase + "fibermap.out";
	std::ofstream mapfile(filename.c_str(),std::ios::out);
	bundle.print_mapping(mapfile);
	mapfile.close();


	float dalpha(float(atof(getenv("drifting_alpha")))*M_PI/nimages);


	bool doublepulse = ((int)atoi(getenv("doublepulse")) > 0 ? true : false);
	double doublepulsedelay(0.);
	if (doublepulse){
		doublepulsedelay = (double)atof( getenv("doublepulsedelay") ) ; // this one gets used directly in atomic units
	}
	double lambda0 = (double)( atof( getenv("lambda0") ) );
	double lambda_width = (double)( atof( getenv("lambda_width") ) );
	double lambda_onoff = (double)( atof( getenv("lambda_onoff") ) );
	double tspan = (double)( atof( getenv("tspan") ) )/fsPau<float>();
	double omega_low = 2*M_PI/(lambda0+lambda_width/2.0)*C_nmPfs<float>()*fsPau<float>();
	double omega_high = 2*M_PI/(lambda0-lambda_width/2.0)*C_nmPfs<float>()*fsPau<float>();
	double omega_width = omega_high-omega_low;
	double omega_onoff = 2*M_PI/(lambda0+lambda_width/2.0-lambda_onoff)*C_nmPfs<float>()*fsPau<float>() - omega_low;

	double omega0 = (omega_high+omega_low)/2.0;

	unsigned ngroupsteps = (unsigned)( atoi( getenv("ngroupsteps") ) );
	double groupdelay = (double)(atof(getenv("groupdelay"))); //this one gets used in fs
	double backdelay = (double)(atof(getenv("backdelay")));
	double groupstep = groupdelay/ngroupsteps ;
	double backstep = backdelay/ngroupsteps;
	const unsigned netalon = (unsigned)(atoi(getenv("netalon")));

	// work on the scan parameters: start with just delay //
	// ultimately, build a matrix, rows are delay, cols are intensity //

	/*
	   const double tstart=(double)( atof( getenv("tstart")));
	   const double tstep=(double)( atof( getenv("tstep")));
	   const double tstop=(double)( atof( getenv("tstop")));
	   const unsigned ntsteps=(unsigned)( (tstop - tstart) / tstep);
	 */

	/*
HERE: use nfibers to do this.  Hijack the tstart and tstep tstop for that.  This will break all other things so be careful.
	 */

	// now build scan parameters into an array of pulses... remember there is only need of one reference //

	// start with building reference and etalon replicas //
	PulseFreq masterpulse(omega0,omega_width,omega_onoff,tspan);
	std::vector<double> chirpvec(4,0.);
	std::vector<double> chirpnoisevec(4,0.);
	chirpvec[0] = (double)( atof( getenv("chirp") ) ) / std::pow(fsPau<float>(),int(2));// the difference in slopes at omega_low versus omega_high must equal tspan
	chirpvec[1] = (double)( atof( getenv("TOD") ) ) / std::pow(fsPau<float>(),int(3));
	chirpvec[2] = (double)( atof( getenv("FOD") ) ) / std::pow(fsPau<float>(),int(4));
	chirpvec[3] = (double)( atof( getenv("fifthOD") ) ) / std::pow(fsPau<float>(),int(5));


	if (atof(getenv("addrandomphase"))>0){
		//masterpulse.addrandomphase();
		std::string filename = filebase + "_spectralphaseFTpower.dat";
		std::ofstream outfile(filename,std::ios::out);
		masterpulse.print_phase_powerspectrum(outfile);
		outfile.close();
		filename = filebase + "_spectralphase.dat";
		outfile.open(filename);
		masterpulse.print_phase(outfile);
		outfile.close();
		filename = filebase + "_spectralamp.dat";
		outfile.open(filename);
		masterpulse.print_amp(outfile);
		outfile.close();
	}
	masterpulse.addchirp(chirpvec);							// chirp that ref pulse
	//HERE HERE HERE HERE 
	// print the power spectrum of the phase //
	PulseFreq * pulsearray[bundle.get_nfibers()];					// An array of pointers to PulseFreq objects
	PulseFreq * crosspulsearray[bundle.get_nfibers()];				// An array of pointers to PulseFreq objects


	MatResponse masterresponse(
			0.0,    // step delay doesn't matter in the reference, only the etalon relative to itself         // stepdelay
			(double)( atof( getenv("stepwidth") ) ),                                        // stepwidth
			((double)( atof( getenv("attenuation") ) ) - 1.0) / ngroupsteps + 1.0,          // attenuation
			(double)( atof( getenv("phase") ) ) / ngroupsteps                               // phase
			);
	masterresponse.aalphabbeta(
			(double)( atof( getenv("a") ) ),                // a
			(double)( atof( getenv("alpha" ) ) ),           // alpha
			(double)( atof( getenv("b") ) ),                //b
			(double)( atof( getenv("beta") ) )              // beta
			);
	masterresponse.setreflectance((double)(atof(getenv("etalon"))));
	masterresponse.setetalondelay((double)(atof(getenv("etalondelay"))));


	for (unsigned n=0;n<nimages;++n){ // outermost loop for nimages to produce //
		double t0 = uni_distribution(rng);
		std::cout << "Running for t0 = " << t0 << std::endl;

		bundle.Ixray(float(xray_distribution(rng)));
		bundle.Ilaser(float(laser_distribution(rng)));
		bundle.delay_angle(dalpha*float(n));
		bundle.center_Ixray(xray_pos_distribution(rng),xray_pos_distribution(rng));
		bundle.center_Ilaser(laser_pos_distribution(rng),laser_pos_distribution(rng));

		unsigned nthreads = (unsigned)atoi( getenv("nthreads") );
		unsigned i,tid,e,g;

		/*** Spawn a parallel region explicitly scoping all variables ***/
		double startdelay;
		bool refbuilt;
		refbuilt = false;
		MatResponse * oneresponsePtr = NULL;
		PulseFreq * onepulsedelayPtr = NULL;
		PulseFreq * crosspulsedelayPtr = NULL;

		PulseFreq pulseref(masterpulse);
		PulseFreq refdelay(masterpulse);

		for (e=0;e<netalon;e++){
			refdelay.delay(masterresponse.getetalondelay());
			refdelay.attenuate(pow(masterresponse.getreflectance(),(int)2));// change this to pow(R,2)
			refdelay.fft_totime();
			pulseref.fft_totime();
			pulseref += refdelay;
			pulseref.fft_tofreq();
			refdelay.fft_tofreq();
		}


#pragma omp parallel default(shared) private(tid,i,e,g,startdelay,oneresponsePtr,onepulsedelayPtr,crosspulsedelayPtr) num_threads(nthreads)
		{ // begin parallel region
			tid = omp_get_thread_num();

#pragma omp for // ordered 
			for(i = 0; i< bundle.get_nfibers(); i++){ // begin fibers loop
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

				//startdelay = (double)(atof( getenv("tstart")) + i*atof( getenv("tstep")));
				startdelay = boost::lexical_cast<double>(atof( getenv("tstart"))) + bundle.delay(i);

				startdelay += t0; // here ware are assuming a random jitter to t0
				pulsearray[i] = new PulseFreq(masterpulse);
				crosspulsearray[i] = new PulseFreq(masterpulse);
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
				oneresponsePtr = new MatResponse(
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

				onepulsedelayPtr = new PulseFreq(masterpulse);
				crosspulsedelayPtr = new PulseFreq(masterpulse);
				crosspulsedelayPtr->delay(interferedelay); // delay and attenuate in frequency domain

				onepulsedelayPtr->fft_totime();
				crosspulsedelayPtr->fft_totime();

				for(g=0;g<ngroupsteps;g++){ // begin groupsteps loop
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

				for (e=0;e<netalon;e++){ // begin etalon loop
					// back propagation step //
					startdelay -= (oneresponsePtr->getetalondelay()); // at front surface, x-rays see counter-propagating light from one full etalon delay
					for(g=0;g<ngroupsteps;g++){
						oneresponsePtr->setdelay(startdelay + g*backstep); // counterpropagating, x-rays work backwards through the optical
						oneresponsePtr->setstepvec_amp(onepulsedelayPtr);
						oneresponsePtr->setstepvec_phase(onepulsedelayPtr);
						oneresponsePtr->setstepvec_amp(crosspulsedelayPtr);
						oneresponsePtr->setstepvec_phase(crosspulsedelayPtr);
						if (doublepulse){
							oneresponsePtr->addstepvec_amp(onepulsedelayPtr,doublepulsedelay);
							oneresponsePtr->addstepvec_phase(onepulsedelayPtr,doublepulsedelay);
							oneresponsePtr->addstepvec_amp(crosspulsedelayPtr,doublepulsedelay);
							oneresponsePtr->addstepvec_phase(crosspulsedelayPtr,doublepulsedelay);
						}
						oneresponsePtr->buffervectors(onepulsedelayPtr); // this pulls down the tail of the response so vector is periodic on nsamples
						oneresponsePtr->buffervectors(crosspulsedelayPtr); // this pulls down the tail of the response so vector is periodic on nsamples
						onepulsedelayPtr->modulateamp_time();
						onepulsedelayPtr->modulatephase_time();
						crosspulsedelayPtr->modulateamp_time();
						crosspulsedelayPtr->modulatephase_time();
					}
					// forward propagation //
					for(g=0;g<ngroupsteps;g++){
						oneresponsePtr->setdelay(startdelay - g*groupstep); // forward propagating, x-rays advance on the optical
						oneresponsePtr->setstepvec_amp(onepulsedelayPtr);
						oneresponsePtr->setstepvec_phase(onepulsedelayPtr);
						oneresponsePtr->setstepvec_amp(crosspulsedelayPtr);
						oneresponsePtr->setstepvec_phase(crosspulsedelayPtr);
						if (doublepulse){
							oneresponsePtr->addstepvec_amp(onepulsedelayPtr,doublepulsedelay);
							oneresponsePtr->addstepvec_phase(onepulsedelayPtr,doublepulsedelay);
							oneresponsePtr->addstepvec_amp(crosspulsedelayPtr,doublepulsedelay);
							oneresponsePtr->addstepvec_phase(crosspulsedelayPtr,doublepulsedelay);

						}
						oneresponsePtr->buffervectors(onepulsedelayPtr); // this pulls down the tail of the response so vector is periodic on nsamples
						oneresponsePtr->buffervectors(crosspulsedelayPtr); // this pulls down the tail of the response so vector is periodic on nsamples
						onepulsedelayPtr->modulateamp_time();
						onepulsedelayPtr->modulatephase_time();
						crosspulsedelayPtr->modulateamp_time();
						crosspulsedelayPtr->modulatephase_time();
					}
					onepulsedelayPtr->fft_tofreq();
					crosspulsedelayPtr->fft_tofreq();
					onepulsedelayPtr->delay(oneresponsePtr->getetalondelay()); // delay and attenuate in frequency domain
					onepulsedelayPtr->attenuate(pow(oneresponsePtr->getreflectance(),(int)2));
					crosspulsedelayPtr->delay(oneresponsePtr->getetalondelay()); // delay and attenuate in frequency domain
					crosspulsedelayPtr->attenuate(pow(oneresponsePtr->getreflectance(),(int)2));
					onepulsedelayPtr->fft_totime();
					crosspulsedelayPtr->fft_totime();

					*(pulsearray[i]) += * onepulsedelayPtr;
					*(crosspulsearray[i]) += * crosspulsedelayPtr;

				} // end etalon loop

				pulsearray[i]->fft_tofreq();
				crosspulsearray[i]->fft_tofreq();

				pulsearray[i]->delay(interferedelay); // expects this in fs // time this back up to the crosspulse

				//#pragma omp flush(pulseref)

				delete onepulsedelayPtr;
				delete crosspulsedelayPtr;
				delete oneresponsePtr;
				onepulsedelayPtr = NULL;
				crosspulsedelayPtr = NULL;
				oneresponsePtr = NULL;
			} // end nfibers loop


		} // end parallel region


		pulseref.delay((double)atof(getenv("interferedelay"))); // expects this in fs
		pulseref.phase((double)atof(getenv("interferephase"))); // expects in units of pi

		filename = filebase + "interference.out." + std::to_string(n);
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

		// file for delay bins
		filename = filebase + "delaybins.out";
		ofstream outbins(filename.c_str(),ios::out); 
		for (i=0;i<bundle.get_nfibers();i++){
			outbins << bundle.delay(i) << "\n";
			*(pulsearray[i]) -= *(crosspulsearray[i]);
			*(pulsearray[i]) *= bundle.Ilaser(size_t(i)); 
			pulsearray[i]->appendwavelength(&interferestream);
		}
		//outstream.close();
		interferestream.close();
		outbins.close();

		for (unsigned i=0;i<bundle.get_nfibers();i++){
			delete pulsearray[i];
			delete crosspulsearray[i];
		}

		filename = filebase + "spectrum_xvals.out";
		std::ofstream outstream(filename.c_str(),ios::out);// app); use app to append delays to same file.
		pulseref.printfrequencybins(&outstream);
		outstream.close();


	} // outermost loop for nimages to produce //

	return 0;
}


