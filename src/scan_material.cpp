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
        << "\t\t===== " << std::asctime(std::localtime(&tstart)) << " ====\n"
        << "\t\t===== on host " << getenv("HOSTNAME") << " ====\n" 
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
        std::vector< uint16_t > calibdata(size_t(calibration.get_ndelays() * masterpulse.get_freqsamples()));
        std::vector<double> x;
        masterpulse.getLamVec(x);

#pragma omp parallel num_threads(nthreads) default(shared) shared(masterpulse,calpulsearray,calibdata,x)
        { // begin parallel region 1

            size_t tid = omp_get_thread_num();

            std::vector<double> y(x.size());
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
                calpulse = masterpulse;
                calcrosspulse = masterpulse;

                double startdelay(calibration.get_delay(d));
                calcrosspulse.delay(scanparams.interferedelay()); // delay in the frequency domain

                calpulse.fft_totime();
                calcrosspulse.fft_totime();

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

                    for(size_t g=0;g<scanparams.ngroupsteps();g++){
                        calibresponse.setdelay(etalondelay + g*scanparams.backstep()); 
                        // counterpropagating, x-rays work backwards through the optical

                        calibresponse.setstepvec_both_carriers(etalonpulse);
                        calibresponse.setstepvec_both_carriers(crossetalonpulse);
                        if (scanparams.doublepulse()){
                            calibresponse.addstepvec_both_carriers(etalonpulse,scanparams.doublepulsedelay());
                            calibresponse.addstepvec_both_carriers(crossetalonpulse,scanparams.doublepulsedelay());
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
                        calibresponse.setstepvec_both_carriers(etalonpulse);
                        calibresponse.setstepvec_both_carriers(crossetalonpulse);
                        if (scanparams.doublepulse()){
                            calibresponse.addstepvec_both_carriers(etalonpulse,scanparams.doublepulsedelay());
                            calibresponse.addstepvec_both_carriers(crossetalonpulse,scanparams.doublepulsedelay());
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
                    etalonpulse.delay(calibresponse.getetalondelay() * (1. + float(d)/float(calpulsearray.size())) ); // delay and attenuate in frequency domain
                    etalonpulse.attenuate(pow(calibresponse.getreflectance(),(int)2));
                    crossetalonpulse.delay(calibresponse.getetalondelay() * (1. + float(d)/float(calpulsearray.size())) ); // delay and attenuate in frequency domain
                    crossetalonpulse.attenuate(pow(calibresponse.getreflectance(),(int)2));
                    etalonpulse.fft_totime();
                    crossetalonpulse.fft_totime();
                    calpulse += etalonpulse;
                    calcrosspulse += crossetalonpulse;
                } // end etalon loop


                calpulse.fft_tofreq();
                calcrosspulse.fft_tofreq();
                calpulse.delay(scanparams.interferedelay()); // expects this in fs // time this back up to the crosspulse

                calpulse -= calcrosspulse;
                size_t f = calpulsearray.size()-d-1; // reversing order for sake of chirp calib matrix
                calpulsearray[f] = calpulse;

                calpulse.fillfrequency_bytes(x,
                        calpulse.getSigVec(y),
                        calibdata,
                        f);

            } // end of loop calibration.get_ndelays() to produce //


#pragma omp master
            {
                std::cout << "\t\t############ ending parallel region 1 ###########\n" << std::flush;
            }
        } // end parallel region 1

        if (bool(getenv("H5OUTPUT"))){
            H5::H5File * calibfilePtr;
            std::string calfilename = scanparams.calfilebase() + "interference.calibration.h5";
            calibfilePtr = new H5::H5File ( calfilename , H5F_ACC_TRUNC );
            size_t dims[1] = {calibdata.size()};
            size_t rank = 1;
            H5::DataSpace * dataspace = new H5::DataSpace( rank , dims ); 	


            std::vector< float > delayvec;
            calibration.fill_delays<float>(delayvec);

            H5::DataSet * calibsetPtr;
            std::string calname = "/calibration";
            calibsetPtr = new H5::DataSet( calibfilePtr->createDataSet( calname, H5::PredType::NATIVE_UINT16, *dataspace ) );
            calibsetPtr->write( calibdata.data(), H5::PredType::NATIVE_USHORT);

            dims[0] = delayvec.size();

            H5::DataSpace * delayspace = new H5::DataSpace( rank , dims ); 	
            std::string delname = "/delays";
            H5::DataSet * delaysetPtr;
            delaysetPtr = new H5::DataSet( calibfilePtr->createDataSet( delname, H5::PredType::NATIVE_FLOAT, *delayspace ) );
            delaysetPtr->write( delayvec.data(), H5::PredType::NATIVE_FLOAT );
            delete delayspace;
            delete delaysetPtr;

            uint16_t calimshape[2] = {calibration.get_ndelays() , masterpulse.get_freqsamples()};
            dims[0] = 2;
            H5::DataSpace * shapespace = new H5::DataSpace( 1 , dims );
            H5::Attribute * calshapePtr = new H5::Attribute( calibsetPtr->createAttribute( "calshape", H5::PredType::NATIVE_UINT16, *shapespace) );
            calshapePtr->write(H5::PredType::NATIVE_UINT16,calimshape);
            delete calshapePtr;
            delete shapespace;

            dims[0] = 4; 
            H5::DataSpace * chirpspace = new H5::DataSpace( rank , dims );
            H5::Attribute * chirpPtr = new H5::Attribute( calibsetPtr->createAttribute( "chirp", H5::PredType::NATIVE_FLOAT, *chirpspace) );
            std::vector<float> chirpvec = { (float) atof( getenv("chirp") ) 
                , (float) atof( getenv("TOD") )
                , (float) atof( "FOD" )
                , (float) atof( "fifthOD" )};
            chirpPtr->write(H5::PredType::NATIVE_FLOAT,chirpvec.data());
            delete chirpPtr;
            delete chirpspace;


            delete calibsetPtr;
            delete calibfilePtr;
            delete dataspace;
            std::cout << "Finished with the h5 file for calibration image/matrix\n" << std::flush;
        } else {
            std::cout << "|\t done with calibration delays\n" << std::flush;
            std::string calfilename = scanparams.calfilebase() + "interference.calibration";
            std::string derivfilename = scanparams.calfilebase() + "interference.calibration.derivative";
            std::string calfilename_delays = scanparams.calfilebase() + "interference.calibration.delays";
            std::string calfilename_wavelengths = scanparams.calfilebase() + "interference.calibration.wavelengths";
            ofstream calibrationstream(calfilename.c_str(),ios::out); 
            ofstream derivstream(derivfilename.c_str(),ios::out); 
            ofstream calibrationstream_delays(calfilename_delays.c_str(),ios::out); 
            ofstream calibrationstream_wavelengths(calfilename_wavelengths.c_str(),ios::out); 
            calibrationstream << "# wavelengths\n#";
            calpulsearray[0].printwavelengthbins(&calibrationstream);
            calpulsearray[0].printwavelengthbins(&derivstream);
            calpulsearray[0].printwavelengthbins(&calibrationstream_wavelengths);
            calibrationstream << "# delays\n#";
            calibrationstream_delays << "# delays\n";
            for (size_t i = 0 ; i< calibration.get_ndelays(); ++i){
                calibrationstream << calibration.get_delay(i) << "\t";
                calibrationstream_delays << calibration.get_delay(i) << "\t";
            }
            calibrationstream << "\n";
            calibrationstream_delays << "\n";

            for (size_t n=0;n<calpulsearray.size();++n){
                calpulsearray[n].appendwavelength(&calibrationstream);
            }

            calibrationstream.close();
            derivstream.close();
            calibrationstream_delays.close();
            calibrationstream_wavelengths.close();
            //bin_calibrationstream.close();
            std::cout << "Finished with the calibration image/matrix\n" << std::flush;
        }



    } // end if (!getenv("skipcalibration"))


    //############################################//
    //############## Images section ##############//
    //############################################//


    auto localnow = std::chrono::system_clock::now();
    std::time_t ttime = std::chrono::system_clock::to_time_t(localnow);
    std::cout << "chrono localtime = " << std::ctime(& ttime) << std::endl;
    std::tm * local_time = std::localtime(& ttime);


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
            //std::cout << "\t\t############ entering parallel image: tid = " << int(tid) << " ###########\n" << std::flush;


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
#pragma omp master
            {
                std::cout << "========================================================================="
                    <<   "\n\t\t ==== http://www.fftw.org/fftw3_doc/Advanced-Complex-DFTs.html ===="
                    <<   "\n\t\t ====         use this for defining multiple fibers as         ===="
                    <<   "\n\t\t ====         contiguous blocks for row-wise FFT as 2D         ===="
                    <<   "\n\t\t ==================================================================\n" << std::flush;
            }

            for (size_t n=0;n<scanparams.nimages();++n)
            { // outermost loop for nimages to produce //

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


                //DebugOps::pushout(std::string("Running image " + std::to_string(n) + " for t0 = " + std::to_string(t0) + " in threaded for loop, thread " + std::to_string(tid)));
                if (bool(getenv("printASCIIimages"))){
                    std::string mapfilename = scanparams.filebase() + "fibermap.out." + std::to_string(n);
                    //std::cout << "fibermap file = " << mapfilename << std::endl << std::flush;
                    std::ofstream mapfile(mapfilename.c_str(),std::ios::out);
                    if (!parabundle.print_mapping(mapfile,t0))
                        std::cerr << "Something failed in printing this fibermapping out\n" << std::flush;
                    mapfile.close();
                }


                std::vector<double> x,y;
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
                } // end nfibers loop

                std::complex<double> z_laser = parabundle.center_Ilaser();
                std::complex<double> z_xray = parabundle.center_Ixray();

                laserposition[scanparams.nimages()*tid + n]=std::array<double,2>({z_laser.real(),z_laser.imag()});
                xrayposition[scanparams.nimages()*tid + n]=std::array<double,2>({z_xray.real(),z_xray.imag()});
                delays[scanparams.nimages()*tid + n] = double(t0);
                alphas[scanparams.nimages()*tid + n] = double(parabundle.delay_angle());
                ilaser[scanparams.nimages()*tid + n] = double(parabundle.Ilaser());
                ixray[scanparams.nimages()*tid + n] = double(parabundle.Ixray());

                if (bool(getenv("printASCIIimages"))){
                    ofstream interferestream;
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
        std::cout << "Opening h5 file for saving:\t" << imfilename << std::endl << std::flush;
        h5filePtr = new H5::H5File ( imfilename , H5F_ACC_TRUNC );
        std::string name = "/params";
        paramgrp = new H5::Group( h5filePtr->createGroup( name ) );
        std::string dsetname = "/datasets";
        dsetgrp = new H5::Group( h5filePtr->createGroup( dsetname ) );


        const int rank(1);
        size_t imshape[2] = {masterbundle.get_nfibers(), masterpulse.get_freqsamples()};
        size_t flatimgsize = std::accumulate( imshape, imshape+1, 1LL, std::multiplies<size_t>() ); 
                                                                                                    
        size_t dims[rank] = {flatimgsize};
        size_t posdims[1] = {2};
        size_t shapedims[1] = {2};
        size_t chirpdims[1] = {4};
        H5::StrType vls_type(0, H5T_VARIABLE); // variable length string type in H5
        H5::DataSpace * dataspace = new H5::DataSpace( rank , dims ); 	
        H5::DataSpace * scalarspace = new H5::DataSpace( H5S_SCALAR );
        H5::DataSpace * tuplespace = new H5::DataSpace( 1 , posdims );
        H5::DataSpace * shapespace = new H5::DataSpace( 1 , shapedims );
        H5::DataSpace * chirpspace = new H5::DataSpace( 1 , chirpdims );

        H5::DataSpace * stringspace = new H5::DataSpace( H5S_SCALAR );
        H5::DataSpace * floatspace = new H5::DataSpace( H5S_SCALAR );
        H5::DataSpace * doublespace = new H5::DataSpace( H5S_SCALAR );
        H5::DataSpace * uintspace = new H5::DataSpace( H5S_SCALAR );

        H5::Attribute * imshapePtr = new H5::Attribute( paramgrp->createAttribute( "imshape", H5::PredType::NATIVE_UINT16, *shapespace) );
        imshapePtr->write(H5::PredType::NATIVE_UINT16,imshape);

        H5::Attribute * chirpPtr = new H5::Attribute( paramgrp->createAttribute( "chirp", H5::PredType::NATIVE_FLOAT, *chirpspace) );
        std::vector<float> chirpvec = { (float) atof( getenv("chirp") ) 
            , (float) atof( getenv("TOD") )
                , (float) atof( "FOD" )
                , (float) atof( "fifthOD" )};
        chirpPtr->write(H5::PredType::NATIVE_FLOAT,chirpvec.data());
        delete chirpPtr;
        delete chirpspace;

        H5::Attribute * carriernamePtr;
        carriernamePtr = new H5::Attribute( paramgrp->createAttribute( "carriersfile" , vls_type, *stringspace) );
        std::string s;
        if (bool(getenv("usediamond"))){
            s = std::string("UsingDiamond_no_carriers_file_only_xray_photonenergy");
        } else {
            s = std::string(getenv("carriersfile"));
        }
        carriernamePtr->write(vls_type, s);
        delete carriernamePtr;

        H5::Attribute * calnamePtr = new H5::Attribute( paramgrp->createAttribute( "calibfile" , vls_type, *stringspace) );
        s = std::string(getenv("calfilebase"));
        calnamePtr->write(vls_type, s);
        delete calnamePtr;

        s = std::string( getenv("HOSTNAME") );
        H5::Attribute * hostPtr = new H5::Attribute( paramgrp->createAttribute( "hostname", vls_type, *stringspace ) );
        hostPtr->write( vls_type, s );
        delete hostPtr;

        float flt;
        uint8_t uint;

        H5::Attribute * floatPtr;
        floatPtr = new H5::Attribute( paramgrp->createAttribute( "xray_photonenergy", H5::PredType::NATIVE_FLOAT, *floatspace ) );
        flt = (float) xrayphoton_energy;
        floatPtr->write( H5::PredType::NATIVE_FLOAT, &flt);
        delete floatPtr;
        floatPtr = new H5::Attribute( paramgrp->createAttribute( "interferedelay" , H5::PredType::NATIVE_FLOAT, *floatspace ) );
        flt = (float) scanparams.interferedelay();
        floatPtr->write( H5::PredType::NATIVE_FLOAT, &flt);
        delete floatPtr;
        floatPtr = new H5::Attribute( paramgrp->createAttribute( "interferephase" , H5::PredType::NATIVE_FLOAT, *floatspace ) );
        flt = (float) scanparams.interferephase();
        floatPtr->write( H5::PredType::NATIVE_FLOAT, &flt);
        delete floatPtr;

        H5::Attribute * netalonPtr;
        netalonPtr = new H5::Attribute( paramgrp->createAttribute( "netalon", H5::PredType::NATIVE_UINT8, *uintspace ) );
        uint = (uint8_t) scanparams.netalon();
        netalonPtr->write( H5::PredType::NATIVE_UINT8, &uint);
        delete netalonPtr;

        // fibermap as dataset inside /params
        std::vector<std::vector<float>> map;
        masterbundle.get_map<float>(map);
        std::vector<uint8_t> keys;
        masterbundle.get_keys<uint8_t>(keys);

        size_t keydims[1] = {keys.size()};
        size_t mapdims[2] = {map.size(),map[0].size()};
        H5::DataSpace * keyspace = new H5::DataSpace( 1 , keydims );
        H5::DataSpace * mapspace = new H5::DataSpace( 2 , mapdims );
        H5::DataSet * fibermapPtr = new H5::DataSet( paramgrp->createDataSet( "/params/fibermap", H5::PredType::NATIVE_FLOAT, *mapspace ) );
        fibermapPtr->write( map.data(), H5::PredType::NATIVE_FLOAT );
        H5::DataSet * fiberkeysPtr = new H5::DataSet( paramgrp->createDataSet( "/params/fiberkeys", H5::PredType::NATIVE_UINT16, *keyspace ) );
        fiberkeysPtr->write( keys.data(), H5::PredType::NATIVE_UINT16 );
        delete fibermapPtr;
        delete fiberkeysPtr;
        delete keyspace;
        delete mapspace;


        delete imshapePtr;
        delete shapespace;
        delete stringspace;
        delete floatspace;
        delete doublespace;
        delete uintspace;

        std::cerr << "Made it here in H5 output\n" << std::flush;
        for (size_t im=0;im<datablock.size();im++){
            H5::DataSet * datasetPtr;
            std::string imname = "/datasets/im_" + std::to_string((int)im);
            datasetPtr = new H5::DataSet( dsetgrp->createDataSet( imname, H5::PredType::NATIVE_UINT16, *dataspace ) );
            datasetPtr->write( datablock[im].data(), H5::PredType::NATIVE_USHORT);

            H5::Attribute * ilaserPtr = new H5::Attribute( datasetPtr->createAttribute( "Ilaser", H5::PredType::NATIVE_DOUBLE, *scalarspace) );
            ilaserPtr->write(H5::PredType::NATIVE_DOUBLE, &(ilaser[im]) );
            H5::Attribute * ixrayPtr = new H5::Attribute( datasetPtr->createAttribute( "Ixray", H5::PredType::NATIVE_DOUBLE, *scalarspace) );
            ixrayPtr->write(H5::PredType::NATIVE_DOUBLE, &(ixray[im]) );
            H5::Attribute * delayPtr = new H5::Attribute( datasetPtr->createAttribute( "delay", H5::PredType::NATIVE_DOUBLE, *scalarspace) );
            delayPtr->write(H5::PredType::NATIVE_DOUBLE, &(delays[im]) );
            H5::Attribute * anglePtr = new H5::Attribute( datasetPtr->createAttribute( "phaseangle", H5::PredType::NATIVE_DOUBLE, *scalarspace) );
            anglePtr->write(H5::PredType::NATIVE_DOUBLE, &(alphas[im]) );

            H5::Attribute * poslaserPtr = new H5::Attribute( datasetPtr->createAttribute( "positionlaser", H5::PredType::NATIVE_DOUBLE, *tuplespace) );
            poslaserPtr->write(H5::PredType::NATIVE_DOUBLE,laserposition[im].data());
            H5::Attribute * posxrayPtr = new H5::Attribute( datasetPtr->createAttribute( "positionxray", H5::PredType::NATIVE_DOUBLE, *tuplespace) );
            posxrayPtr->write(H5::PredType::NATIVE_DOUBLE,xrayposition[im].data());

            delete datasetPtr;

            delete ilaserPtr;
            delete ixrayPtr;
            delete delayPtr;
            delete anglePtr;
            delete poslaserPtr;
            delete posxrayPtr;
            // The images of e.g. the xray and laser intensity and individual delay can be used in the future for image specific material/delay/shadow flutuations
        }
        //std::cerr << "Made it past H5 image fill" << std::endl << std::flush;
        delete dataspace;

        delete scalarspace;
        delete tuplespace;

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
        << "\t\t======= scan_material stopped ========\n"
        << "\t\t===== " << std::asctime(std::localtime(&tstop)) << " ====\n"
        << "\t\t===== in " << tstop << " s ===========\n"
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


