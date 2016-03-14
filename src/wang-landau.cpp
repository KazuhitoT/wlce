#include <iomanip>
#include <iterator>
#include <numeric>
#include <random>
#include <chrono>

#include "./Parser.hpp"
#include "./Input.hpp"
#include "./WLconf.hpp"

// #include "./WLconf.hpp"
// #include "update.hpp"
// #include "wlInput.hpp"
// #include "output.hpp"
// #include "wlStep.hpp"
// #include "enum.hpp"

bool checkHistogramFlat(const std::vector<double>& histogram, const std::vector<int>& index_neglect_bin, double lflat, double low_cutoff=0.5){
	double ave = accumulate(histogram.begin(), histogram.end(), 0) / (double)histogram.size();
	double limit = ave * lflat ;
	for(int i=0, imax=histogram.size(); i<imax; ++i){
		/*  存在しないbinはスルーする */
		auto it = find( index_neglect_bin.begin(), index_neglect_bin.end() , i);
		if( it != index_neglect_bin.end() ) continue;

		if(histogram.at(i) < (ave*(1.0-low_cutoff))) continue;
		else if(histogram.at(i) < limit) return false;
	}
	return true;
}


int main(int argc, char* argv[]){
	Input in("wang-landau.ini");

	std::vector<double> spince, spinposcar;
	std::vector<double> chemical_potential;
	int mcstep, setrandom, bin, flatcheck;
	double factor, flimit, emin, emax, lflat;
	std::string input_spin_filename;

	in.setData("SPINCE", spince, true);
	in.setData("SPINPOSCAR", spinposcar, true);
	in.setData("CHEMIPOT", chemical_potential);
	in.setData("MCSTEP", mcstep, true);
	in.setData("SPININPUT", input_spin_filename);
	in.setData("SETRANDOM", setrandom);

	in.setData("BIN", bin, true);
	in.setData("FLATCHECK", flatcheck, true);
	in.setData("FACTOR", factor, true);
	in.setData("FLIMIT", flimit, true);
	in.setData("EMIN", emin, true);
	in.setData("EMAX", emax, true);
	in.setData("LFLAT", lflat, true);

	bool spin_exchange = false;
	for(int i=1; i<argc; i++) {
		std::string str(argv[i]);
		if( str == "-exchange" ){
			spin_exchange = true;
		} else {
			std::cerr << " ERROR : invalid commandline argument [" << str << "]" << std::endl;
			exit(1);
		}
	}

	const ParseEcicar ecicar("./ecicar");
	const ParseMultiplicityIn multiplicity_in("./multiplicity.in", ecicar.getIndex());
	const ParseClusterIn  cluster_in("./clusters.in", ecicar.getIndex(), multiplicity_in.getMultiplicityIn());

	WLconf PoscarSpin("./poscar.spin",
			spinposcar,
			spince,
			cluster_in.getCluster(),
			ecicar.getEci(),
			nullptr,
			nullptr,
			spin_exchange,
			chemical_potential
	);

	PoscarSpin.dispCorr();

//
// 	const int N = PoscarSpin.getSpins().size();
// 	std::vector<int>     index_neglect_bin;
//
// 	if( in.setrandom and in.input_spin_filename.size() ){
// 		std::cout << "ERROR : Both SETRANDOM and SPININPUT are available." << std::endl;
// 		std::cout << "Plese make available only one configuration." << std::endl;
// 		exit(1);
// 	} else if(in.input_spin_filename.size()){
// 		PoscarSpin.setSpinsFromDat(in.input_spin_filename, in.emin, in.emax);
// 		index_neglect_bin = PoscarSpin.getNeglectBinIndex();
// 		PoscarSpin.setCorrelationFunctionFromClucar();
// 	} else if(in.setrandom){
// 		PoscarSpin.setSpinsRandom();
// 	}
// 	PoscarSpin.setIndex();
//
// 	double delta = ( in.emax - in.emin )/ (double)(in.bin);
//
// 	std::cout << "----------  initial correlation functions" << std::endl;
// 	PoscarSpin.dispCorr();
// 	std::cout << std::endl;
// 	std::cout << std::endl;
//
// 	#ifdef DEBUG_SPIN_DISP
// 	std::cout << "----------  initial spins set (DEBUG)" << std::endl;
// 	PoscarSpin.dispSpin();
// 	std::cout << std::endl;
// 	std::cout << std::endl;
// 	#endif
//
// 	std::cout << "----------  Information about [wang-landau.ini]" << std::endl;
// 	in.disp();
//
// 	if( index_neglect_bin.size() ){
// 		std::cout << "NEGLECT BIN INDEX :" << std::endl;
// 		std::cout << " ";
// 		for(auto i : index_neglect_bin) std::cout << i << " ";
// 		std::cout << std::endl;
// 	}
// 	std::cout << std::endl;
// 	std::cout << std::endl;
//
// 	auto start = std::chrono::system_clock::now();
//
// 	if(in.macromode == MAHALANOBIS){
// 		PoscarSpin.setCovarianceMatrix(in.input_spin_filename);
// 		auto end = std::chrono::system_clock::now();
// 		auto sec = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
// 		std::cout << "Time for setting covariance matrix is " << sec << " seconds. "<< std::endl;
// 		std::cout << std::endl;
// 	} else if(in.macromode == MAHALANOBIS_FOR_ALLBIN){
// 		std::cout << "This mode is not be implemented."<< std::endl;
// 		exit(1);
// 	}
//
// 	std::cout << "----------  Log  ----------" << std::endl;
// 	std::cout << std::endl;
//
// 	std::vector<double>  dos(in.bin, 0);
// 	std::vector<std::vector<double>> corr_dos_2d = std::vector<std::vector<double>>(in.bin, std::vector<double>(in.bin, 0));
// 	std::vector<std::vector<double>> corr_histograms = std::vector<std::vector<double>>(in.corr_histogram_vec_index.size(), std::vector<double>(in.bin, 0));
// 	// std::vector<std::vector<double>> transition_matrix = std::vector<std::vector<double>>(in.bin, std::vector<double>(in.bin, 0));
//
// 	random_device rnd;
// 	mt19937 mt(rnd());
// 	uniform_int_distribution<> randN(0, N-1);
//
// 	double factor = in.factor;
// 	unsigned int fstep = 0;
// 	unsigned int tstep = 0;
// 	unsigned int mc_time = 0;
// 	bool isConverged = false;
//
// 	if( in.profile ){
// 		ofstream ofs_profile("profile.out");
// 		ofs_profile.setf(std::ios_base::fixed, std::ios_base::floatfield);
// 		ofs_profile.close();
// 	}
//
// 	if( in.corr_profile ){
// 		ofstream ofs_corr_profile("corr-profile.out");
// 		ofs_corr_profile.setf(std::ios_base::fixed, std::ios_base::floatfield);
// 		ofs_corr_profile.close();
// 	}
//
// 	ofstream ofs_tt("tunneling_time.out");
// 	ofs_tt.setf(std::ios_base::fixed, std::ios_base::floatfield);
// 	ofs_tt.close();
//
// 	start = std::chrono::system_clock::now();
//
// 	auto f_start = std::chrono::system_clock::now();
// 	auto f_end   = std::chrono::system_clock::now();
//
// 	/*
// 	 * tunneling_time_flag = -1 -> attained lowest energy and searching highest energy state
//  	 * tunneling_time_flag = 0  -> not attained lowest nor highest energy state
// 	 * tunneling_time_flag = 1  -> attained highest energy and searching lowest energy state
// 	 */
//
// 	int tunneling_time_flag = 0;
// 	unsigned int tmp_tunneling_time_step = 0;
//
//  	int final_mcsweep = 0;
//  	while(in.flimit < factor){
// 		std::vector<double> histogram(in.bin, 0);
// 		std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
// 		std::cout << setprecision(10) << "--- fstep " << fstep << " --  log(factor) = " << factor << std::endl;
// 		for(int i=1; ; ++i, ++mc_time){
// 			for(int j=0; j<N ; ++j){  // スピン数だけステップ回す これで1MCSweep
//
// 				int before_index = PoscarSpin.getIndex();
//
// 				if(in.method == CONTINUOUS){
// 					double cindex = wl_step(PoscarSpin, dos, in);
// 					updateDos(dos, PoscarSpin, factor, in);
// 					updateDos(histogram, cindex, 1, in);
// 				} else if(in.method == DISCRETE){
// 					wl_step(PoscarSpin, dos, in);
// 					int index = PoscarSpin.getIndex();
// 					dos.at(index) += factor;
// 					histogram.at(index) += 1;
// 				}
//
// 				int after_index = PoscarSpin.getIndex();
//
// 				// transition_matrix[before_index][after_index] += 1.0;
//
// 				double ene = PoscarSpin.getTotalEnergy();
// 				if( ene <= in.lowest_energy_tt ){
// 					if( tunneling_time_flag == 1 ) {
// 						outputTunnelingTime(fstep, mc_time, mc_time-tmp_tunneling_time_step);
// 					}
// 					tunneling_time_flag = -1;
// 					tmp_tunneling_time_step = mc_time;
// 				} else if( ene >= in.highest_energy_tt ) {
// 					if( tunneling_time_flag == -1 ) {
// 						outputTunnelingTime(fstep, mc_time, mc_time-tmp_tunneling_time_step);
// 					}
// 					tunneling_time_flag = 1;
// 					tmp_tunneling_time_step = mc_time;
// 				}
//
// 				if( in.corr2ddos_vec_index.size() ) updateCorrDos(corr_dos_2d, PoscarSpin, in);
//
// 				if(corr_histograms.size()){
// 					for(int k=0; k<corr_histograms.size(); ++k){
// 						exit(1);
// 						// updateHistogram(corr_histograms.at(k), PoscarSpin.getCorrelationFunctions().at(in.corr_histogram_vec_index.at(k)), factor, in.kwidth, in.cutoff, in.bin);
// 					}
// 				}
//
// 			}  /*  end MC sweep */
// 			if( in.profile )     PoscarSpin.outputCorrSpin(mc_time, ENERGY, "profile.out");
// 			if( in.output )      PoscarSpin.outputCorrSpin(mc_time, in.output);
// 			if( in.corr_profile )PoscarSpin.outputCorr();
//
// 			if((i % in.mcstep) == 0){
// 				formatDos(dos);
// 				// outputTransitionMatrix(transition_matrix);
// 				outputHistogram(dos,  histogram, in.emin, in.edelta, fstep);
// 				#ifdef DEBUG_CORR_ANIMATION
// 				outputHistogram(dos,  histogram, in.emin, in.edelta, fstep, "-animation"+to_string(i/in.mcstep)+"-");
// 				#endif
// 				if(corr_histograms.size())  outputCorrHistogram(corr_histograms, in.bin, fstep);
// 				if( in.corr2ddos_vec_index.size() ) outputCorr2DDos(corr_dos_2d, in.bin, fstep, "-2d-step-"+to_string(i/in.mcstep)+"-f");
// 				std::cout << i << "sweep done " << std::endl;
// 			}
// 			if((i % in.flatcheck) == 0 and (
// 				(   in.method==DISCRETE   and checkHistogramFlat(histogram, in.lflat, in.low_cutoff) )
// 				or (in.method==CONTINUOUS and in.criterion == FLAT_CRITERION and checkHistogramFlat(histogram, index_neglect_bin, in.lflat, in.low_cutoff) )
// 				or (in.method==CONTINUOUS and in.criterion == MIN_STEP_CRITERION and checkHistogramFlat(histogram, index_neglect_bin, in.lflatstep) )
// 				)){
// 				outputHistogram(dos,  histogram, in.emin, in.edelta, fstep);
// 				final_mcsweep = i;
// 				break;
// 			}
// 		}
// 		f_end = std::chrono::system_clock::now();
// 		auto f_sec = std::chrono::duration_cast<std::chrono::seconds>(f_end-f_start).count();
// 		std::cout << " fstep = " << fstep << " " << f_sec << " seconds "<< final_mcsweep << " sweeps" << std::endl;
// 		f_start = std::chrono::system_clock::now();
// 		factor/=2.0;
// 		++fstep;
// 	}
//
// 	auto end = std::chrono::system_clock::now();
// 	auto sec = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
// 	std::cout << "Time for calculating the DOS is " << sec << " seconds. "<< std::endl;
// 	std::cout << std::endl;
// 	std::cout << std::endl;
// 	std::cout << "check final correlation functions" << std::endl;
// 	PoscarSpin.dispCorr();
// 	std::cout << "---------------------------" << std::endl;
// 	PoscarSpin.setCorrelationFunctionFromClucar();
// 	PoscarSpin.dispCorr();
//
// 	return 0;
}
