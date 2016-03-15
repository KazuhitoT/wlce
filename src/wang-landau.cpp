#include <iomanip>
#include <iterator>
#include <numeric>
#include <random>
#include <chrono>

#include "./Parser.hpp"
#include "./Input.hpp"
#include "./WLconf.hpp"

bool checkHistogramFlat(const std::vector<double>& histogram, const std::vector<int>& index_neglect_bin, double lflat, double low_cutoff=0.5){
	double ave = accumulate(histogram.begin(), histogram.end(), 0) / (double)histogram.size();
	double limit = ave * lflat ;
	for(int i=0, imax=histogram.size(); i<imax; ++i){

		auto it = find( index_neglect_bin.begin(), index_neglect_bin.end() , i);
		if( it != index_neglect_bin.end() ) continue;

		if(histogram.at(i) < (ave*(1.0-low_cutoff))) continue;
		else if(histogram.at(i) < limit) return false;
	}
	return true;
}


int main(int argc, char* argv[]){
	std::shared_ptr<Input> in(new Input("wang-landau.ini"));

	int mcstep, bin, flatcheck_step;
	double logfactor, logflimit, emin, emax, flat_criterion;

	in->setData("BIN",  bin, true);
	in->setData("MCSTEP", mcstep, true);

	in->setData("FLATCHECKSTEP", flatcheck_step, true);
	in->setData("LOGFACTOR", logfactor, true);
	in->setData("LOGFLIMIT", logflimit, true);
	in->setData("FLATCRITERION", flat_criterion, true);

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

	WLconf PoscarSpin("./poscar.spin", in, cluster_in.getCluster(), ecicar.getEci(), nullptr, nullptr, spin_exchange);
	PoscarSpin.dispCorr();

	const int N = PoscarSpin.getSpins().size();

	std::vector<int> index_neglect_bin = PoscarSpin.getNeglectBinIndex();

	std::cout << "----------  initial correlation functions" << std::endl;
	PoscarSpin.dispCorr();
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << "----------  Information about [wang-landau.ini]" << std::endl;
	PoscarSpin.dispInput();

	if( index_neglect_bin.size() ){
		std::cout << "NEGLECT BIN INDEX :" << std::endl;
		std::cout << " ";
		for(auto i : index_neglect_bin) std::cout << i << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << std::endl;

	auto start = std::chrono::system_clock::now();

	std::cout << "----------  Log  ----------" << std::endl;
	std::cout << std::endl;

	std::vector<double>  dos(bin, 0);


	unsigned int fstep = 0;
	unsigned int tstep = 0;
	unsigned int mc_time = 0;
	bool isConverged = false;

	start = std::chrono::system_clock::now();

	auto f_start = std::chrono::system_clock::now();
	auto f_end   = std::chrono::system_clock::now();
	//
 // 	int final_mcsweep = 0;
 // 	while(logflimit < logfactor){
	// 	std::vector<double> histogram(bin, 0);
	// 	std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
	// 	std::cout << std::setprecision(10) << "--- fstep " << fstep << " --  log(factor) = " << logfactor << std::endl;
	// 	for(int i=1; ; ++i, ++mc_time){
	// 		for(int j=0; j<N ; ++j){  // スピン数だけステップ回す これで1MCSweep
	//
	// 			int before_index = PoscarSpin.getIndex();
	//
	// 			wl_step(PoscarSpin, dos);
	// 			int index = PoscarSpin.getIndex();
	// 			dos.at(index) += factor;
	// 			histogram.at(index) += 1;
	//
	//
	// 		}  /*  end MC sweep */
	//
	// 		if((i % mcstep) == 0){
	// 			outputHistogram(dos,  histogram, emin, edelta, fstep);
	// 			std::cout << i << "sweep done " << std::endl;
	// 		}
	// 		if((i % in.flatcheck) == 0 and checkHistogramFlat(histogram, in.lflat, in.low_cutoff)){
	// 			outputHistogram(dos,  histogram, in.emin, in.edelta, fstep);
	// 			final_mcsweep = i;
	// 			break;
	// 		}
	// 	}
	// 	f_end = std::chrono::system_clock::now();
	// 	auto f_sec = std::chrono::duration_cast<std::chrono::seconds>(f_end-f_start).count();
	// 	std::cout << " fstep = " << fstep << " " << f_sec << " seconds "<< final_mcsweep << " sweeps" << std::endl;
	// 	f_start = std::chrono::system_clock::now();
	// 	factor/=2.0;
	// 	++fstep;
	// }

	auto end = std::chrono::system_clock::now();
	auto sec = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
	std::cout << "Time for calculating the DOS is " << sec << " seconds. "<< std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "check final correlation functions" << std::endl;
	PoscarSpin.dispCorr();
	std::cout << "---------------------------" << std::endl;
	PoscarSpin.setInitialCorrelationFunction();
	PoscarSpin.dispCorr();

	return 0;
}
