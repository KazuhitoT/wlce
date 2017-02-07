#include <iomanip>
#include <iterator>
#include <numeric>
#include <random>
#include <chrono>

#include "./parser.hpp"
#include "./input.hpp"
#include "./wlconfT.hpp"

constexpr double kb = 0.0000861733; // [eV/K]

void dispInput(std::shared_ptr<Input> in){
	//  [INPUT]
	std::cout << "TMIN          : " << in->getDataByString("TMIN") << std::endl;
	std::cout << "TMAX          : " << in->getDataByString("TMAX") << std::endl;
	std::cout << "TMOVELIMIT    : " << in->getDataByString("TMOVELIMIT") << std::endl;
	std::cout << "BIN           : " << in->getDataByString("BIN")  << std::endl;
	std::cout << "MCSTEP        : " << in->getDataByString("MCSTEP")        << std::endl;
	std::cout << "FLATCHECKSTEP : " << in->getDataByString("FLATCHECKSTEP") << std::endl;
	std::cout << "LOGFACTOR     : " << in->getDataByString("LOGFACTOR")     << std::endl;
	std::cout << "FLATCRITERION : " << in->getDataByString("FLATCRITERION") << std::endl;

	//  [OPTION]
	std::cout << "CHEMIPOT      : " << in->getDataByString("CHEMIPOT")  << std::endl;
	std::cout << "SPININPUT     : " << in->getDataByString("SPININPUT") << std::endl;
	std::cout << "SETRANDOM     : " << in->getDataByString("SETRANDOM") << std::endl;
}

double wl_step_by_T(WLconfT::WLconfT& conf, const std::vector<double>& dos){

	conf.setTotalEnergy();
	double energy = conf.getTotalEnergy();
	int index_before = conf.getIndex();
	int index_after  = index_before;
	double temperature_before = conf.calcTemperature();

	conf.setNewConfByTemperature();
	index_after  = conf.getIndex();
	double temperature_after  = conf.calcTemperature();

	double b_wl = exp( - energy*(1.0/temperature_after-1.0/temperature_before)/kb + dos.at(index_before)-dos.at(index_after));
	if(b_wl<conf.RandReal()){
		conf.Memento();
	}
	return conf.getTotalEnergy();
}

double wl_step_by_E(WLconfT::WLconfT& conf, const std::vector<double>& dos){

	conf.setTotalEnergy();
	double ene_before = conf.getTotalEnergy();

	conf.setNewConfByEnergy();
	double ene_after  = conf.getTotalEnergy();

	double b_metropolis = exp( -(ene_after-ene_before)/kb/conf.calcTemperature() );
	if(b_metropolis<conf.RandReal()){
		conf.Memento();
	}

	return conf.getTotalEnergy();
}


int main(int argc, char* argv[]){
	std::shared_ptr<Input> in(new Input("directZ.ini"));

	int mcstep, bin, flatcheck_step, minstep=0;
	double logfactor, logflimit, tmin, tmax, tdelta, flat_criterion;
	double low_cutoff = 1.0;
	double probability_move_by_T = 0.0001;
	std::string filename_spin_input;

	in->setData("BIN",  bin, true);
	in->setData("MCSTEP", mcstep, true);

	in->setData("TMIN", tmin, true);
	in->setData("TMAX", tmax, true);
	tdelta = (tmax - tmin) / (double)bin;

	in->setData("LOGFACTOR", logfactor, true);
	in->setData("LOGFLIMIT", logflimit, true);
	in->setData("FLATCRITERION", flat_criterion, true);
	in->setData("PROBABILITYBYT" ,probability_move_by_T, true);

	in->setData("SPININPUT", filename_spin_input);

	in->setData("FLATCHECKSTEP", flatcheck_step);
	in->setData("MINSTEP",   minstep);
	in->setData("LOWCUTOFF", low_cutoff);


	const ParseEcicar ecicar("./ecicar");
	const ParseClusterOut cluster("./cluster.out", ecicar.getIndex());

	WLconfT::WLconfT PoscarSpin("./poscar.spin", in, cluster.getLabel(), cluster.getCluster(), ecicar.getEci(), nullptr, nullptr);

	const int N = PoscarSpin.getSpins().size();

	std::vector<int> index_neglect_bin = PoscarSpin.getNeglectBinIndex();

	std::cout << "----------  initial correlation functions" << std::endl;
	PoscarSpin.dispCorr();
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << "----------  Information about [directZ.ini]" << std::endl;
	dispInput(in);

	if( index_neglect_bin.size() and  filename_spin_input.size() ){
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

	std::string filename_rep_macrostate = "macrostate.out";
	std::ofstream ofs_macrostate(filename_rep_macrostate);
	ofs_macrostate.setf(std::ios_base::fixed, std::ios_base::floatfield);
	ofs_macrostate.close();

	unsigned int fstep = 0;
	unsigned int tstep = 0;
	unsigned int mc_time = 0;
	bool isConverged = false;

	start = std::chrono::system_clock::now();

	auto f_start = std::chrono::system_clock::now();
	auto f_end   = std::chrono::system_clock::now();

 	int final_mcsweep = 0;
 	while(logflimit < logfactor){
		std::vector<double> histogram(bin, 0);
		std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
		std::cout << std::setprecision(10) << "--- fstep " << fstep << " --  log(factor) = " << logfactor << std::endl;
		for(int i=1; ; ++i, ++mc_time){
			for(int j=0; j<N ; ++j){  // スピン数だけステップ回す これで1MCSweep

				int before_index = PoscarSpin.getIndex();

				if( PoscarSpin.RandReal() < probability_move_by_T ) {
					wl_step_by_T(PoscarSpin, dos);
				} else {
					wl_step_by_E(PoscarSpin, dos);
				}

				int index = PoscarSpin.getIndex();

				if( dos[index] == 0 ){
					auto it = find( index_neglect_bin.begin(), index_neglect_bin.end() , index);
					if( it != index_neglect_bin.end() ){
						dos[index] = dos[PoscarSpin.getBeforeIndex()];
						histogram[index] = histogram[PoscarSpin.getBeforeIndex()];
						index_neglect_bin.erase(it);
						PoscarSpin.outputEnergySpin(index, filename_rep_macrostate);
					}
				}

				dos.at(index) += logfactor;
				histogram.at(index) += 1;

			}  /*  end MC sweep */

			if((i % mcstep) == 0){
				WLconfT::outputHistogram(dos,  histogram, tmin, tdelta, fstep);
				std::cout << i << "sweep done " << std::endl;
			}
			if((i % flatcheck_step) == 0 and WLconfT::checkHistogramFlat(histogram, index_neglect_bin, flat_criterion, minstep, low_cutoff)){
				WLconfT::outputHistogram(dos,  histogram, tmin, tdelta, fstep);
				final_mcsweep = i;
				break;
			}
		}
		f_end = std::chrono::system_clock::now();
		auto f_sec = std::chrono::duration_cast<std::chrono::seconds>(f_end-f_start).count();
		std::cout << " fstep = " << fstep << " " << f_sec << " seconds "<< final_mcsweep << " sweeps" << std::endl;
		f_start = std::chrono::system_clock::now();
		logfactor/=2.0;
		++fstep;
	}

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
