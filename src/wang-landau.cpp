#include <iomanip>
#include <iterator>
#include <numeric>
#include <random>
#include <chrono>

#include "./parser.hpp"
#include "./input.hpp"
#include "./wlconf.hpp"

bool checkHistogramFlat(const std::vector<double>& histogram, const std::vector<int>& index_neglect_bin, double lflat, int minstep, double low_cutoff=1){
	double ave = accumulate(histogram.begin(), histogram.end(), 0) / (double)(histogram.size()-index_neglect_bin.size());
	double limit = ave * lflat ;
	for(int i=0, imax=histogram.size(); i<imax; ++i){

		auto it = find( index_neglect_bin.begin(), index_neglect_bin.end() , i);
		if( it != index_neglect_bin.end() ) continue;

		if( minstep>0 and histogram.at(i)<minstep ) return false;

		if(histogram.at(i) < (ave*(1.0-low_cutoff))) continue;
		else if(histogram.at(i) < limit) return false;
	}
	return true;
}

void outputHistogram(const std::vector<double>& dos, const std::vector<double>& histogram, double emin, double delta, int fstep, std::string prefix = ""){
		std::ofstream ofs("out-wl"+prefix+std::to_string(fstep)+".dat");
		ofs.setf(std::ios_base::fixed, std::ios_base::floatfield);
		for(int i=0, imax=dos.size(); i<imax; ++i){
			ofs << i << " " << std::setprecision(10) << emin + i*delta << " " << emin + (i+1)*delta << " " <<
			std::setprecision(10) << dos[i] << " " << histogram[i] << std::endl;
		}
		ofs.close();
}

double wl_step(WLconf& conf, const std::vector<double>& dos){
	conf.setTotalEnergy();
	conf.setIndex();
	int index_before = conf.getIndex();
	int index_after  = index_before;

	conf.setNewConf();
	index_after  = conf.getIndex();

	double b = exp(dos.at(index_before) - dos.at(index_after));
	if(b>=1.0 or b>conf.RandReal()){
	} else {
		conf.Memento();
	}
	return conf.getTotalEnergy();
}

int main(int argc, char* argv[]){
	std::shared_ptr<Input> in(new Input("wang-landau.ini"));

	int mcstep, bin, flatcheck_step, minstep=0;
	double logfactor, logflimit, emin, emax, edelta, flat_criterion;
	double low_cutoff = 1.0;
	std::string filename_spin_input;

	in->setData("BIN",  bin, true);
	in->setData("MCSTEP", mcstep, true);

	in->setData("EMIN", emin, true);
	in->setData("EMAX", emax, true);
	edelta = (emax - emin) / (double)bin;

	in->setData("FLATCHECKSTEP", flatcheck_step, true);
	in->setData("LOGFACTOR", logfactor, true);
	in->setData("LOGFLIMIT", logflimit, true);
	in->setData("FLATCRITERION", flat_criterion, true);

	in->setData("SPININPUT", filename_spin_input);

	in->setData("MINSTEP",   minstep);
	in->setData("LOWCUTOFF", low_cutoff);


	const ParseLabels label("./labels.in");
	const ParseEcicar ecicar("./ecicar");
	const ParseMultiplicityIn multiplicity_in("./multiplicity.in", ecicar.getIndex());
	const ParseClusterIn  cluster_in("./clusters.in", ecicar.getIndex(), multiplicity_in.getMultiplicityIn());

	WLconf PoscarSpin("./poscar.spin", in, label.getLabels(), cluster_in.getCluster(), ecicar.getEci(), nullptr, nullptr);

	const int N = PoscarSpin.getSpins().size();

	std::vector<int> index_neglect_bin = PoscarSpin.getNeglectBinIndex();

	std::cout << "----------  initial correlation functions" << std::endl;
	PoscarSpin.dispCorr();
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << "----------  Information about [wang-landau.ini]" << std::endl;
	PoscarSpin.dispInput();

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

				wl_step(PoscarSpin, dos);
				PoscarSpin.setIndex();
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
				outputHistogram(dos,  histogram, emin, edelta, fstep);
				std::cout << i << "sweep done " << std::endl;
			}
			if((i % flatcheck_step) == 0 and checkHistogramFlat(histogram, index_neglect_bin, flat_criterion, minstep, low_cutoff)){
				outputHistogram(dos,  histogram, emin, edelta, fstep);
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
