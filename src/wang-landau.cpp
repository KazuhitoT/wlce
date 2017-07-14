#include <iomanip>
#include <iterator>
#include <numeric>
#include <random>
#include <chrono>

#include "./parser.hpp"
#include "./input.hpp"
#include "./wlconf.hpp"

void dispInput(std::shared_ptr<Input> in){
	//  [INPUT]
	std::cout << "EMIN          : " << in->getDataByString("EMIN") << std::endl;
	std::cout << "EMAX          : " << in->getDataByString("EMAX") << std::endl;
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

double wl_step(WLconf::WLconf& conf, const std::vector<double>& dos){
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

	int mcstep, bin, flatcheck_step, num_ignore_edge_index=0, minstep=0;
	double logfactor, logflimit, emin, emax, edelta, flat_criterion;
	double low_cutoff = 1.0;
	int is_restart=0, is_output_tunnelingtime=0;

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
	in->setData("NUMIGNOREEDGEINDEX",   minstep);
	in->setData("LOWCUTOFF", low_cutoff);
	in->setData("OUTTUNNELINGTIME", is_output_tunnelingtime);

	in->setData("RESTART", is_restart);


	const ParseEcicar ecicar("./ecicar");
	const ParseClusterOut cluster("./cluster.out", ecicar.getIndex());

	WLconf::WLconf PoscarSpin("./poscar.spin", in, cluster.getLabel(), cluster.getCluster(), ecicar.getEci());

	const int N = PoscarSpin.getSpins().size();

	std::vector<int> index_neglect_bin = PoscarSpin.getNeglectBinIndex();

	std::cout << "----------  initial correlation functions" << std::endl;
	PoscarSpin.dispCorr();
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << "----------  Information about [wang-landau.ini]" << std::endl;
	in->disp();

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
	unsigned int tunneling_time_start = 0;
	unsigned int tunneling_time_end   = 0;

	bool isConverged = false;

	start = std::chrono::system_clock::now();

	auto f_start = std::chrono::system_clock::now();
	auto f_end   = std::chrono::system_clock::now();

	/*  set lowest/highest index for tunneling-time calculation  */
	int index_lowest = bin-1, index_highest = 0;
	if( is_restart>0 ){
		WLconf::restart(fstep, logfactor, dos, index_neglect_bin);

		if( index_neglect_bin.size()>0 ){

			for(int i=0; i<bin; ++i){
				if( i != index_neglect_bin[i] ){
					index_lowest = i;
					break;
				}
			}

			for(int i=bin; i>=0; --i){
				if( i != index_neglect_bin[i] ){
					index_highest = i;
					break;
				}
			}

		}

	}

	enum class TunnelingTimeStartPoint
	{
		Lowest,
		Highest,
	};
	TunnelingTimeStartPoint randomwalk_start_point = TunnelingTimeStartPoint::Lowest;
	std::ofstream ofs_tunneling_time("tunneling_time.out");
	ofs_tunneling_time.setf(std::ios_base::fixed, std::ios_base::floatfield);


 	int final_mcsweep = 0;
 	while(logflimit < logfactor){
		std::vector<double> histogram(bin, 0);
		std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
		std::cout << std::setprecision(10) << "--- fstep " << fstep << " --  log(factor) = " << logfactor << std::endl;
		for(int i=1; ; ++i, ++mc_time){
			for(int j=0; j<N ; ++j){  // start a MC sweep

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

						if( index < index_lowest  ) index_lowest  = index;
						if( index > index_highest ) index_highest = index;

					}
				}

				if( index <= (index_lowest+num_ignore_edge_index) and randomwalk_start_point == TunnelingTimeStartPoint::Highest ){
					tunneling_time_end   = mc_time;
					ofs_tunneling_time << fstep << " " << index << " " << tunneling_time_end - tunneling_time_start << std::endl;
					tunneling_time_start = tunneling_time_end;
					randomwalk_start_point =  TunnelingTimeStartPoint::Lowest;
				} else if ( index >= (index_highest-num_ignore_edge_index) and randomwalk_start_point == TunnelingTimeStartPoint::Lowest ){
					tunneling_time_end   = mc_time;
					ofs_tunneling_time << fstep << " " << index << " " << tunneling_time_end - tunneling_time_start << std::endl;
					tunneling_time_start = tunneling_time_end;
					randomwalk_start_point =  TunnelingTimeStartPoint::Highest;
				}

				dos[index] += logfactor;
				histogram[index] += 1;

			}  /*  end a MC sweep */

			if((i % mcstep) == 0){
				WLconf::outputHistogram(dos,  histogram, emin, edelta, fstep);
				PoscarSpin.outputPoscar("latest");
				std::cout << i << "sweep done " << std::endl;
			}
			if((i % flatcheck_step) == 0 and WLconf::checkHistogramFlat(histogram, index_neglect_bin, flat_criterion, minstep, low_cutoff)){
				WLconf::outputHistogram(dos,  histogram, emin, edelta, fstep);
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

	ofs_tunneling_time.close();
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
