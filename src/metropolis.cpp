#include <iomanip>
#include <iterator>
#include <numeric>
#include <random>
#include <chrono>

#include "./parser.hpp"
#include "./input.hpp"
#include "./metroconf.hpp"


constexpr double kb = 0.0000861733; // [eV/K]

double metropolis_step(Metroconf& conf, double temperature){
	conf.setTotalEnergy();
	double ene_before = conf.getTotalEnergy();

	conf.setNewConf();
	double ene_after  = conf.getTotalEnergy();

	double b = exp( -(ene_after-ene_before)/kb/temperature );
	if(b>=1.0 or b>conf.RandReal()){
	} else {
		conf.Memento();
	}
	return conf.getTotalEnergy();
}

int main(int argc, char* argv[]){
	std::shared_ptr<Input> in(new Input("metropolis.ini"));

	int mcstep, samplestep;
	std::vector<double> temperature_input;

	in->setData("MCSTEP",      mcstep, true);
	in->setData("SAMPLESTEP",  samplestep, true);
	in->setData("TEMPERATURE", temperature_input, true);

	std::vector<double> vec_temperature;
	if( temperature_input.size() == 1 ){
		vec_temperature = temperature_input;
	} else if ( temperature_input.size() == 3 ) {
		if( temperature_input[0] > temperature_input[2] ){ // 900 100 100
			for( double iniT=temperature_input[0]; iniT>temperature_input[2]; iniT-=temperature_input[1]){
				vec_temperature.push_back(iniT);
			}
		} else if( temperature_input[0] < temperature_input[2] ){  // 100 100 900
			for( double iniT=temperature_input[0]; iniT<temperature_input[2]; iniT+=temperature_input[1]){
				vec_temperature.push_back(iniT);
			}
		} else {
			vec_temperature.push_back(temperature_input[0]);
		}
	} else {
		std::cerr << " ERROR : TEMPERATURE in metropolis.ini must be 1 or 3 arguments." << std::endl;
		exit(1);
	}

	const ParseLabels label("./labels.in");
	const ParseEcicar ecicar("./ecicar");
	const ParseMultiplicityIn multiplicity_in("./multiplicity.in", ecicar.getIndex());
	const ParseClusterIn  cluster_in("./clusters.in", ecicar.getIndex(), multiplicity_in.getMultiplicityIn());

	Metroconf PoscarSpin("./poscar.spin", in, label.getLabels(), cluster_in.getCluster(), ecicar.getEci(), nullptr, nullptr);

	const int N = PoscarSpin.getSpins().size();

	std::cout << "----------  initial correlation functions" << std::endl;
	PoscarSpin.dispCorr();
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << "----------  Information about [metropolis.ini]" << std::endl;
	PoscarSpin.dispInput();
	std::cout << std::endl;
	std::cout << std::endl;

	auto start = std::chrono::system_clock::now();

	std::cout << "----------  Log  ----------" << std::endl;
	std::cout << std::endl;

	std::vector<double> vec_ave;
	std::vector<double> vec_var;

	std::ofstream ofs("metropolis.out");
	ofs.setf(std::ios_base::fixed, std::ios_base::floatfield);

	int loop_count = 0;
	for( const auto& temperature : vec_temperature ){
		std::vector<double> vec_ene;

		for(int i=1; i<=mcstep; ++i){
			for(int j=0; j<N; ++j){
				metropolis_step(PoscarSpin, temperature);
			}
			std::cout << temperature << "[K] : " << i << " step done." << std::endl;
		}

		PoscarSpin.outputPoscar("T="+std::to_string((int)temperature));

		for(int i=1; i<=mcstep; ++i){
			for(int j=0; j<N; ++j){
				double e = metropolis_step(PoscarSpin, temperature);
				vec_ene.push_back(e);
			}
			std::cout << "SAMPLE " << temperature << "[K] : " << i << " step done." << std::endl;
		}

		double average  = std::accumulate(vec_ene.begin(), vec_ene.end(), 0.0 ) / double(vec_ene.size());
		double variance = std::accumulate(vec_ene.begin(), vec_ene.end(), 0.0,
			[average]( double sum, double val ){ return sum+(val-average)*(val-average);} ) / double(vec_ene.size());

		vec_ave.push_back(average);
		vec_var.push_back(variance);

		ofs << vec_temperature[loop_count] << " " << vec_ave[loop_count] << " " << vec_var[loop_count] << std::endl;
		++loop_count;
	}

	ofs.close();


	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "check final correlation functions" << std::endl;
	PoscarSpin.dispCorr();
	std::cout << "---------------------------" << std::endl;
	PoscarSpin.setInitialCorrelationFunction();
	PoscarSpin.dispCorr();


	return 0;
}
