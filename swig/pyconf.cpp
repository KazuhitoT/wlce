#include <iomanip>
#include <iterator>
#include <numeric>
#include <random>
#include <chrono>

#include "../src/parser.hpp"
#include "../src/input.hpp"
#include "../src/metroconf.hpp"

std::shared_ptr<Input> in(new Input("pyconf.ini"));
static ParseEcicar ecicar("./ecicar");
static ParseClusterOut cluster("./cluster.out", ecicar.getIndex());
static Metroconf PoscarSpin("./poscar.spin", in, cluster.getLabel(), cluster.getCluster(), ecicar.getEci(), nullptr, nullptr);


constexpr double kb = 0.0000861733; // [eV/K]

double metropolis_step(double temperature){
	PoscarSpin.setTotalEnergy();
	double ene_before = PoscarSpin.getTotalEnergy();

	PoscarSpin.setNewConf();
	double ene_after  = PoscarSpin.getTotalEnergy();

	double b = exp( -(ene_after-ene_before)/kb/temperature );
	if(b>=1.0 or b>PoscarSpin.RandReal()){
	} else {
		PoscarSpin.Memento();
	}
	return PoscarSpin.getTotalEnergy();
}

const std::vector<int> getSpins(){
	return PoscarSpin.getSpins();
}

const double getTotalEnergy(){
	return PoscarSpin.getTotalEnergy();
}

void setEquilibriate(double temperature, int num_steps){
	for(int i=0; i<num_steps; ++i) metropolis_step(temperature);
}

void outputPoscar(){
	PoscarSpin.outputPoscar("latest");
}
