#include <iomanip>
#include <iterator>
#include <numeric>
#include <random>
#include <chrono>

#include "./parser.hpp"
#include "./input.hpp"
#include "./wlconf.hpp"

int main(int argc, char* argv[]){
	std::shared_ptr<Input> in(new Input("totalene.ini"));

	const ParseLabels labels_in("./labels.in");
	const ParseEcicar ecicar("./ecicar");
	const ParseMultiplicityIn multiplicity_in("./multiplicity.in", ecicar.getIndex());
	const ParseClusterIn  cluster_in("./clusters.in", ecicar.getIndex(), multiplicity_in.getMultiplicityIn());

	WLconf PoscarSpin("./poscar.expand.spin", in, labels_in.getLabels(), cluster_in.getCluster(), ecicar.getEci(), nullptr, nullptr);
	PoscarSpin.dispCorr();
	PoscarSpin.setTotalEnergy();
	std::cout << std::endl;
	std::cout << PoscarSpin.getTotalEnergy() << std::endl;

}
