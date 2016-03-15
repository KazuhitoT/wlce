
#include <iomanip>
#include <iterator>
#include <numeric>
#include <random>
#include <chrono>

#include "./Parser.hpp"
#include "./Input.hpp"
#include "./WLconf.hpp"

int main(){
	const ParseMultiplicityIn multiplicity_in("./multiplicity.in");
	const ParseClusterIn  cluster_in("./clusters.in", multiplicity_in.getMultiplicityIn());

	Input in("corrdump.ini");
	std::vector<double> spince;
	in.setData("SPINCE", spince, true);

	Conf2corr PoscarSpin("./poscar.spin", spince, spince, cluster_in.getCluster());
	PoscarSpin.dispCorr();
}
