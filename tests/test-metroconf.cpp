#include "../src/parser.hpp"
#include "../src/metroconf.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE(metroconf)

BOOST_AUTO_TEST_CASE(setNewConf)
{
	std::shared_ptr<Input> in(new Input("test.ini"));
	const ParseEcicar ecicar("./ecicar.2dising");
	const ParseClusterOut cluster("./cluster.out.2dising", ecicar.getIndex());

	Metroconf PoscarSpin("./poscar.in.2dising", in, cluster.getLabel(), cluster.getCluster(), ecicar.getEci(), nullptr, nullptr);

	for( int i=0; i<100; ++i ) {
		PoscarSpin.setNewConf();

		auto corr_2dising  = PoscarSpin.getCorrelationFunctions();
		PoscarSpin.setInitialCorrelationFunction();
		auto corr_2dising_confirm = PoscarSpin.getCorrelationFunctions();

		for( int j=0; j<corr_2dising.size(); ++j){
			for( int k=0; k<corr_2dising[j].size(); ++k){
				BOOST_CHECK_CLOSE(corr_2dising[j][k], corr_2dising_confirm[j][k], 0.00001);
			}
		}
	}

	for( int i=0; i<100; ++i ) {
		PoscarSpin.setNewConf();
		PoscarSpin.Memento();

		auto corr_2dising  = PoscarSpin.getCorrelationFunctions();
		PoscarSpin.setInitialCorrelationFunction();
		auto corr_2dising_confirm = PoscarSpin.getCorrelationFunctions();

		for( int j=0; j<corr_2dising.size(); ++j){
			for( int k=0; k<corr_2dising[j].size(); ++k){
				BOOST_CHECK_CLOSE(corr_2dising[j][k], corr_2dising_confirm[j][k], 0.00001);
			}
		}
	}


}

BOOST_AUTO_TEST_SUITE_END()
