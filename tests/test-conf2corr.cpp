#include "../src/parser.hpp"
#include "../src/conf2corr.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE(conf2corr)

BOOST_AUTO_TEST_CASE(setcorr)
{
	std::shared_ptr<Input> in(new Input("test.ini"));

	const ParseLabels labels_in("./labels.in.2dising");
	const ParseEcicar ecicar("./ecicar.2dising");
	const ParseMultiplicityIn multiplicity_in_2dising("./multiplicity.in.2dising", ecicar.getIndex());
	const ParseClusterIn  cluster_in_2dising("./clusters.in.2dising", ecicar.getIndex(), multiplicity_in_2dising.getMultiplicityIn());

	Conf2corr conf2corr_2dising("./poscar.in.2dising", in, labels_in.getLabels(), cluster_in_2dising.getCluster());

	for( int i=0; i<1000; ++i ) {
		conf2corr_2dising.setMemento();
		conf2corr_2dising.setCorrelationFunction_flip();
	}

	auto corr_2dising  = conf2corr_2dising.getCorrelationFunctions();
	conf2corr_2dising.setInitialCorrelationFunction();
	auto corr_2dising_confirm = conf2corr_2dising.getCorrelationFunctions();

	for( int i=0; i<corr_2dising.size(); ++i){
		for( int j=0; j<corr_2dising[i].size(); ++j){
			BOOST_CHECK_CLOSE(corr_2dising[i][j], corr_2dising_confirm[i][j], 0.00001);
		}
	}

	for( int i=0; i<1000; ++i ){
		conf2corr_2dising.setMemento();
		conf2corr_2dising.setCorrelationFunction_exchange();
	}

	auto corr_2dising_exchange  = conf2corr_2dising.getCorrelationFunctions();
	conf2corr_2dising.setInitialCorrelationFunction();
	auto corr_2dising_exchange_confirm = conf2corr_2dising.getCorrelationFunctions();

	for( int i=0; i<corr_2dising_exchange.size(); ++i){
		for( int j=0; j<corr_2dising_exchange[i].size(); ++j){
			BOOST_CHECK_CLOSE(corr_2dising_exchange[i][j], corr_2dising_exchange_confirm[i][j], 0.00001);
		}
	}

}

BOOST_AUTO_TEST_SUITE_END()
