#include "../src/parser.hpp"
#include "../src/conf2corr.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE(conf2corr)

BOOST_AUTO_TEST_CASE(setcorr_2dising_L16)
{
	std::shared_ptr<Input> in(new Input("test.ini"));
	const ParseClusterOut  cluster_out_2dising("./cluster.out.2dising");
	Conf2corr conf2corr_2dising("./poscar.in.2dising", in, cluster_out_2dising.getLabel(), cluster_out_2dising.getCluster());

	for( int i=0; i<100; ++i ) {
		conf2corr_2dising.setMemento();
		conf2corr_2dising.setCorrelationFunction_flip();

		auto corr_2dising  = conf2corr_2dising.getCorrelationFunctions();
		conf2corr_2dising.setInitialCorrelationFunction();
		auto corr_2dising_confirm = conf2corr_2dising.getCorrelationFunctions();

		for( int j=0; j<corr_2dising.size(); ++j){
			for( int k=0; k<corr_2dising[j].size(); ++k){
				BOOST_CHECK_CLOSE(corr_2dising[j][k], corr_2dising_confirm[j][k], 0.00001);
			}
		}
	}

	for( int i=0; i<100; ++i ) {
		conf2corr_2dising.setMemento();
		conf2corr_2dising.setCorrelationFunction_flip();

		auto corr_2dising_exchange  = conf2corr_2dising.getCorrelationFunctions();
		conf2corr_2dising.setInitialCorrelationFunction();
		auto corr_2dising_exchange_confirm = conf2corr_2dising.getCorrelationFunctions();

		for( int j=0; j<corr_2dising_exchange.size(); ++j){
			for( int k=0; k<corr_2dising_exchange[j].size(); ++k){
				BOOST_CHECK_CLOSE(corr_2dising_exchange[j][k], corr_2dising_exchange_confirm[j][k], 0.00001);
			}
		}
	}

}


BOOST_AUTO_TEST_CASE(setcorr_fcc_unit)
{
	std::shared_ptr<Input> in(new Input("test.ini"));
	const ParseClusterOut  cluster_out_fcc_unit("./cluster.out.fcc.unit");
	Conf2corr conf2corr_fcc_unit("./poscar.in.fcc.unit", in, cluster_out_fcc_unit.getLabel(), cluster_out_fcc_unit.getCluster());

	for( int i=0; i<100; ++i ) {
		conf2corr_fcc_unit.setMemento();
		conf2corr_fcc_unit.setCorrelationFunction_flip();

		auto corr_fcc_unit = conf2corr_fcc_unit.getCorrelationFunctions();
		conf2corr_fcc_unit.setInitialCorrelationFunction();
		auto corr_fcc_unit_confirm = conf2corr_fcc_unit.getCorrelationFunctions();

		for( int j=0; j<corr_fcc_unit.size(); ++j){
			for( int k=0; k<corr_fcc_unit[j].size(); ++k){
				BOOST_CHECK_SMALL( (corr_fcc_unit_confirm[j][k]-corr_fcc_unit[j][k]), 0.00001);
			}
		}

	}

	for( int i=0; i<100; ++i ){
		 conf2corr_fcc_unit.setMemento();
		 conf2corr_fcc_unit.setCorrelationFunction_exchange();

		 auto corr_fcc_unit_exchange  =  conf2corr_fcc_unit.getCorrelationFunctions();
		 conf2corr_fcc_unit.setInitialCorrelationFunction();
		 auto corr_fcc_unit_exchange_confirm =  conf2corr_fcc_unit.getCorrelationFunctions();

		 for( int j=0; j<corr_fcc_unit_exchange.size(); ++j){
			 for( int k=0; k<corr_fcc_unit_exchange[j].size(); ++k){
				 BOOST_CHECK_SMALL( (corr_fcc_unit_exchange_confirm[j][k]-corr_fcc_unit_exchange[j][k]), 0.00001);
			 }
		 }

	}

}

BOOST_AUTO_TEST_SUITE_END()
