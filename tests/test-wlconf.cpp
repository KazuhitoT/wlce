// #include "../src/parser.hpp"
// #include "../src/site.hpp"
// #include "../src/wlconf.hpp"
//
// #include <boost/test/unit_test.hpp>
// #include <boost/test/floating_point_comparison.hpp>
//
// BOOST_AUTO_TEST_SUITE(wlconf)
//
// BOOST_AUTO_TEST_CASE(setcorr)
// {
// 	const ParseLabels labels_in("./labels.in.2dising");
// 	const ParseEcicar ecicar("./ecicar.2dising");
// 	const ParseMultiplicityIn multiplicity_in_2dising("./multiplicity.in.2dising", ecicar.getIndex());
// 	const ParseClusterIn  cluster_in_2dising("./clusters.in.2dising", ecicar.getIndex(), multiplicity_in_2dising.getMultiplicityIn());
//
// 	WLconf wl2dising("poscar.in.2dising", std::vector<double> {-1, 1}, std::vector<double> {-1, 1}, labels_in.getLabels(), cluster_in_2dising.getCluster(), ecicar.getEci()); //, nullptr, nullptr);
//
// 	for( int i=0; i<1000; ++i ){
// 		wl2dising.setCorrelationFunction();
// 	}
// 	auto corr_2dising  = wl2dising.getCorrelationFunctions();
// 	// wl2dising.dispCorr();
// 	wl2dising.setInitialCorrelationFunction();
// 	auto _corr_2dising = wl2dising.getCorrelationFunctions();
// 	// wl2dising.dispCorr();
//
// 	for( int i=0; i<corr_2dising.size(); ++i){
// 		for( int j=0; j<corr_2dising[i].size(); ++j){
// 			BOOST_CHECK_CLOSE(corr_2dising[i][j], _corr_2dising[i][j], 0.00001);
// 		}
// 	}
//
// 	WLconf wl2dising_exchange("poscar.in.2dising", std::vector<double> {-1, 1}, std::vector<double> {-1, 1}, labels_in.getLabels(), cluster_in_2dising.getCluster(), ecicar.getEci(), nullptr, nullptr, std::vector<double> {0, 0});
//
// 	// wl2dising_exchange.dispCorr();
// 	for( int i=0; i<10000; ++i ){
// 		wl2dising_exchange.setCorrelationFunction();
// 	}
// 	auto corr_2dising_exchange  = wl2dising_exchange.getCorrelationFunctions();
// 	// wl2dising_exchange.dispCorr();
// 	wl2dising_exchange.setInitialCorrelationFunction();
// 	auto _corr_2dising_exchange = wl2dising_exchange.getCorrelationFunctions();
// 	// wl2dising_exchange.dispCorr();
//
// 	for( int i=0; i<corr_2dising_exchange.size(); ++i){
// 		for( int j=0; j<corr_2dising_exchange[i].size(); ++j){
// 			BOOST_CHECK_CLOSE(corr_2dising_exchange[i][j], _corr_2dising_exchange[i][j], 0.00001);
// 		}
// 	}
//
// }
//
// BOOST_AUTO_TEST_SUITE_END()
