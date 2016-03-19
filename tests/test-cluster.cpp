// #include <memory>
// #include <iterator>
// #include <Eigen/Core>
// #include <Eigen/LU>
#include "../src/parser.hpp"
#include "../src/site.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// class Conf2corrTest : public ::Conf2corr {
// public:
// 	Conf2corrTest(char* filename, vector<double> spinposcar, vector<double> spince, Ensemble ensemble = Ensemble::ALLOY) : Conf2corr(filename, spince, spinposcar, ensemble) {};
// // protected:
// 	// T1 get_private_member1(Conf2corr* obj) {
// 	vector<vector<vector<vector<int>>>> getIndexOrders() {
// 		return Conf2corr::index_orders;
// 		// return obj->private_member1_;
// 	}
// };

BOOST_AUTO_TEST_SUITE(cluster)

BOOST_AUTO_TEST_CASE(multiplicity)
{

	// const ParseEcicar ecicar("./ecicar");
	const ParseMultiplicityIn multiplicity_in_2dising("./multiplicity.in.2dising");
	const ParseClusterIn  cluster_in_2dising("./clusters.in", multiplicity_in_2dising.getMultiplicityIn());

	std::map<int, std::pair<int,int>> m2dising;
	for( const auto& i : multiplicity_in_2dising.getMultiplicityIn() ){
		/*  i.first == numpoints_in_cluster  i.second == numclusters_for_a_site */

	}
	// Conf2corrTest obj("./POSCAR.spin_s3", vector<double> {-1, 0, 1}, vector<double> {-1, 0, 1}, Ensemble::SPIN);
	// auto index_orders = obj.getIndexOrders();
	// cout << "order list" << endl;
	// for(const auto& i : index_orders){
	// 	for(const auto& j : i){
	// 		for(const auto& k : j){
	// 			for(const auto& l : k){
	// 				cout << l;
	// 			}
	// 			cout << " ";
	// 		}
	// 	}
	// 	cout << endl;
	// }
	// cout << "checking random walk..." << endl;
	// for(int i=0; i<1000; ++i){
	// 	obj.setCorrelationFunction();
	// }
	// auto corr  = obj.getCorrelationFunctions();
	// obj.setCorrelationFunctionFromClucar();
	// auto corr_clucar = obj.getCorrelationFunctions();
	//
	// for(int i=0; i<corr.size(); ++i){
	// 	for(int j=0; j<corr[i].size(); ++j){
	// 		BOOST_CHECK_CLOSE(corr[i][j], corr_clucar[i][j], 0.00001);
	// 	}
	// }

}

BOOST_AUTO_TEST_SUITE_END()
