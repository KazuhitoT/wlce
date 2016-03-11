#ifndef CONF2CORR_HPP
#define CONF2CORR_HPP

#include <string>
#include <map>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <typeinfo>
#include <fstream>
#include <iostream>
#include <functional>
#include <boost/algorithm/string.hpp>
#include <Eigen/Core>
#include <Eigen/LU>
#include <memory>
#include <cfloat>
#include "./parsePoscar.hpp"

using allclusters = std::vector<std::vector<std::vector<std::vector<int>>>>;
using indexorders = std::vector<std::vector<std::vector<std::vector<int>>>>;
using basisfunc   = std::vector<std::vector<double>>;

class Conf2corr {
	private:
		friend class Conf2corrTest;

		std::vector<double> spins;
		std::vector<double> spins_before;
		std::vector<std::vector<double>> correlation_functions;
		std::vector<std::vector<double>> correlation_functions_before;

		/* vec[basis][degree] vec[basis]->polynomials */
		// std::vector<std::vector<double>> basis_functions;
		std::shared_ptr<std::vector<std::vector<double>>> pbasis_functions;
		std::shared_ptr<basisfunc> pbasis_clusters;
		/* vec[num_in_cluster][index][combinations][order] */
		// std::vector<std::vector<std::vector<std::vector<int>>>> index_orders;
		std::shared_ptr<indexorders> pindex_orders;

		/* --------- Input parameter for Conf2corr --------- */
		std::vector<double> spinposcar;
		std::vector<double> spince;
		/*  allclusters index->site->multiplicity->clusters */
		std::shared_ptr<allclusters> pall_clusters;
		/* ------------------------------------------------- */

	public:
		Conf2corr(){};
		// Conf2corr(char* filename,
		// 	std::vector<double> _spinposcar,
		// 	std::vector<double> _spince,
		// 	std::shared_ptr<allclusters> _pall_clusters,
		// 	std::shared_ptr<basisfunc>   _pbasis_functions = nullptr,
		// 	std::shared_ptr<indexorders> _pindex_orders = nullptr
		// );
		Conf2corr(std::vector<double> _spins,
			std::vector<double> _spinposcar,
			std::vector<double> _spince,
			std::shared_ptr<allclusters> _pall_clusters,
			std::shared_ptr<basisfunc>   _pbasis_functions = nullptr,
			std::shared_ptr<indexorders> _pindex_orders = nullptr
		);
		~Conf2corr(){};

		void setSpins(char* filename);

		void setBasisCoefficient();
		void setIndexOrders();

		double getBasisFunction(int, int);
		void setInitialCorrelationFunction();

		//
		//
		// std::vector<double>             getSpins()                const {return spins;};
		// std::vector<double>             getBeforeSpins()                const {return spins_before;};
		// std::vector<std::vector<double>>             getCorrelationFunctions() const {return correlation_functions;};
		// std::vector<std::vector<double>>             getBeforeCorrelationFunctions() const {return correlation_functions_before;};
		// double getBasisFunction(/* order */  int, /* spince_num */ int);
		// std::vector<std::vector<double>> getBasisFunctions() const {return basis_functions;}
		// std::vector<std::vector<std::vector<std::vector<int>>>> getIndexOrders() const {return index_orders;}
		// // double getNumInCluster(int line_num_dimcar) const { return dimcar.getDimcar(line_num_dimcar).first;}
		// // double getNumOfCluster(int line_num_dimcar) const { return dimcar.getDimcar(line_num_dimcar).second;}
		//
		void dispCorr();
		// void dispSpin();
		// void dispIndexOrder();
		// void dispBasisCoefficient();
		//
		// void outputCorrSpin(int num=-1, int mode=0, std::string filename="corrspin.out");
		// void outputCorr(std::string filename = "corr-profile.out");
		// void outputPoscar(std::string filename = "poscar.out");

  	Conf2corr &operator=(const Conf2corr&);    // 代入演算子
};

#endif
