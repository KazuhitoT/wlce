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

#include "./parser.hpp"
#include "./input.hpp"

using allclusters = std::vector<std::vector<std::vector<std::vector<int>>>>;
using indexorders = std::vector<std::vector<std::vector<std::vector<int>>>>;
using basisfunc   = std::vector<std::vector<double>>;

class Conf2corr {
	private:
		std::vector<double> compositions;
		std::vector<double> compositions_before;
		std::vector<double> spins;
		std::vector<double> spins_before;
		std::vector<std::vector<double>> correlation_functions;
		std::vector<std::vector<double>> correlation_functions_before;

		/* vec[basis][degree] vec[basis]->polynomials */
		std::shared_ptr<std::vector<std::vector<double>>> pbasis_functions;
		std::shared_ptr<basisfunc> pbasis_clusters;
		/* vec[num_in_cluster][index][combinations][order] */
		std::shared_ptr<indexorders> pindex_orders;

		/* --------- Input parameter for Conf2corr --------- */
		std::vector<double> spinposcar;
		std::vector<double> spince;
		/*  allclusters index->site->multiplicity->clusters */
		std::shared_ptr<allclusters> pall_clusters;
		/* ------------------------------------------------- */

		std::random_device rd;
		std::mt19937 mt;
		std::function<int ()> rnd_int_N;
		std::function<double ()> rnd_real;

	public:
		Conf2corr(){};
		Conf2corr(char* filename,
			std::vector<double> _spinposcar,
			std::vector<double> _spince,
			std::shared_ptr<labels> _plabels,
			std::shared_ptr<allclusters> _pall_clusters,
			std::shared_ptr<basisfunc>   _pbasis_functions = nullptr,
			std::shared_ptr<indexorders> _pindex_orders = nullptr
		);
		Conf2corr(char* filename,
			std::shared_ptr<Input> _in,
			std::shared_ptr<labels> _plabels,
			std::shared_ptr<allclusters> _pall_clusters,
			std::shared_ptr<basisfunc>   _pbasis_functions = nullptr,
			std::shared_ptr<indexorders> _pindex_orders = nullptr
		);
		~Conf2corr(){};

		void setSpins(char* filename, std::shared_ptr<labels> plabels);
		void setSpins(std::vector<double> _spins){ this->spins = _spins; }
		void setSpins(int i, double spin){ this->spins.at(i) = spin; }
		void setSpinsBefore(std::vector<double> _spins){ this->spins_before = _spins;}
		void setSpinsRandom(){
			std::random_shuffle(this->spins.begin(), this->spins.end());
			this->spins_before = this->spins;
			this->setInitialCorrelationFunction();
		}
		void setCompositions();

		void setBasisCoefficient();
		void setIndexOrders();
		double getBasisFunction(int, int);

		void setInitialCorrelationFunction();
		void setCorrelationFunction_flip();
		void setCorrelationFunction_exchange();

		virtual void setMemento()
		{
			this->spins_before = this->spins;
			this->correlation_functions_before = this->correlation_functions;
			this->compositions_before = this->compositions;
		};

		virtual void Memento()
		{
			this->spins = this->spins_before;
			this->correlation_functions = this->correlation_functions_before;
			this->compositions = this->compositions_before;
		};

		void dispCorr();

		int RandN(){return rnd_int_N();};
		double RandReal(){ return rnd_real();};

		std::vector<double> getSpins(){ return spins; };
		std::vector<double> getCompositions(){ return compositions; };
		double getCompositions(int i){ return compositions[i]; };
		std::vector<double> getSpinCE(){ return spince; };
		std::vector<std::vector<double>> getCorrelationFunctions(){ return correlation_functions; };
		double getCorrelationFunctions(int i, int j){ return correlation_functions[i][j]; };


  	// Conf2corr &operator=(const Conf2corr&);    // 代入演算子
};

#endif
