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
#include <unordered_map>
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
		/* spins consists of index of spince */
		std::vector<int> spins;
		std::vector<std::vector<double>> correlation_functions;
		/*  tuple[spin_index, spin_before, spin_after] */
		std::vector<std::tuple<int, int, int>> vec_changed_spins;
		/* all_calculated_basis_functions[order][spin -> basisfunc] */
		std::vector<std::vector<double>> all_calculated_basis_functions;

		/* vec[basis][degree] vec[basis]->polynomials */
		std::shared_ptr<basisfunc> pbasis_functions;
		/*
		 * vec[num_in_cluster][index][combinations][order]
		 * 1 <= order <= compositions.size()-1
		*/
		std::shared_ptr<indexorders> pindex_orders;
		/*  allclusters index->site->multiplicity->clusters */
		std::shared_ptr<allclusters> pall_clusters;

		const std::shared_ptr<ParsePoscar> pposcar_spin;

		/* --------- Input parameter for Conf2corr --------- */
		std::vector<double> spinposcar;
		std::vector<double> spince;  // if 24 8 in poscar at line t, [0]-> -1 for 24 atoms and [1] -> 1 for 8 atoms
		std::vector<int> index_for_check_inside;
		double spin_for_check_inside = -10000;
		/* ------------------------------------------------- */

		std::mt19937 mt;
		std::function<int ()> rnd_int_N;
		std::function<int ()> rnd_int_spince_index;
		std::function<double ()> rnd_real;

	public:
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

		void setSpins(std::shared_ptr<labels> plabels);
		void setSpins(std::vector<int> _spins){ this->spins = _spins; }
		void setSpins(int i, double spin){ this->spins.at(i) = spin; }
		void setSpinsRandom(){
			std::shuffle(this->spins.begin(), this->spins.end(), mt);
			this->setInitialCorrelationFunction();
		}
		void setCompositions();

		void setBasisCoefficient();
		void setIndexOrders();
		void setBasisFunction();

		void setInitialCorrelationFunction();
		void setCorrelationFunction_flip(int lattice_point = -1, int after_spin = -1);
		void setCorrelationFunction_exchange();

		double calcCorrelationFunctionNorm(double p);

		virtual void setMemento()
		{
		};

		virtual void Memento()
		{
			for(const auto& changed_spin : this->vec_changed_spins){
				this->setCorrelationFunction_flip(std::get<0>(changed_spin), std::get<1>(changed_spin));
			}
		};

		bool isInNthNearestNeighborPair(int lattice_point);

		void dispCorr();
		void outputPoscar(std::string prefix="out");

		int RandN(){return rnd_int_N();};
		double RandReal(){ return rnd_real();};

		std::vector<int> const& getSpins() const { return spins; };
		std::vector<double> const& getCompositions() const{ return compositions; };
		double const& getCompositions(int i) const{ return compositions[i]; };
		std::vector<double> const& getSpinCE() const{ return spince; };
		std::vector<double> const& getSpinPoscar() const{ return spinposcar; };
		std::vector<std::vector<double>> const& getCorrelationFunctions() const { return correlation_functions; };
		double const& getCorrelationFunctions(int i, int j) const{ return correlation_functions[i][j]; };

		double getBasisFunction(int order, int _spin) {return all_calculated_basis_functions[order][_spin];};

		std::vector<double> const& getBeforeCompositions() const{ return compositions_before; };
		// std::vector<double> getBeforeSpins(){ return spins_before; };
		std::vector<std::vector<double>> const& getBeforeCorrelationFunctions() const{ return correlation_functions; };
		std::shared_ptr<basisfunc> const& getBasisFunc() const{ return pbasis_functions;}
		std::shared_ptr<indexorders> const& getIndexOrders() const{ return pindex_orders;}
		std::shared_ptr<allclusters> const& getAllClusters() const{ return pall_clusters;}
		std::mt19937 const& getMt() const{ return mt;}
		std::function<int ()> const& getRndIntN() const{ return rnd_int_N;}
		std::function<double ()> const& getRndReal() const{ return rnd_real;}

};

#endif
