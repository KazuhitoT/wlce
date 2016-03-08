// TODO : spinからPOSCARを吐き出す関数
// TODO : 多元系に対応させる

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
#include <random>
#include <chrono>
#include <functional>
#include <boost/algorithm/string.hpp>
#include <Eigen/Core>
#include <Eigen/LU>
#include <memory>
#include <cfloat>
#include "./enum.hpp"
#include "./parsePoscar.hpp"
#include "./Parser.hpp"
// #include "./ConfInput.hpp"
#include "./repeated_combination_generator.hpp"

using namespace std;
using namespace Eigen;

class Conf2corr {
	private:
		friend class Conf2corrTest;

		vector<double> spins;
		vector<double> spins_before;
		vector<pair<int, vector<double>>> eci;
		vector<vector<double>> correlation_functions;
		vector<vector<double>> correlation_functions_before;
		vector<vector<double>> radial_distributions;  // for each spins

		double totalEnergy;
		double totalEnergy_before;
		double mahalanobis_distance;
		double ene_mahalanobis_distance;
		double growth_factor;

		const static ParseEcicar ecicar;
		const static ParseDimcar dimcar;
		const static ParseMccar  mccar;
		const static ParseClucar clucar;

		// static unique_ptr<WLInput> in;

		// const static ConfInput *in;

		/*  random number */
		static std::random_device rd;
		static std::mt19937 mt;

		static MatrixXd cov_matrix;
		static MatrixXd cov_matrix_before;
		static VectorXd ave_vector;
		static VectorXd ave_vector_before;
		static int num_allconf;
		/* vec[basis][degree] vec[basis]->polynomials */
		static vector<vector<double>> basis_functions;
		/* vec[num_in_cluster][index][combinations][order] */
		static vector<vector<vector<vector<int>>>> index_orders;

		/* --------- Input parameter for Conf2corr --------- */
		static vector<double> spinposcar;
		static vector<double> spince;
		static vector<double> chemical_potential;
		static Ensemble ensemble;
		static ParsePoscar poscarCluster;
		/* ------------------------------------------------- */

	public:
		Conf2corr(){};

		// Conf2corr(char*);             /* POSCAR.spinをパースする場合 */
		// Conf2corr(char*, int, char**);      /* POSCAR.spinをパースする場合 */

		// Conf2corr(char*);             /* POSCAR.spinをパースする場合 */

		Conf2corr(vector<double>);    /* スピンベクトルを渡した場合にコピーする*/
		Conf2corr(const Conf2corr &obj); /* copy construct */
		Conf2corr(char* filename,
			 vector<double> _spinposcar = vector<double>(),
			 vector<double> _spince = vector<double>(),
			 Ensemble ensemble = Ensemble::ALLOY,
			 vector<double> _chemical_potential = vector<double>(),
			 int _setrandom = -1,
			 std::string _input_spin_filename = ""
		 );
		// Conf2corr(char* filename, vector<double> _spinposcar, vector<double> _spince);
		~Conf2corr(){};
		virtual void setMemento();
		virtual void Memento();

		static std::function<int ()> rnd_int_N;
		static std::function<double ()> rnd_real;

		void setInitialCondition();

		static void setSpinPoscar(vector<double> _spins){spinposcar = _spins;}
		static void setSpinCE(vector<double> _spins){spince = _spins;}
		static void setBasisFunction(){};

		// void setInitialize(char*);
		void setSpins(char* filename);
		void setSpins(vector<double> _spins){ this->spins = _spins; };
		void setSpinsBefore(){ spins = spins_before; };
		void setBeforeSpins(){ spins_before = spins; };
		virtual void setSpinsRandom();
		void setConfRandom();
		bool setSpinsFromDat(string, double /* emin */, double /* emax */, bool isSkipNoSpins = false);
		void setEci();
		void setCorrelationFunctionFromClucar();  	/* 一回だけ行う */
		void setCorrelationFunction_flip(int);       /* mccarから計算する(single spin flip) */
		void setCorrelationFunction_exchange();      /* mccarから計算する(single spin swap)*/
		void setCorrelationFunction();
		void setCorrelationFunction_withMultiAtomSwap();
		void setCorrelationFunctionBefore(){    /* エネルギー指定からズレた時にバックアップに戻す*/
			this->correlation_functions = this->correlation_functions_before;
		};
		void setTotalEnergy();
		virtual void setNewConf();
		void setRadialDistributions();
		void setGrowthFactor();
		static void setBasisCoefficient();
		static double trace(const vector<double>&, const vector<double>&);
		static void setIndexOrders();

		/*  必ずこの順番で  */
		static void setNumAllConf(){++num_allconf;}
		void setAverageVector();
		void updateCovarianceMatrix();
		void setCovarianceMatrix(string);
		void setMahalanobisDistance();			/* 全体用 */
		void setMahalanobisDistance(const int);

		static void setInitializeCov(int num_corr){
			Conf2corr::cov_matrix = MatrixXd::Zero(num_corr, num_corr);
			Conf2corr::ave_vector = VectorXd::Zero(num_corr);
			Conf2corr::ave_vector(0) = 1;
			Conf2corr::ave_vector_before = Conf2corr::ave_vector;
			Conf2corr::num_allconf = 0;
		};

		vector<double>             getSwapSpins();
		vector<int>                getSwapSpinsIndex();
		vector<double>             getSpins()                const {return spins;};
		vector<double>             getBeforeSpins()                const {return spins_before;};
		vector<pair<int, vector<double>>> getEci()                  const {return eci;};
		vector<vector<double>>             getCorrelationFunctions() const {return correlation_functions;};
		vector<vector<double>>             getBeforeCorrelationFunctions() const {return correlation_functions_before;};
		double                     getTotalEnergy()          const {return totalEnergy;};
		double                     getProperty(int);
		double getGrowthFactor(){return this->growth_factor;};
		double getBasisFunction(/* order */  int, /* spince_num */ int);
		vector<vector<double>> getBasisFunctions() const {return basis_functions;}
		bool     isPropertyInside();
		double   getMahalanobisDistance() const {return mahalanobis_distance;};
		double getTotalCorrelationFunctionAbs(int order=1){
			double result = 0;
			for( const auto& corr_funcs : this->correlation_functions ){
				for( const auto& corr_func : corr_funcs ){
					result += pow(abs(corr_func), order);
				}
			}
			return result-1; /* for empty cluster */
		}
		vector<vector<vector<vector<int>>>> getIndexOrders() const {return index_orders;}
		double getNumInCluster(int line_num_dimcar) const { return dimcar.getDimcar(line_num_dimcar).first;}
		double getNumOfCluster(int line_num_dimcar) const { return dimcar.getDimcar(line_num_dimcar).second;}

		static double   getCovarianceMatrix(int i, int j){return Conf2corr::cov_matrix(i, j);}
		static MatrixXd getCovarianceMatrix(){return Conf2corr::cov_matrix;}
		static double   getAverageVector(int i){return Conf2corr::ave_vector(i);}
		static int      getNumAllConf(){return Conf2corr::num_allconf;}

		void dispCorr();
		void dispSpin();
		void dispEci();
		void dispIndexOrder();
		void dispCovarianceCoefficient();
		void dispBasisCoefficient();

		void outputCorrSpin(int num=-1, int mode=0, string filename="corrspin.out");
		void outputCorr(string filename = "corr-profile.out");
		void outputPoscar(string filename = "poscar.out");

  	Conf2corr &operator=(const Conf2corr&);    // 代入演算子
};

#endif
