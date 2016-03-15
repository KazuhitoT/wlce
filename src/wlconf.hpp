#ifndef WLCONF_HPP
#define WLCONF_HPP

#include "./Conf2corr.hpp"

class WLconf : public Conf2corr{
	private:
		double emin, emax, edelta;
		int bin, index;

		std::vector<double> chemical_potential;

	public:
		WLconf(char* filename,
			 std::vector<double> _spinposcar,
			 std::vector<double> _spince,
			 std::shared_ptr<allclusters> _pall_clusters,
			 std::map<int /*index*/ , double /*eci*/> ecicar,
			 std::shared_ptr<basisfunc>   _pbasis_functions = nullptr,
			 std::shared_ptr<indexorders> _pindex_orders = nullptr,
			 bool is_exchange = false,
			 std::vector<double> _chemical_potential = std::vector<double>()
		 ):Conf2corr(filename, _spinposcar, _spince, _pall_clusters, _pbasis_functions, _pindex_orders){}

		WLconf(char* filename,
			std::shared_ptr<Input> _in,
			std::shared_ptr<allclusters> _pall_clusters,
			std::map<int /*index*/ , double /*eci*/> ecicar,
			std::shared_ptr<basisfunc>   _pbasis_functions = nullptr,
			std::shared_ptr<indexorders> _pindex_orders = nullptr,
			bool is_exchange = false
		);



		// WLconf(char* filename,
		// 	 vector<double> _spinposcar = vector<double>(),
		// 	 vector<double> _spince = vector<double>(),
		// 	 Ensemble ensemble = Ensemble::ALLOY,
		// 	 vector<double> _chemical_potential = vector<double>()
		//  	):Conf2corr(filename, _spinposcar, _spince, ensemble, _chemical_potential){
		// 	if( in == nullptr )
		// 		WLconf::in.reset(new WLInput());
 	// 		emin       = in->emin;
 	// 		emax       = in->emax;
 	// 		if(in->corrDosIndex > 0){
 	// 			emin   = -1;
 	// 			emax   = 1;
 	// 		}
 	// 		bin        = in->bin;
 	// 		edelta     = (emax - emin) / (double)bin;
 	// 		corrdelta  = 2 / (double)bin;
		//  };
		//
		// WLconf(vector<double> spins) : Conf2corr(spins){};    /* スピンベクトルを渡した場合にコピーする*/
		// WLconf(const WLconf &obj) : Conf2corr(obj){
		// 	this->bin_index = obj.getIndex();
		// }
		//
		// WLconf &operator=(const WLconf &obj){
		// 	if (this == &obj)
		// 		return *this;
		// 	Conf2corr::operator=(obj);
		// 	this->bin_index = obj.getIndex();
		// 	return *this;
		// }
		//
		// static double   getEmax(){return emax;}
		// static double   getEmin(){return emin;}
		// static double   getEdelta(){return edelta;}
		// static int      getBin(){return bin;}
		// static double   getCorrDelta(){return corrdelta;}
		//
		// double getProperty(int a = -1){ return this->getTotalEnergy();}
		//
		// ~WLconf(){};
		//
		// /*  macroscopic states 用  */
		// void setEquivalentMacroState();
		// // void setCovarianceMatrixForAllMacroStates(string);
		// // void setCovarianceMatrixInBin(int, string, bool is_thisspin=false);
		//
		// void setSpinsRandom();
		// void setProperty(int corrDosIndex=-1){
		// 	if(corrDosIndex > 0)
		// 		this->setTotalEnergy();
		// 	else
		// 		this->setTotalEnergy();
		// }
		//
		// void setMemento(){
		// 	Conf2corr::setMemento();
		// 	this->bin_index_before = this->bin_index;
		// }
		// void Memento(){
		// 	Conf2corr::Memento();
		// 	this->bin_index = this->bin_index_before;
		// }
		//
		// void setIndex(int corrDosIndex = -1){
		// 	double tmp_index;
		// 	// if( corrDosIndex > 0 )
		// 	// 	tmp_index = (this->getCorrelationFunctions().at(corrDosIndex) - (-1)) / this->getCorrDelta();
		// 	// else
		// 		tmp_index = (this->getTotalEnergy() - this->getEmin()) /  this->getEdelta();
		// 	this->bin_index = this->indexValidation(tmp_index);
		// }
		// void setIndex(double property){
		// 	double tmp_index = (property - emin) / edelta;
		// 	this->bin_index = this->indexValidation(tmp_index);
		// }
		//
		// int indexValidation(int _index){
		// 	if( _index<0 )
		// 		_index = 0;
		// 	else if( _index>(bin-1) )
		// 		_index = bin-1;
		// 	return _index;
		// };
		//
		// int             getIndex()			 const {return bin_index;};
		//
		// vector<int> getNeglectBinIndex();
		// void dispCovarianceCoefficient();
		// void setNewConf();
		// void setNewConf(WLconf&);
		// void setNewConfInBin();

};


#endif
