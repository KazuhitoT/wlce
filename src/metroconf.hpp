#ifndef METROCONF_HPP
#define METROCONF_HPP

#include "./conf2corr.hpp"

class Metroconf : public Conf2corr{
	private:
		double totalEnergy, totalEnergy_before;

		std::vector<std::pair<int, std::vector<double>>> eci;
		std::vector<double> chemical_potential;

	public:
		Metroconf(char* filename,
			 std::vector<double> _spinposcar,
			 std::vector<double> _spince,
			 std::shared_ptr<labels> _plabels,
			 std::shared_ptr<allclusters> _pall_clusters,
			 std::map<int /*index*/ , std::vector<double> /*eci*/> ecicar,
			 std::shared_ptr<basisfunc>   _pbasis_functions = nullptr,
			 std::shared_ptr<indexorders> _pindex_orders = nullptr
		 ):Conf2corr(filename, _spinposcar, _spince, _plabels, _pall_clusters, _pbasis_functions, _pindex_orders){
			 setEci(ecicar);
		 }

		Metroconf(char* filename,
			std::shared_ptr<Input> _in,
			std::shared_ptr<labels> _plabels,
			std::shared_ptr<allclusters> _pall_clusters,
			const std::map<int /*index*/ , std::vector<double> /*eci*/>& ecicar,
			std::shared_ptr<basisfunc>   _pbasis_functions = nullptr,
			std::shared_ptr<indexorders> _pindex_orders = nullptr
		);

		void dispInput();

		void setEci(const std::map<int, std::vector<double>>& ecicar);
		void setTotalEnergy();
		void setCorrelationFunction(){
			if( this->chemical_potential.size()>0 )
				 this->setCorrelationFunction_flip();
			else
				this->setCorrelationFunction_exchange();
		};

		virtual void setNewConf();
		virtual void setMemento(){
			Conf2corr::setMemento();
			this->totalEnergy_before = this->totalEnergy;
		}
		virtual void Memento(){
			Conf2corr::Memento();
			this->totalEnergy = this->totalEnergy_before;
		}

		double getTotalEnergy(){ return totalEnergy;}
		std::vector<double> getChemicalPotential(){ return chemical_potential;}
		void outputEnergySpin(int, std::string);

};


#endif
