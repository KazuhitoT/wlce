#include "./metroconf.hpp"

Metroconf::Metroconf(char* filename,
	std::shared_ptr<Input> _in,
	std::shared_ptr<labels> _plabels,
	std::shared_ptr<allclusters> _pall_clusters,
	const std::map<int /*index*/ , std::vector<double> /*eci*/>& _ecicar,
	std::shared_ptr<basisfunc>   _pbasis_functions,
	std::shared_ptr<indexorders> _pindex_orders
):
	Conf2corr(filename, _in, _plabels, _pall_clusters, _pbasis_functions, _pindex_orders)
	{

	setEci(_ecicar);
	_in->setData("CHEMIPOT", this->chemical_potential);

}

void Metroconf::dispInput(void){
	if( this->chemical_potential.size()>0 ) {
		std::cout << "CHEMIPOT      : ";
		for(auto i : this->chemical_potential )
			std::cout << i << " ";
		std::cout << std::endl;
	}
}

void Metroconf::setEci(const std::map<int ,std::vector<double>>& ecicar){
	if(this->eci.size() > 0) return ;
	for(const auto& i : ecicar){
		this->eci.push_back(std::pair<int, std::vector<double>>(i.first, i.second));
	}
}

void Metroconf::setTotalEnergy(){
	totalEnergy = 0;
	for(int i=0, imax=eci.size(); i<imax; ++i){
		for(int j=0, jmax=eci[i].second.size(); j<jmax; ++j){
			totalEnergy += this->getCorrelationFunctions(i,j) * eci[i].second[j];
		}
	}
	if( chemical_potential.size()>0 ){
		for(int i=0; i<this->getCompositions().size(); ++i){
			totalEnergy -= chemical_potential[i] * this->getCompositions(i);
		}
	}
}

void Metroconf::setNewConf(){
	this->setCorrelationFunction();
	this->setTotalEnergy();
}

void Metroconf::outputEnergySpin(int index, std::string filename){
	std::ofstream ofs(filename, std::ios::app);
	ofs << index << " ";
	ofs.precision(10);
	ofs << this->getTotalEnergy() << " ";
	for( auto j : this->getSpins() )  ofs << j << " ";
	ofs << std::endl;
	ofs.close();
}
