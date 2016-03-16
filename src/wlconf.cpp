#include "./wlconf.hpp"

WLconf::WLconf(char* filename,
	std::shared_ptr<Input> _in,
	std::shared_ptr<allclusters> _pall_clusters,
	const std::map<int /*index*/ , std::vector<double> /*eci*/>& _ecicar,
	std::shared_ptr<basisfunc>   _pbasis_functions,
	std::shared_ptr<indexorders> _pindex_orders,
	bool _is_exchange
):Conf2corr(filename, _in, _pall_clusters, _pbasis_functions, _pindex_orders),
	is_exchange(_is_exchange){

	setEci(_ecicar);

	_in->setData("EMIN", this->emin, true);
	_in->setData("EMAX", this->emax, true);
	_in->setData("BIN",  this->bin, true);

	_in->setData("MCSTEP", mcstep, true);
	_in->setData("FLATCHECKSTEP", flatcheck_step, true);
	_in->setData("LOGFACTOR", logfactor, true);
	_in->setData("LOGFLIMIT", logflimit, true);
	_in->setData("FLATCRITERION", flat_criterion, true);

	_in->setData("CHEMIPOT", this->chemical_potential);

	_in->setData("SPININPUT", input_spin_filename);
	_in->setData("SETRANDOM", setrandom);

	if( emin >= emax ){
		std::cerr << "ERROR : input EMIN > EMAX " << std::endl;
		exit(1);
	}

	if( input_spin_filename.size()>0 ){
		this->setSpinsFromDat();
		this->setInitialCorrelationFunction();
		// index_neglect_bin = PoscarSpin.getNeglectBinIndex();
	}
	// else if( setrandom>0 ){
	// 	this->setSpinsRandom();
	// }

	this->edelta = (emax - emin) / (double)bin;
}

void WLconf::dispInput(void){
	std::cout << "EMIN          : " << emin       << std::endl;
	std::cout << "EMAX          : " << emax       << std::endl;
	std::cout << "BIN           : " << bin        << std::endl;

	std::cout << "MCSTEP        : " << mcstep       << std::endl;
	std::cout << "FLATCHECKSTEP : " << flatcheck_step  << std::endl;
	std::cout << "LOGFACTOR     : " << logfactor       << std::endl;
	std::cout << "FLATCRITERION : " << flat_criterion       << std::endl;

	if( this->chemical_potential.size()>0 ) {
		std::cout << "CHEMIPOT      : ";
		for(auto i : this->chemical_potential )
			std::cout << i << " ";
		std::cout << std::endl;
	}

	if( this->input_spin_filename.size()>0 )
		std::cout << "SPININPUT      : " << input_spin_filename << std::endl;
	if( this->setrandom>0 )
		std::cout << "SETRANDOM      : " << setrandom << std::endl;

}

bool WLconf::setSpinsFromDat(){
	std::ifstream ifs(input_spin_filename.c_str());

	if(!ifs){
		ifs.close();
		std::cout << "ERROR : flle [" << input_spin_filename << "] does not exist." << std::endl;
		exit(1);
	}

	bool isSpinSearched = false;
	std::string buf;
	std::vector<std::string> v;
	while(ifs && getline(ifs, buf)){
		if(buf.size() == 0) continue;
		boost::algorithm::trim(buf);
		boost::algorithm::split(v, buf, boost::is_space());
		double e = stod(v.at(1));
		v.erase(v.begin());
		v.erase(v.begin());

		if(v.size() != this->getSpins().size()){
			std::cout << "ERROR : input spin size in [" << input_spin_filename << "] differs from [POSCAR.spin]" << std::endl;
			exit(1);
		}
		if(emin <= e and e <= emax){
			for(int i=0; i<v.size(); ++i){
				this->setSpins(i,stod(v.at(i)));
			}
			isSpinSearched = true;
			break;
		}
	}
	ifs.close();

	if(!isSpinSearched){
		std::cout << "ERROR : no spin configuration satisfies the condition in [wang-landau.ini]" << std::endl;
		exit(1);
	}

	this->setSpinsBefore(this->getSpins());
	return true;
}

std::vector<int> WLconf::getNeglectBinIndex(){
	if( input_spin_filename.size() == 0 ) return std::vector<int>();

	std::vector<double> e_vec;
	std::vector<bool>   is_bin_exist_vec(this->bin, false);

	std::ifstream ifs(input_spin_filename);
	if(!ifs){
		std::cerr << "ERROR : flle [" << input_spin_filename << "] does not exist." << std::endl;
		exit(1);
	}

	std::string buf;
	std::vector<std::string> v;
	while(ifs && getline(ifs, buf)){
		if(buf.size() == 0) continue;
		boost::algorithm::trim(buf);
		boost::algorithm::split(v, buf, boost::is_space());
		double e = stod(v.at(1));
		e_vec.push_back(e);
	}
	ifs.close();

	for(auto e : e_vec){
		int index = (e - this->emin) / this->edelta;
		if( index<0 ) index = 0;
		else if(index>(bin-1)) index = bin-1;
		is_bin_exist_vec.at(index) = true;
	}

	std::vector<int> neglect_bin_index;
	for(int i=0; i<is_bin_exist_vec.size(); ++i){
		if(!is_bin_exist_vec.at(i)) neglect_bin_index.push_back(i);
	}

	return neglect_bin_index;
}

void WLconf::setEci(const std::map<int ,std::vector<double>>& ecicar){
	if(this->eci.size() > 0) return ;
	for(const auto& i : ecicar){
		this->eci.push_back(std::pair<int, std::vector<double>>(i.first, i.second));
	}
}

void WLconf::setTotalEnergy(){
	totalEnergy = 0;
	for(int i=0, imax=eci.size(); i<imax; ++i){
		for(int j=0, jmax=eci[i].second.size(); j<jmax; ++j){
			totalEnergy += this->getCorrelationFunctions(i,j) * eci[i].second[j];
			// std:: cout << i <<" " << this->getCorrelationFunctions(i,j) << " " << eci[i].second[j] << " " << totalEnergy << std::endl;
		}
	}
	/*  setCorrに組み込んだほうが早い  */
	if( chemical_potential.size()>0 ){
		std::vector<double> compositions;
		std::vector<double> spins = this->getSpins();
		for(const auto& spin : this->getSpinCE()){
			 double composition = std::count(spins.begin(), spins.end(), spin) / double(spins.size());
			 compositions.push_back(composition);
		}
		assert( compositions.size() == chemical_potential.size() );
		for(int i=0; i<compositions.size(); ++i){
			totalEnergy -= chemical_potential[i] * compositions[i];
		}
	}
}

void WLconf::setNewConf(){
	int count              = 0;
	double Econf           = this->getTotalEnergy();
	this->setIndex();
	int index_conf         = this->getIndex();
	const int index_before = index_conf;
	double m_distance_before = 0;
	double m_distance_after  = 0;
	while(index_conf == index_before || Econf > this->emax || Econf < this->emin){
		this->setCorrelationFunction();
		this->setTotalEnergy();
		// std::cout << this->getTotalEnergy() << std::endl;
		this->setIndex();
		Econf = this->getTotalEnergy();
		index_conf  = this->getIndex();

		if( index_conf != index_before and Econf <= this->emax and Econf >= this->emin )
			break;

		if( Econf > this->emax || Econf < this->emin ) {
			this->Memento();
		}

	}
}
