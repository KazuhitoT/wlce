#include "./wlconfT.hpp"

namespace WLconfT{

	bool checkHistogramFlat(const std::vector<double>& histogram, const std::vector<int>& index_neglect_bin, double lflat, int minstep, double low_cutoff){
		double ave = accumulate(histogram.begin(), histogram.end(), 0) / (double)(histogram.size()-index_neglect_bin.size());
		double limit = ave * lflat ;
		for(int i=0, imax=histogram.size(); i<imax; ++i){

			auto it = find( index_neglect_bin.begin(), index_neglect_bin.end() , i);
			if( it != index_neglect_bin.end() ) continue;

			if( minstep>0 and histogram.at(i)<minstep ) return false;

			if(histogram.at(i) < (ave*(1.0-low_cutoff))) continue;
			else if(histogram.at(i) < limit) return false;
		}
		return true;
	}

	void outputHistogram(const std::vector<double>& dos, const std::vector<double>& histogram, double tmin, double delta, int fstep, std::string prefix){
			std::ofstream ofs("out-wl"+prefix+std::to_string(fstep)+".dat");
			ofs.setf(std::ios_base::fixed, std::ios_base::floatfield);
			for(int i=0, imax=dos.size(); i<imax; ++i){
				ofs << i << " " << std::setprecision(10) << tmin + i*delta << " " << tmin + (i+1)*delta << " " <<
				std::setprecision(10) << dos[i] << " " << histogram[i] << std::endl;
			}
			ofs.close();
	}

	WLconfT::WLconfT(char* filename,
		std::shared_ptr<Input> _in,
		std::shared_ptr<labels> _plabels,
		std::shared_ptr<allclusters> _pall_clusters,
		const std::map<int /*index*/ , std::vector<double> /*eci*/>& _ecicar,
		std::shared_ptr<basisfunc>   _pbasis_functions,
		std::shared_ptr<indexorders> _pindex_orders
	):
		Metroconf(filename, _in, _plabels, _pall_clusters, _ecicar, _pbasis_functions, _pindex_orders),
		input_spin_filename(""),
		setrandom(-1)
		{
		this->index = 0;

		_in->setData("TMIN", this->tmin, true);
		_in->setData("TMAX", this->tmax, true);
		_in->setData("TMOVELIMIT", this->tmovelimit, true);
		_in->setData("BIN",  this->bin, true);

		_in->setData("SPININPUT", input_spin_filename);
		_in->setData("SETRANDOM", setrandom);

		if( tmin >= tmax ){
			std::cerr << "ERROR : input TMIN > TMAX " << std::endl;
			exit(1);
		}

		rnd_int_tmove = std::bind( std::uniform_int_distribution<int>(-this->tmovelimit, this->tmovelimit), mt);

		if( input_spin_filename.size()>0 ){
			this->setSpinsFromDat();
			this->setInitialCorrelationFunction();
		}	else if( setrandom>0 ){
			this->setSpinsRandom();
			this->setInitialCorrelationFunction();
		}

		this->tdelta = (tmax - tmin) / (double)bin;
	}

	void WLconfT::dispInput(void){
		std::cout << "TMIN          : " << tmin       << std::endl;
		std::cout << "TMAX          : " << tmax       << std::endl;
		std::cout << "BIN           : " << bin        << std::endl;

		if( this->getChemicalPotential().size()>0 ) {
			std::cout << "CHEMIPOT      : ";
			for(auto i : this->getChemicalPotential() )
				std::cout << i << " ";
			std::cout << std::endl;
		}

		if( this->input_spin_filename.size()>0 )
			std::cout << "SPININPUT     : " << input_spin_filename << std::endl;
		if( this->setrandom>0 )
			std::cout << "SETRANDOM     : " << "YES" << std::endl;

	}

	bool WLconfT::setSpinsFromDat(){
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
			if(tmin <= e and e <= tmax){
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

	std::vector<int> WLconfT::getNeglectBinIndex(){
		if( input_spin_filename.size() == 0 ) {
			std::vector<int> result;
			for(int i=0; i<this->bin; ++i) result.push_back(i);
			return result;
		}

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
			int index = (e - this->tmin) / this->tdelta;
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

	void WLconfT::setNewConfByEnergy(){
		this->setMemento();
		this->setCorrelationFunction();
		this->setTotalEnergy();
	}


	void WLconfT::setNewConfByTemperature(){

		const int index_before = index;
		while(this->index == index_before || this->index >= 0 || this->index < this->bin){
			this->setMemento();
			this->index += rnd_int_tmove();

			if( this->index != index_before and this->index < this->bin and this->index >= 0 )
				break;

			if( this->index >= this->bin || this->index < 0 ) {
				this->Memento();
			}

		}
	}

}
