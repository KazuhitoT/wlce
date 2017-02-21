#include <regex>
#include <cstdlib>
#include <dirent.h>

#include "./wlconf.hpp"

namespace WLconf{

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

	void outputHistogram(const std::vector<double>& dos, const std::vector<double>& histogram, double emin, double delta, int fstep, std::string prefix){
			std::ofstream ofs("out-wl"+prefix+std::to_string(fstep)+".dat");
			ofs.setf(std::ios_base::fixed, std::ios_base::floatfield);
			for(int i=0, imax=dos.size(); i<imax; ++i){
				ofs << i << " " << std::setprecision(10) << emin + i*delta << " " << emin + (i+1)*delta << " " <<
				std::setprecision(10) << dos[i] << " " << histogram[i] << std::endl;
			}
			ofs.close();
	}

	//  NOTE available format : out-wl*.dat
	void restart(unsigned int& num_fstep, double& logfactor, std::vector<double>& dos, std::vector<int>& index_neglect_bin){

		DIR *dp;
		dirent* entry;

		dp = opendir("./");
		if (dp==NULL) exit(1);

		entry = readdir(dp);
		std::string filename;

		while( entry != NULL ){

			std::string tmp_filename_restart = entry->d_name;

			std::regex re( R"(out\-wl(\d+)\.dat)" );
			std::smatch m ;

			if( std::regex_match( tmp_filename_restart, m, re ) ){

				auto result = std::regex_replace( tmp_filename_restart, re, "$1" ) ;

				int n = std::atoi(result.c_str());
				if( n >= num_fstep ){
					num_fstep = n;
					filename = tmp_filename_restart;
				}

			}

			entry = readdir(dp);

		}

		std::cout << "restart from " << filename << std::endl;
		std::cout << std::endl;

		logfactor *= std::pow(0.5, num_fstep);

		std::ifstream ifs(filename);
		std::string tmp_str;

		for (int i=0; i<dos.size(); ++i){
			int index;
			double emin_bin, emax_bin, entropy, histogram;
			ifs >> index >> emin_bin >> emax_bin >> entropy >> histogram;
			dos[index] = entropy;

			if( entropy > 0 ){
				auto it = find( index_neglect_bin.begin(), index_neglect_bin.end() , index);
				if( it != index_neglect_bin.end() ){
					index_neglect_bin.erase(it);
				}
			}

			std::getline(ifs, tmp_str);
		}

	}


	WLconf::WLconf(char* filename,
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

		_in->setData("EMIN", this->emin, true);
		_in->setData("EMAX", this->emax, true);
		_in->setData("BIN",  this->bin, true);

		_in->setData("SPININPUT", input_spin_filename);
		_in->setData("SETRANDOM", setrandom);

		if( emin >= emax ){
			std::cerr << "ERROR : input EMIN > EMAX " << std::endl;
			exit(1);
		}

		if( input_spin_filename.size()>0 ){
			this->setSpinsFromDat();
			this->setInitialCorrelationFunction();
		}	else if( setrandom>0 ){
			this->setSpinsRandom();
			this->setInitialCorrelationFunction();
		}

		this->edelta = (emax - emin) / (double)bin;
	}

	void WLconf::dispInput(void){
		std::cout << "EMIN          : " << emin       << std::endl;
		std::cout << "EMAX          : " << emax       << std::endl;
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

	void WLconf::setNewConf(){

		int count              = 0;
		double Econf           = this->getTotalEnergy();
		this->setIndex();
		int index_conf         = this->getIndex();
		const int index_before = index_conf;
		while(index_conf == index_before || Econf > this->emax || Econf < this->emin){
			this->setMemento();
			this->setCorrelationFunction();
			this->setTotalEnergy();
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

}
