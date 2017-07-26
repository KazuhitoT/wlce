#ifndef WLCONF_HPP
#define WLCONF_HPP

#include "./metroconf.hpp"

namespace WLconf{

	bool checkHistogramFlat(const std::vector<double>& histogram, const std::vector<int>& index_neglect_bin, double lflat, int minstep, double low_cutoff=1);
	void outputHistogram(const std::vector<double>& dos, const std::vector<double>& histogram, double emin, double delta, int fstep, std::string prefix = "");
	void restart(unsigned int& num_fstep, double& logfactor, std::vector<double>& dos, std::vector<int>& index_neglect_bin);

	class WLconf : public Metroconf{
		private:
			double emin, emax, edelta;
			int setrandom;
			std::string input_spin_filename;

			int bin, index, index_before;

		public:
			// WLconf(char* filename,
			// 	 std::vector<double> _spinposcar,
			// 	 std::vector<double> _spince,
			// 	 std::shared_ptr<labels> _plabels,
			// 	 std::shared_ptr<allclusters> _pall_clusters,
			// 	 std::map<int /*index*/ , std::vector<double> /*eci*/> ecicar,
			// 	 std::shared_ptr<basisfunc>   _pbasis_functions = nullptr,
			// 	 std::shared_ptr<indexorders> _pindex_orders = nullptr,
			// 	 std::vector<double> _chemical_potential = std::vector<double>()
			//  ):Conf2corr(filename, _spinposcar, _spince, _plabels, _pall_clusters, _pbasis_functions, _pindex_orders),
			//  chemical_potential(_chemical_potential){
			// 	 setEci(ecicar);
			//  }

			WLconf(char* filename,
				std::shared_ptr<Input> _in,
				std::shared_ptr<labels> _plabels,
				std::shared_ptr<allclusters> _pall_clusters,
				const std::map<int /*index*/ , std::vector<double> /*eci*/>& ecicar,
				std::shared_ptr<basisfunc>   _pbasis_functions = nullptr,
				std::shared_ptr<indexorders> _pindex_orders = nullptr
			);

			void dispInput();
			bool setSpinsFromDat();

			std::vector<int> getNeglectBinIndex();

			void setNewConf();
			void setIndex(){
				index = (this->getTotalEnergy() - emin) / edelta;
				this->indexValidation();
			}

			void setMemento(){
				Metroconf::setMemento();
				this->index_before = this->index;
			}
			void Memento(){
				Metroconf::Memento();
				this->index = this->index_before;
			}

			void indexValidation(){
				if( index<0 )
					index = 0;
				else if( index>(bin-1) )
					index = bin-1;
			};

			int getIndex(){ return index; }
			int getBeforeIndex(){ return index_before; }

	};

}

#endif
