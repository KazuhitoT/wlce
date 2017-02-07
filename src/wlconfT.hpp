#ifndef WLCONFT_HPP
#define WLCONFT_HPP

#include "./metroconf.hpp"

namespace WLconfT{

	bool checkHistogramFlat(const std::vector<double>& histogram, const std::vector<int>& index_neglect_bin, double lflat, int minstep, double low_cutoff=1);
	void outputHistogram(const std::vector<double>& dos, const std::vector<double>& histogram, double tmin, double delta, int fstep, std::string prefix = "");

	class WLconfT : public Metroconf{
		private:
			double temperature, temperature_before, tmin, tmax, tdelta;
			int tmovelimit;
			int setrandom;
			std::string input_spin_filename;

			int bin, index, index_before;

			std::mt19937 mt;
			std::function<int ()> rnd_int_tmove;

		public:

			WLconfT(char* filename,
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

			void setNewConfByTemperature();
			void setNewConfByEnergy();

			double calcTemperature(){ return this->tmin + ((double(index)+0.5) * this->tdelta);}
			void setIndex(){
				index = (this->temperature - tmin) / tdelta;
				this->indexValidation();
			}

			void setMemento(){
				Metroconf::setMemento();
				this->index_before = this->index;
				this->temperature_before = this->temperature;
			}
			void Memento(){
				Metroconf::Memento();
				this->index = this->index_before;
				this->temperature = this->temperature_before;
			}

			void indexValidation(){
				if( index<0 )
					index = 0;
				else if( index>(bin-1) )
					index = bin-1;
			};

			int getIndex(){ return index; }
			int getBeforeIndex(){ return index_before; }
			double getTemperature(){ return this->temperature;}
			double getBeforeTemperature(){ return this->temperature_before;}


	};

}

#endif
