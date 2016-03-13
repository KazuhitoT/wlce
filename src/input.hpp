#ifndef _INPUT_HPP
#define _INPUT_HPP

#include <string>
#include <map>
#include <cmath>
#include <cstdlib>
#include <typeinfo>
#include <fstream>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/optional.hpp>
#include <boost/lexical_cast.hpp>

// using namespace boost;
// using namespace boost::property_tree;
// using namespace std;


class Input {
private:
	boost::property_tree::ptree pt;

public:
	Input(char* filename){
		read_ini(filename, pt);
	};
	template<typename T> void setData(std::string dataName, T &data, bool required=false);
	template<typename T> void setData(std::string dataName, std::vector<T> &data, bool required=false);
	void disp(void);
	//
	// char* ensemble_name   [myenum_util::underlying_value(Ensemble::SIZE)]   = {"SPIN", "ALLOY"};
	// char* property_name   [NULL_PROPERTY]   = {"ENERGY"};
	// char* output_name     [NULL_OUTPUT]     = {"NO OUTPUT", "ENERGY", "ENERGY_AND_SPIN"};
	// int sampling_step, mcstep;
	// double temperature;
	// string input_spin_filename;
	// vector<double> spince, spinposcar, chemical_potential;
	// Ensemble  ensemble     = Ensemble::ALLOY;
	// int  property     = ENERGY;
	// int  output       = OUTPUT_NONE;
	// bool debug        = false;
	// bool profile      = false;
	// bool corr_profile = false;
	//
	// int setrandom = -1;
	// int max_swap = 1;
	//
	// /* for WL */
	// int  method       = CONTINUOUS;
	// int  algorithm    = DEFAULT;
	//
	// int bin, flatcheck,  corrDosIndex_disp;
	// double emin, emax, edelta, factor, flimit;
	//
	// double lflat = 1;
	// double low_cutoff = 0.5;
	//
	// // WLInput(int argc = 0, char** argv = nullptr, string filename="wang-landau.ini");
	// void setCommadLineArgument(int, char**);
	// void setInput();
	// static vector<double> corrIndex2vecIndex(const vector<pair<int, double>>&, const vector<double>&);
	// void corrDosIndex2vecIndex(const vector<pair<int, double>>&);
};

template<typename T> void Input::setData(std::string dataName, T &data, bool required){
	if ( boost::optional<T> value = pt.get_optional<T>("INPUT."+dataName)) {
		data = value.get();
	} else if( required ) {
        std::cout << dataName + " is not valid." << std::endl;
        exit(1);
    }
}

template<typename T> void Input::setData(std::string dataName, std::vector<T> &data, bool required){
	if ( boost::optional<std::string> value = pt.get_optional<std::string>("INPUT."+dataName)) {
		boost::trim(*value);
		std::vector<std::string> vec_tmp_str;
		boost::split(vec_tmp_str, *value, boost::is_space());
		for(auto i : vec_tmp_str){
			if( i != "" ) data.push_back(boost::lexical_cast<T>(i));
		}
	} else if( required ) {
        std::cout << dataName + " is not valid." << std::endl;
        exit(1);
    }
}
#endif
