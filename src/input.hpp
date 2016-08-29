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

class Input {
private:
	std::string filename;
	boost::property_tree::ptree pt;

public:
	Input(char* _filename){ read_ini(_filename, pt); filename = _filename;};
	template<typename T> void setData(std::string data_key, T &data, bool is_required=false, std::string label="INPUT");
	template<typename T> void setData(std::string data_key, std::vector<T> &data, bool is_required=false, std::string label="INPUT");
	std::string getDataByString(std::string data_key, std::string label="INPUT");
	void disp(void);
};

template<typename T> void Input::setData(std::string data_key, T &data, bool is_required, std::string label){
	if ( boost::optional<T> value = pt.get_optional<T>(label+"."+data_key)) {
		data = value.get();
	} else if( is_required ) {
		std::cout << data_key + " is required." << std::endl;
		exit(1);
	}
}

template<typename T> void Input::setData(std::string data_key, std::vector<T> &data, bool is_required, std::string label){
	if ( boost::optional<std::string> value = pt.get_optional<std::string>(label+"."+data_key)) {
		boost::trim(*value);
		std::vector<std::string> vec_tmp_str;
		boost::split(vec_tmp_str, *value, boost::is_space());
		for(auto i : vec_tmp_str){
			if( i != "" ) data.push_back(boost::lexical_cast<T>(i));
		}
	} else if( is_required ) {
		std::cout << data_key + " is required." << std::endl;
		exit(1);
	}
}

#endif
