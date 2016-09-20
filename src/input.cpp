#include "input.hpp"

std::string Input::getDataByString(std::string data_key, std::string label){
	if( boost::optional<std::string> value = pt.get_optional<std::string>(label+"."+data_key) ) {
		return value.get();
	} else {
		return "";
	}
}
