#include <boost/foreach.hpp>
#include "input.hpp"

std::string Input::getDataByString(std::string data_key, std::string label){
	if( boost::optional<std::string> value = pt.get_optional<std::string>(label+"."+data_key) ) {
		return value.get();
	} else {
		return "";
	}
}

void Input::disp(std::string label){
	BOOST_FOREACH (const boost::property_tree::ptree::value_type& child, pt.get_child(label)) {
		std::cout << child.first.data() << " : " << child.second.data() << std::endl;
	}
}
