#ifndef __PARSER_HPP
#define __PARSER_HPP

#include <string>
#include <map>
#include <cmath>
#include <ctime>
#include <typeinfo>
#include <fstream>
#include <iostream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <memory>

// #define DEBUG

class Parser{
	private:
		std::ifstream ifs;
		const char* filename;
		std::vector<std::vector<double> > content;
	public:
		Parser(const char*);						/* ecicarとか */
		Parser(const char*, const std::vector<int>&);	/* index指定した行のみ読み込み */
		~Parser(){};
		std::vector<std::vector<double> > getContent()      const { return content; };
		std::vector<double>          getContent(int i) const { return content[i]; };
		double			        getContent(int i, int j) const { return content[i][j]; };

		void clearContent() { std::vector<std::vector<double> > ().swap(content); };
};


class ParseEcicar : public Parser{
	private:
		std::map<int /*index*/ , double /*eci*/> ecicar;
		std::vector<int> index;
	public:
		// ParseEcicar(){};
		ParseEcicar(const char*);
		~ParseEcicar(){};

		std::vector<int>     getIndex() 			const {return index;};
		int 			getEci(int index)	{return ecicar[index];};
		std::map<int, double> 	getEci()			const {return ecicar;};
};

using allclusters = std::vector<std::vector<std::vector<std::vector<int>>>>;
class ParseMultiplicityIn : public Parser{
	private:
		std::map<int /* = index */, std::pair<int /* = numPoints */, int /* = numClusters */> > multiplicity;
	public:
		// ParseMultiplicityIn(){};
		ParseMultiplicityIn(const char*, const std::vector<int>&);
		~ParseMultiplicityIn(){};
		std::pair<int, int> getMultiplicityIn(int index) const {return multiplicity.at(index);};
		std::map<int, std::pair<int, int> > getMultiplicityIn() const {return multiplicity;};
};

class ParseClusterIn : public Parser{
	private:
		// ClusterOut[index][lattice_point]~
		std::shared_ptr<allclusters> pclusters;
		std::vector<int> index;
	public:
		ParseClusterIn(const char*, const std::vector<int>&, const std::map<int , std::pair<int, int> >&);
		std::shared_ptr<allclusters> getCluster() const {return pclusters;};

 		/* ちゃんとパースできてるかの確認*/
		void checkClusterIn();
};

#endif
