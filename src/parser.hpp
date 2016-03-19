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
#include <Eigen/Core>
#include <memory>
#include <random>
#include <iomanip>

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
		std::map<int /*index*/ ,std::vector<double> /*eci*/> ecicar;
		std::vector<int> index;
	public:
		// ParseEcicar(){};
		ParseEcicar(const char*);
		~ParseEcicar(){};

		std::vector<int>     getIndex() 			const {return index;};
		std::vector<double> 	getEci(int index)	{return ecicar[index];};
		std::map<int, std::vector<double>> 	getEci()			const {return ecicar;};
};

using allclusters = std::vector<std::vector<std::vector<std::vector<int>>>>;
class ParseMultiplicityIn : public Parser{
	private:
		std::map<int /* = index */, std::pair<int /* = numPoints */, int /* = numClusters */> > multiplicity;
	public:
		// ParseMultiplicityIn(){};
		ParseMultiplicityIn(const char*);
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
		ParseClusterIn(const char*, const std::map<int , std::pair<int, int> >& multiplicity);
		ParseClusterIn(const char*, const std::vector<int>&, const std::map<int , std::pair<int, int> >&);
		std::shared_ptr<allclusters> getCluster() const {return pclusters;};

 		/* for Confirm */
		void checkClusterIn();
};


void outputPoscar(double lattice[3][3], double position[][3], int N, std::string prefix);

class ParsePoscar {
private:
	std::string comment;
	double unit;
	Eigen::Matrix3d lattice_basis;
	std::vector<int> atom_types;
	std::vector<std::pair<int, Eigen::Vector3d>> atoms;
	std::string coordinate_type;

public:
	ParsePoscar(){};
	ParsePoscar(const char*);

	Eigen::Matrix3d getLatticeBasis();
	std::vector<int> getAtomTypes();
	std::vector<std::pair<int, Eigen::Vector3d>> getAtoms();
	std::string getComment() const {return comment;}
	std::string getCoordinateType() const {return coordinate_type;}
};

#endif
