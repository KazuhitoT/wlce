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

using allclusters = std::vector<std::vector<std::vector<std::vector<int>>>>;
using labels      = std::vector<std::pair<int, Eigen::Vector3d>>;

class ParseClusterOut {
	public:
		ParseClusterOut(const char* f, std::vector<int> index = {});
		std::shared_ptr<allclusters> getCluster() const {return this->pall_clusters;}
		std::shared_ptr<labels> getLabel() const {return this->plabels;}

	private:
		std::ifstream ifs;
		const char* filename;
		std::shared_ptr<allclusters> pall_clusters;
		std::shared_ptr<labels> plabels;
};


class ParseEcicar{
public:
	ParseEcicar(const char*);
	~ParseEcicar(){};

	std::vector<int>     getIndex() 			const {return index;};
	std::vector<double> 	getEci(int index)	{return ecicar[index];};
	std::map<int, std::vector<double>> 	getEci()			const {return ecicar;};

private:
	std::ifstream ifs;
		std::map<int /*index*/ ,std::vector<double> /*eci*/> ecicar;
		std::vector<int> index;
};

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

	void expandPoscar(int, int, int);

	Eigen::Matrix3d getLatticeBasis() const;
	std::vector<int> getAtomTypes() const;
	std::vector<std::pair<int, Eigen::Vector3d>> getAtoms() const;
	std::string getComment() const {return comment;}
	std::string getCoordinateType() const {return coordinate_type;}
};

using labels = std::vector<std::pair<int, Eigen::Vector3d>>;

#endif
