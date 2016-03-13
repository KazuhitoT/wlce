#ifndef PARSEPOSCAR_HPP
#define PARSEPOSCAR_HPP

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <Eigen/Core>

void outputPoscar(double lattice[3][3], double position[][3], int N, std::string prefix);

class ParsePoscar {
private:
	std::vector<int> atom_types;
	std::vector<std::pair<int, Eigen::Vector3d>> atoms;
	Eigen::Matrix3d lattice_basis;
	std::string coordinate_type;
	std::string comment;

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
