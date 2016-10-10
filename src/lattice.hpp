#ifndef _LATTICE_HPP
#define _LATTICE_HPP

#include <Eigen/Core>
#include <vector>
#include <iostream>
#include "./site.hpp"

class Site;

class Lattice {
private :

	const int num_lattice;
	int num_symmetry;
	std::vector<std::shared_ptr<Site>> sites;

	Eigen::Matrix3d basis;
	std::vector<Eigen::Matrix3d> vec_mat3d_rotations;
	std::vector<Eigen::Vector3d> vec_vec3d_translations;

public :

	Lattice(int _num_lattice):num_lattice(_num_lattice) {};

	int getLatticeNum() const { return num_lattice;};
	Eigen::Matrix3d getLatticeBasis() { return basis;};

	void setSite(std::shared_ptr<Site> _site){
		this->sites.push_back(_site);
	};

	std::vector<std::shared_ptr<Site>> getSites(){ return this->sites; }

	void setBasis(double _basis[3][3]){
		this->basis =  Eigen::Map<Eigen::Matrix3d>(&_basis[0][0]);
	}

	void setBasis(const Eigen::Matrix3d& _basis){
		this->basis =  _basis;
	}

	void setRotation(const int _rotation[][3][3]){
			//
			//
			// std::vector<Eigen::Matrix3d> vec_mat_rotations;
			// std::vector<Eigen::Vector3d> vec_vec_translations;
			// for(int i=0; i<num_sym; ++i){
			// 	Eigen::Matrix3d tmp_rotation;
			// 	tmp_rotation <<
			// 			rotation[i][0][0], rotation[i][0][1], rotation[i][0][2],
			// 			rotation[i][1][0], rotation[i][1][1], rotation[i][1][2],
			// 			rotation[i][2][0], rotation[i][2][1], rotation[i][2][2];
			// 	vec_mat_rotations.push_back(tmp_rotation);
			//
			// 	Eigen::Vector3d tmp_translation;
			// 	tmp_translation << translation[i][0], translation[i][1], translation[i][2];
			// 	vec_vec_translations.push_back(tmp_translation);
			// }

	}

	void setTranslation(const double _translation[][3]){

	}
	void setNumSymmetry(const int _num_symmetry){
		this->num_symmetry = _num_symmetry;
	}

	void dispInfo(){
		std::cout << " -- sub lattice : " << this->num_lattice << "--" << std::endl;
		std::cout << " Number of symmetry operations : " << num_symmetry  << std::endl;
	}

};

#endif
