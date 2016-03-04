// read poscar.in
// command line argument -2=10 -3=11

// requiremnt Eigen boost
#include <unordered_map>
#include <memory>
#include <iterator>
#include <Eigen/Core>
#include "./parsePoscar.hpp"
extern "C" {
	#include "../spglib/src/spglib.h"
}

#define MAX_BODY 6

constexpr int max_symmetry_num = 48;

#define DEBUG

void outputPoscar(Eigen::Matrix3d lattice, std::vector<Eigen::Vector3d> position, int N, std::string prefix){
	std::ofstream ofs("poscar."+prefix);
	ofs << "POSCAR" << std::endl;
	ofs << "1.0" << std::endl;
	ofs.setf(std::ios::right, std::ios::adjustfield);
	ofs.setf(std::ios_base::fixed, std::ios_base::floatfield);
	ofs << std::setprecision(15);
	ofs << std::setw(20);

	for(int i=0; i<3; ++i){
		for(int j=0; j<3; ++j){
			ofs << lattice(i,j) << " ";
		}
		ofs << std::endl;
	}
	ofs << N << std::endl;
	ofs << "Direct" << std::endl;
	for(int i=0; i<N; ++i){
		for(int j=0; j<3; ++j){
			ofs << position[i](j) << " ";
		}
		ofs << std::endl;
	}
};

void outputSymmetry(int rotation[][3][3], double translation[][3], int num_sym){
	std::ofstream ofs("sym.out");
	for( int i=0; i<num_sym; ++i){
		ofs << "--- sym " << i+1 << " --- " << std::endl;
		for( int j=0; j<3; ++j){
			for( int k=0; k<3; ++k){
				ofs << rotation[i][j][k] << " ";
			}
			ofs << std::endl;
		}
		assert( translation[i][0]==0 );
		ofs << translation[i][0] << " ";
		ofs << translation[i][1] << " ";
		ofs << translation[i][2] << " ";
		ofs << std::endl;
	}
}

Eigen::Vector3d validDistance(const Eigen::Vector3d& lhs, const Eigen::Vector3d& rhs){
	Eigen::Vector3d d_vec = lhs - rhs;
	for( int i=0; i<3; ++i ){ /* boundary condition */
		if( d_vec(i) >0.5 )        {d_vec(i) = 1-d_vec(i);}
		else if( d_vec(i) < -0.5 ) {d_vec(i) = 1+d_vec(i);}
	}
	return d_vec;
};

int main(int argc, char* argv[]){
	double d2, d3, d4 = 0;
	// for (int i=1;i<argc;i++) {
		// if( argv[i][1] == '2' )
		// if (!strcmp(argv[i],"-S")) ensemble=Ensemble::SPIN;
		// else if(!strcmp(argv[i],"-p"))  profile=true;
		// else if(!strcmp(argv[i],"-cp")) corr_profile=true;
		// else if(!strcmp(argv[i],"-d")) debug=true;
		// else{
		// 	cout << "ERROR command line option [" << argv[i] << "] does not exist." << endl;
		// 	exit(1);
		// }
	// }

	#ifdef DEBUG
	d2 = 3;
	d3 = 5;
	d4 = 5;
	#endif

	ParsePoscar poscar("poscar.in");
	const int N = poscar.getAtoms().size();

	double lattice[3][3];
	double lattice_prim[3][3];
	double lattice_unit[3][3];
	Eigen::Map<Eigen::Matrix3d>(&(lattice[0][0]), 3, 3)      = poscar.getAxis();
	Eigen::Map<Eigen::Matrix3d>(&(lattice_prim[0][0]), 3, 3) = poscar.getAxis();
	Eigen::Map<Eigen::Matrix3d>(&(lattice_unit[0][0]), 3, 3) = poscar.getAxis();

	double position[N][3];
	double position_prim[N][3];
	double position_unit[N][3];
	const auto atoms = poscar.getAtoms();
	for( int i=0; i<N; ++i ){
		for( int j=0; j<3; ++j ){
			position[i][j]      = atoms[i].second[j];
			position_prim[i][j] = position[i][j];
			position_unit[i][j] = position[i][j];
		}
 	}

	int count = 0;
	int type  = 1;
	int types[N];
	for( const auto& i : poscar.getNumAtoms() ){
		for( int j=0; j<i; ++j){
			types[count++] = type;
		}
		++type;
	}
	assert( N == count);

	int prim_N = spg_standardize_cell(lattice_prim, position_prim, types, N, 1, 0, 0.00001);
	outputPoscar(lattice_prim, position_prim, prim_N, "prim");

	const int max_size_op = prim_N * max_symmetry_num;
	int rotation[max_size_op][3][3];
	double translation[max_size_op][3];
	int prim_types[1] = {1};
	const int num_sym = spg_get_symmetry(rotation, translation, max_size_op, lattice_prim, position_prim, prim_types, prim_N, 0.0001);
	outputSymmetry(rotation, translation, num_sym);

	int N_unit = spg_standardize_cell(lattice_unit, position_unit, types, N, 0, 0, 0.00001);
	outputPoscar(lattice_unit, position_unit, N_unit, "unit");

	int expand_x = 1;
	int expand_y = 1;
	int expand_z = 1;
	Eigen::Matrix3d lattice_unit_matrix = Eigen::Map<Eigen::Matrix3d>(&lattice_unit[0][0]);
	{
		Eigen::Vector3d unit_x, unit_y, unit_z;
		unit_x << 1, 0, 0;
		unit_y << 0, 1, 0;
		unit_z << 0, 0, 1;
		double unit_length_x = (lattice_unit_matrix * unit_x).norm();
		double unit_length_y = (lattice_unit_matrix * unit_y).norm();
		double unit_length_z = (lattice_unit_matrix * unit_z).norm();
		while( unit_length_x <= (d2*2.) ){
			unit_x[0] = ++expand_x;
			unit_length_x = (lattice_unit_matrix * unit_x).norm();
		}
		while( unit_length_y <= (d2*2.) ){
			unit_y[1] = ++expand_y;
			unit_length_y = (lattice_unit_matrix * unit_y).norm();
		}
		while( unit_length_z <= (d2*2.) ){
			unit_z[2] = ++expand_z;
			unit_length_z = (lattice_unit_matrix * unit_z).norm();
		}
	}
	std::cout << " poscar.expand = ";
	std::cout << expand_x << " * " << expand_y << " * " << expand_z << "  ";
	std::cout << "poscar.unit" << std::endl;


	std::vector<Eigen::Vector3d> position_ex;
	Eigen::Matrix3d lattice_ex;
	{
		Eigen::Matrix3d div_mat;
		div_mat << 1./double(expand_x),0,0,
							 0,1./double(expand_y),0,
							 0,0,1./double(expand_z);
		/* !! Eigen::Map<Eigen::MatrixXd>(&(position_unit[0][0]), N_unit,3) does not work */
		std::vector<Eigen::Vector3d> position_unit_vec ;
		for(int i=0; i<N_unit; ++i){
			position_unit_vec.push_back( div_mat * Eigen::Map<Eigen::Vector3d>(&(position_unit[i][0]), 3));
		}
		for(int x=0; x<expand_x; ++x){
			for(int y=0; y<expand_y; ++y){
				for(int z=0; z<expand_z; ++z){
					Eigen::Vector3d block_vec;
					block_vec << double(x)/double(expand_x), double(y)/double(expand_y), double(z)/double(expand_z);
					for(int i=0; i<N_unit; ++i){
						position_ex.push_back( block_vec + position_unit_vec[i]);
					}
				}
			}
		}
		Eigen::Matrix3d max_ex_mat;
		max_ex_mat << expand_x,0,0, 0,expand_y,0, 0,0,expand_z;
		lattice_ex = max_ex_mat * lattice_unit_matrix;
	}
	outputPoscar(lattice_ex, position_ex, position_ex.size(), "expand");

	// std::cout << lattice_ex << std::endl;
	// std::cout << d_vec.transpose() << std::endl;


	int nbody = 0;

	nbody = 2;
	std::unordered_map<double, Eigen::Vector3d> distance_atom;
	for( int i=1; i<position_ex.size(); ++i ){
		Eigen::Vector3d d_vec = validDistance(position_ex[i], position_ex[0]);
		double distance = (lattice_ex * d_vec).norm();
		if( distance <= d2 ){
			distance_atom[distance] = d_vec;
			auto itr = distance_atom.find(distance);
			if( itr == distance_atom.end() ) {distance_atom[distance] = d_vec; }
		}
	}

	//
	// for(auto itr = distance_atom.begin(); itr != distance_atom.end(); ++itr) {
	// 	std::cout << "key = " <<  distance_atom.at(itr->first) << std::endl;
	// 	std::cout << "key = " << itr->first << ", val = " << itr->second << "\n";    // 値を表示
	// }

	std::vector<Eigen::Matrix3d> rotation_matrix;
	std::vector<Eigen::Vector3d> translation_vector;
	for(int i=0; i<max_size_op; ++i) {
		double tmp_rot[3][3];
		double tmp_tra[3];
		for( int tmp_i=0; tmp_i<3; ++tmp_i ){
			for( int tmp_j=0; tmp_j<3; ++tmp_j ){
				tmp_rot[tmp_i][tmp_j] = double(rotation[i][tmp_i][tmp_j]);
			}
			tmp_tra[tmp_i] = translation[i][tmp_i];
		}
		rotation_matrix.push_back(Eigen::Map<Eigen::Matrix3d>(&(tmp_rot[0][0])));
		translation_vector.push_back(Eigen::Map<Eigen::Vector3d>(&(tmp_tra[0])));
	}
		// exit(1);
	// std::unordered_map<double, std::vector<std::pair<int, Eigen::Vector3d>>> clusterlist2;
	std::vector<std::vector<int>> clusterlist2;
	for( int i=0; i<position_ex.size(); ++i ){
		std::vector<int> clusterlist2_tmp;
		for( auto base_atom : distance_atom ) {
			for( int num_op=0; num_op<max_size_op; ++num_op){
				Eigen::Vector3d target = rotation_matrix[num_op] * distance_atom[base_atom.first];
				for( int j=0; j<position_ex.size(); ++j ){
					// std::cout << (position_ex[j] - position_ex[i]) << std::endl;
					if( target == (position_ex[j] - position_ex[i]) ){
						clusterlist2_tmp.push_back(j);
						// clusterlist2.push_back
						// clusterlist2[i.first] =
					}
				}
			}
		}
		clusterlist2.push_back(clusterlist2_tmp);
	}


	std::cout << clusterlist2.size() << std::endl;
	std::cout << clusterlist2.size() << std::endl;
	for( int i=0; i<	clusterlist2.size(); ++i ){
		for( int j=0; j<	clusterlist2[i].size(); ++j ){
			std::cout << clusterlist2[i][j];
		}
		std::cout << std::endl;
	}

	// for(auto itr = distance_atom.begin(); itr != distance_atom.end(); ++itr) {
	// 	std::cout << "key = " << itr->first << ", val = " << itr->second << "\n";    // 値を表示
	// 	std::cout << "key = " << itr->first << ", val = " << itr->second << "\n";    // 値を表示
	// }

	//
	// for( int i=0; i<3; ++i ){
	// 	for( int j=0; j<3; ++j ){
	// 		std::cout << lattice[i][j] << " ";
	// 	}
	// 	std::cout << std::endl;
 // 	}
	//
	// for( int i=0; i<N; ++i ){
	// 	std::cout << types[i] << " ";
	// }
	// std::cout << std::endl;
	//
	// for( int i=0; i<N; ++i ){
	// 	for( int j=0; j<3; ++j ){
	// 		std::cout << position[i][j] << " ";
	// 	}
	// 	std::cout << std::endl;
 // 	}


	// std::cout << num_sym << std::endl;

}
