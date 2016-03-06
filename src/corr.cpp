// read poscar.in
// command line argument -2=10 -3=11

// requiremnt Eigen boost
#include <unordered_map>
#include <memory>
#include <iterator>
#include <Eigen/Core>
#include "./parsePoscar.hpp"
#include "./site.hpp"
extern "C" {
	#include "../spglib/src/spglib.h"
}

using vector2i = std::vector<std::vector<int>>;
using vector3i = std::vector<std::vector<std::vector<int>>>;
// using vector3i = std::vector<std::vector<std::vector<int>>>;
// using vector3i = std::vector<std::vector<std::vector<int>>>;

#define MAX_BODY 6

constexpr int max_symmetry_num = 48;

#define DEBUG

class hash_vecd {  // ハッシュ関数オブジェクト
public:
  double operator()(const std::vector<double> &x) const {
		std::vector<double> copy = x;
		std::sort(copy.begin(), copy.end());//昇順ソート
    const int C = 997;      // 素数
    double t = 0;
    for (int i = 0; i != x.size(); ++i) {
        t = t * C + copy[i];
    }
    return t;
  }
};

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
		ofs << translation[i][0] << " ";
		ofs << translation[i][1] << " ";
		ofs << translation[i][2] << " ";
		ofs << std::endl;
	}
}

void outputInfoCluster(){
	// std::ofstream ofs("info-cluster.out");
	// for( int i=0; i<num_sym; ++i){
	// 	ofs << "--- sym " << i+1 << " --- " << std::endl;
	// 	for( int j=0; j<3; ++j){
	// 		for( int k=0; k<3; ++k){
	// 			ofs << rotation[i][j][k] << " ";
	// 		}
	// 		ofs << std::endl;
	// 	}
	// 	ofs << translation[i][0] << " ";
	// 	ofs << translation[i][1] << " ";
	// 	ofs << translation[i][2] << " ";
	// 	ofs << std::endl;
	// }
}

// bool is_eq_float(const Eigen::Vector3d& lhs, const Eigen::Vector3d& rhs){
bool is_eq_coordinate(const Eigen::Vector3d& lhs, const Eigen::Vector3d& rhs){
	for(int i=0; i<3; ++i){
		double x = lhs[i];
		if( std::abs(x) < 0.00001 ) x = 0;
		else if( std::abs(1-x) < 0.00001 ) x = 1;

		if( x < 0 ) x+=1;
		else if( x >= 1 ) x-=1;
		if( 0.00001 < std::abs(x-rhs[i]) ){
			return false;
		}
	}
	return true;
}

Eigen::Vector3d validDistance(const Eigen::Vector3d& lhs, const Eigen::Vector3d& rhs){
	Eigen::Vector3d d_vec = lhs - rhs;
	for( int i=0; i<3; ++i ){ /* boundary condition */
		if( d_vec(i) > 0.5 )       {d_vec(i) = 1-d_vec(i);}
		else if( d_vec(i) < -0.5 ) {d_vec(i) = 1+d_vec(i);}
	}
	return d_vec;
};

int main(int argc, char* argv[]){
	double d2, d3, d4 = 0;
	double prec = 0.00001;

	#ifdef DEBUG
	d2 = 3;
	d3 = 5;
	d4 = 5;
	#endif

	d2 = atof(argv[1]);

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
			bool is_found = false;
			for( const auto& i : distance_atom ){
				if( std::abs(i.first-distance) < prec ) {
					is_found = true;
					break;
				}
			}
			auto itr = distance_atom.find(distance);
			if( !is_found ) {distance_atom[distance] = d_vec; }
		}
	}
	// std::shared_ptr<int> ptr2(new int(0));

	std::vector<std::shared_ptr<Site>> site_vec;
	for(int i=0; i<position_ex.size(); ++i) {
		// new std::shared_ptr<Site> ptr(i, position_ex[i]);
		site_vec.push_back( std::shared_ptr<Site>( new Site(i, position_ex[i])) );
	}

	std::ofstream clusters_out( "clusters.out", std::ios::out );

	/*  set nbody = 1 */
	std::cout << " -- point cluster" << std::endl;
	for( int i=0; i<position_ex.size(); ++i){
		clusters_out << i << " ";
	}
	clusters_out << std::endl;

	/*  set nbody = 2 */
	std::cout << " -- 2 body cluster" << std::endl;

	std::unordered_map<double, std::vector<std::vector<int>>> distance_site_to_sites;
	{
		for( const auto& d_a : distance_atom ){
			std::vector<std::vector<int>> a_d_a;
			for(int i=0; i<position_ex.size(); ++i){
				std::vector<int> vec_a;
				for(int j=0; j<position_ex.size(); ++j){
					Eigen::Vector3d d_vec = validDistance(position_ex[i], position_ex[j]);
					double distance = (lattice_ex * d_vec).norm();
					if( std::abs( distance - d_a.first ) < prec ) {
						site_vec[i]->setLinkedSite(d_a.first, site_vec[j]);
						vec_a.push_back(j);
					}
				}
				a_d_a.push_back(vec_a);
			}
			distance_site_to_sites[d_a.first] = a_d_a;
		}
	}

	for( const auto& i : distance_site_to_sites ){
		std::cout << i.first << std::endl;
		for( int j=0; j<i.second.size(); ++j) {
			for( const auto& k : i.second[j] ){
				clusters_out << j << " " << k << " ";
			}
		}
		clusters_out << std::endl;
	}

	/*  set nbody = 3 */
	std::cout << " -- 3 body cluster" << std::endl;
	std::vector<double> d_vec;
	for( const auto& i : distance_site_to_sites ){
		d_vec.push_back(i.first);
		// std::cout << i.first << std::endl;
	}
	std::vector<std::vector<double>> d_b3_index;
	for(const auto& linked_site : site_vec[0]->getLinkedSite() ){
		for(const auto& site2 : linked_site.second ){
			for(const auto& linked_site2 : site2->getLinkedSite() ){
				for(const auto& site3 : linked_site2.second ){
					Eigen::Vector3d d3 = validDistance(site3->getCoordinate(), site_vec[0]->getCoordinate());
					double distance = (lattice_ex * d3).norm();
					if( site3->getLinkedSite(distance, 0) ){
						std::vector<double> d_vec = {linked_site.first, linked_site2.first, distance};
						sort(d_vec.begin(), d_vec.end());
						auto it2 = find( d_b3_index.begin(), d_b3_index.end(), d_vec );
						if( it2 == d_b3_index.end() ){
							std::cout << d_vec[0] << " " << d_vec[1] << "  " << d_vec[2] << std::endl;
							d_b3_index.push_back(d_vec);
							// is_searched = true;
						}
					}

					// for(const auto& linked_site3 : site3->getLinkedSite() ){
						// auto it = find_if( linked_site3.second.begin(), linked_site3.second.end(),
						// 		[](const std::shared_ptr<Site> s){ return ( s->getSiteNum() == 0 );} );
						// if( it != linked_site3.second.end() ){
						// 	std::vector<double> d_vec = {linked_site.first, linked_site.first, linked_site3.first};
						// 	sort(d_vec.begin(), d_vec.end());
						// 	auto it2 = find( d_b3_index.begin(), d_b3_index.end(), d_vec );
						// 	if( it2 == d_b3_index.end() ){
						// 		std::cout << d_vec[0] << " " << d_vec[1] << "  " << d_vec[2] << std::endl;
						// 		d_b3_index.push_back(d_vec);
						// 		// is_searched = true;
						// 	}
						// }
					// }
				}
			}
		}
	}

	// for( int i=0; i<d_vec.size(); ++i ){
	// 	for( int j=i; j<d_vec.size(); ++j ){
	// 		for( int k=j; k<d_vec.size(); ++k ){
	// 			std::vector<double> tmp = {d_vec[i], d_vec[j], d_vec[k]};
	//
	// 			/*  search representative cluster */
	// 			bool is_searched = false;
	// 			for(const auto& l : site_vec[0]->getLinkedSiteviaDisncace(d_vec[i]) ){
	// 				for(const auto& m : l->getLinkedSiteviaDisncace(d_vec[j]) ){
	// 					if( m->getLinkedSite(d_vec[k], l->getSiteNum()) ) {
	// 						std::cout << d_vec[i] << " " << d_vec[j] << "  " << d_vec[k] << std::endl;
	// 						d_b3_index.push_back(tmp);
	// 						is_searched = true;
	// 						break;
	// 					}
	// 				}
	// 				if(is_searched) break;
	// 			}
	//
	// 		}
	// 	}
	// }

	std::vector<std::vector<std::shared_ptr<Site>>> clusters_3body;		/* site -> other 2 sites*/
	for(const auto& d_b3 : d_b3_index){
		for( const auto& site : site_vec ){
			for( const auto& i : site->getLinkedSiteviaDisncace(d_b3[0]) ){
				for( const auto& j : i->getLinkedSiteviaDisncace(d_b3[1]) ){
					if( j->getLinkedSite(d_b3[2], site->getSiteNum()) ) {
						clusters_out << site->getSiteNum() << " " << i->getSiteNum() << " " << j->getSiteNum() << " ";
						clusters_3body.push_back( std::vector<std::shared_ptr<Site>>{i ,j} );
					}
				}
			}
		}
		clusters_out << std::endl;
	}

	clusters_out.close();
	/*  set nbody = 4 */
	std::vector<std::vector<std::shared_ptr<Site>>> clusters_4body;		/* site -> other 2 sites*/
	// for(const auto& d_b3 : d_b3_index){

	// std::vector<std::vector<double>> d_b4_index;
	// for( int i=0; i<d_vec.size(); ++i ){
	// 	for( int j=i; j<d_vec.size(); ++j ){
	// 		for( int k=j; k<d_vec.size(); ++k ){
	// 			for( int l=k; l<d_vec.size(); ++l ){
	// 				std::vector<double> tmp = {d_vec[i], d_vec[j], d_vec[k], d_vec[l]};
	// 				d_b4_index.push_back(tmp);
	// 			}
	// 		}
	// 	}
	// }
	//
	// for(const auto& d_b4 : d_b4_index){
	// 	for( const auto& site : site_vec ){
	// 		for( const auto& i : site->getLinkedSiteviaDisncace(d_b4[0]) ){
	// 			for( const auto& j : i->getLinkedSiteviaDisncace(d_b4[1]) ){
	// 				for( const auto& k : j->getLinkedSiteviaDisncace(d_b4[2]) ){
	// 					auto l = k->getLinkedSite(d_b4[3], site->getSiteNum());
	// 					auto m = k->getLinkedSite(d_b4[3], i->getSiteNum());
	// 					auto n = j->getLinkedSite(d_b4[3], site->getSiteNum());
	// 					if( l and m and n ) {
	// 						std::cout << site->getSiteNum() << " " << i->getSiteNum() << " " << j->getSiteNum() << " " << k->getSiteNum() << " ";
	// 					}
	// 				}
	// 			}
	// 		}
	// 	}
	// 	std::cout << std::endl;
	// }

}
