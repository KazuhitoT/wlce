// read poscar.in
// command line argument -2=10 -3=11

// requiremnt Eigen boost
#include <unordered_map>
#include <map>
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

class hash_vecd {
public:
  double operator()(const std::vector<double> &x) const {
    const int C = 997;      // 素数
    double t = 0;
    for (int i = 0; i != x.size(); ++i) {
        t = t * C + x[i];
    }
    return t;
  }
};

class hash_vecd_less {
public:
	bool operator ()(const std::vector<double> &lhs, const std::vector<double> &rhs) const {
		return accumulate(lhs.begin(), lhs.end(), 0) < accumulate(rhs.begin(), rhs.end(), 0);
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

Eigen::Vector3d validDistance(const Eigen::Vector3d& lhs, const Eigen::Vector3d& rhs){
	Eigen::Vector3d d_vec = lhs - rhs;
	for( int i=0; i<3; ++i ){ /* boundary condition */
		if( d_vec(i) > (0.5-0.00001) )       {d_vec(i) = 1-d_vec(i);}
		else if( d_vec(i) < -(0.5+0.00001) ) {d_vec(i) = 1+d_vec(i);}
	}
	return d_vec;
};

bool is_line_cluster(std::vector<std::shared_ptr<Site>> cluster){
	Eigen::Vector3d center;
	center << 0,0,0;
	for( const auto& site : cluster ){
		center += site->getCoordinate();
	}
	center /= cluster.size();
	for( const auto& site : cluster ){
		if( center == site->getCoordinate() ) return true;
	}
	return false;
}


int main(int argc, char* argv[]){
	double d2, d3, d4 = 0;
	double prec = 0.00001;

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
			// for(int i=0; i<1; ++i){
				std::vector<int> vec_a;
				for(int j=0; j<position_ex.size(); ++j){
					Eigen::Vector3d d_vec = validDistance(position_ex[i], position_ex[j]);
					double distance = (lattice_ex * d_vec).norm();
					if( std::abs( distance - d_a.first ) < prec ) {
						site_vec[i]->setLinkedSite(d_a.first, site_vec[j]);
						// std::cout << i << " " << j << " ";
						vec_a.push_back(j);
					}
				}
				a_d_a.push_back(vec_a);
			}
			distance_site_to_sites[d_a.first] = a_d_a;
		}
	}

	for( const auto& i : distance_site_to_sites ){
		std::cout << i.first << " " << i.second[0].size() << std::endl;
		for( int j=0; j<i.second.size(); ++j) {
			assert( i.second[j].size() == i.second[0].size() );
			for( const auto& k : i.second[j] ){
				clusters_out << j << " " << k << " ";
			}
		}
		clusters_out << std::endl;
	}

	/*  set nbody = 3 */
	std::cout << " -- 3 body cluster" << std::endl;
	std::vector<std::vector<double>> d_b3_index;
	for(const auto& linked_site : site_vec[0]->getLinkedSite() ){
		for(const auto& site2 : linked_site.second ){
			for(const auto& linked_site2 : site2->getLinkedSite() ){
				for(const auto& site3 : linked_site2.second ){
					if( is_line_cluster( std::vector<std::shared_ptr<Site>>{site_vec[0], site2 ,site3} ) ){
						continue;
					}
					Eigen::Vector3d valid_ditance = validDistance(site3->getCoordinate(), site_vec[0]->getCoordinate());
					double distance = (lattice_ex * valid_ditance).norm();
					distance = site3->getDistance(distance);
					if( site3->getLinkedSite(distance, 0) ){
						std::vector<double> d_vec = {linked_site.first, linked_site2.first, distance};
						sort(d_vec.begin(), d_vec.end());
						auto it2 = find( d_b3_index.begin(), d_b3_index.end(), d_vec );
						if( it2 == d_b3_index.end() ){
							d_b3_index.push_back(d_vec);
						}
					}
				}
			}
		}
	}

	for(int i=2; i>=0; --i){
		stable_sort(d_b3_index.begin(), d_b3_index.end(), [i](const std::vector<double>& lhs, const std::vector<double>& rhs){ return lhs[i] < rhs[i]; });
	}

	/* index -> site -> other 2 sites*/
	std::unordered_map<std::vector<double>, std::vector<std::vector<std::vector<std::shared_ptr<Site>>>>, hash_vecd> clusters_3body;
	for(const auto& representative_d_b3 : d_b3_index){
		/* construct equivalent d_b3 vec */
		std::vector<std::vector<double>> all_d_b3;
		std::vector<double> data = representative_d_b3;
	  do{ all_d_b3.push_back(data); }while(next_permutation(data.begin(), data.end()));

		std::vector<std::vector<std::vector<std::shared_ptr<Site>>>> index_clusters;
		for( const auto& site : site_vec ){
			std::vector<std::vector<std::shared_ptr<Site>>> site_clusters;
			for(const auto& d_b3 : all_d_b3){
				for( const auto& i : site->getLinkedSiteviaDisncace(d_b3[0]) ){
					for( const auto& j : i->getLinkedSiteviaDisncace(d_b3[1]) ){
						if( j->getLinkedSite(d_b3[2], site->getSiteNum()) ) {
							std::vector<std::shared_ptr<Site>> tmp = {i ,j};
							sort(tmp.begin(), tmp.end());
							if( find(site_clusters.begin(), site_clusters.end(), tmp) == site_clusters.end() ){
								// !! 3*3*3 bug
								// 4 5.65685 5.65685
								// if( std::abs(d_b3[0]-4)<0.00001 and std::abs(d_b3[1]-5.65685 )<0.00001 and std::abs(d_b3[2]-5.65685)<0.00001 ){
								// 	std::cout << site->getSiteNum() << " " << tmp[0]->getSiteNum() <<" "<< tmp[1]->getSiteNum() << std::endl;
								// 	exit(1);
								// }
								site_clusters.push_back( tmp );
							}
						}
					}
				}
			}
			index_clusters.push_back(site_clusters);
		}
		clusters_3body[representative_d_b3] = index_clusters;
	}

	for(const auto& rep_d_b3 : d_b3_index){
		std::cout << rep_d_b3[0] << " " << rep_d_b3[1] << " " << rep_d_b3[2] << " " << clusters_3body[rep_d_b3][0].size() << std::endl;
		int i = 0;
		for(const auto& site_clusters : clusters_3body[rep_d_b3]){
			assert( clusters_3body[rep_d_b3][0].size() == site_clusters.size() );
			for(const auto& site_vec : site_clusters){
				for(const auto& site : site_vec){
					clusters_out << i << " " << site->getSiteNum() << " ";
				}
			}
			++i;
		}
		clusters_out << std::endl;
	}

	/*  set nbody = 4 */
	std::cout << " -- 4 body cluster" << std::endl;
	std::vector<std::vector<double>> d_b4_index;
	for(const auto& sites : clusters_3body){
		for(const auto& linked_site : sites.second[0] ){
			for(const auto& another_site : site_vec ){

				assert( linked_site.size() == 2 );
				if( another_site == site_vec[0] or another_site == linked_site[0] or another_site == linked_site[1])
					continue;

				Eigen::Vector3d valid_ditance1 = validDistance(another_site->getCoordinate(), site_vec[0]->getCoordinate());
				double distance1 = (lattice_ex * valid_ditance1).norm();
				distance1 = another_site->getDistance(distance1);
				Eigen::Vector3d valid_ditance2 = validDistance(another_site->getCoordinate(), linked_site[0]->getCoordinate());
				double distance2 = (lattice_ex * valid_ditance2).norm();
				distance2 = another_site->getDistance(distance2);
				Eigen::Vector3d valid_ditance3 = validDistance(another_site->getCoordinate(), linked_site[1]->getCoordinate());
				double distance3 = (lattice_ex * valid_ditance3).norm();
				distance3 = another_site->getDistance(distance3);

				if( another_site->getLinkedSite(distance1, 0) and another_site->getLinkedSite(distance2, linked_site[0]->getSiteNum()) and another_site->getLinkedSite(distance3, linked_site[1]->getSiteNum())){
					std::vector<double> d_vec = {sites.first[0], sites.first[1], sites.first[2], distance1, distance2, distance3};
					sort(d_vec.begin(), d_vec.end());
					// 2.82843 2.82843 2.82843 2.82843 2.82843 4
					// if( std::abs(d_vec[0]-2.82843)<0.00001 and std::abs(d_vec[1]-2.82843 )<0.00001 and std::abs(d_vec[2]-2.82843)<0.00001 and std::abs(d_vec[3]-2.82843)<0.00001 and std::abs(d_vec[4]-2.82843)<0.00001 and std::abs(d_vec[5]-4)<0.00001 ){
					// 	std::cout << 0 << " " << linked_site[0]->getSiteNum() <<" "<< linked_site[1]->getSiteNum() <<" "<< another_site->getSiteNum() << std::endl;
					// 	exit(1);
					// }
					auto it = find( d_b4_index.begin(), d_b4_index.end(), d_vec);
					if( it == d_b4_index.end() ){
						d_b4_index.push_back(d_vec);
					}
				}
			}
		}
	}

	for(int i=5; i>=0; --i){
		stable_sort(d_b4_index.begin(), d_b4_index.end(), [i](const std::vector<double>& lhs, const std::vector<double>& rhs){ return lhs[i] < rhs[i]; });
	}

	/* index -> site -> other 3 sites*/
	std::unordered_map<std::vector<double>, std::vector<std::vector<std::vector<std::shared_ptr<Site>>>>, hash_vecd> clusters_4body;
	for(const auto& representative_d_b4 : d_b4_index){
		/* construct equivalent d_b4 vec */
		std::vector<std::vector<double>> all_d_b4;
		std::vector<double> data = representative_d_b4;
	  do{ all_d_b4.push_back(data); }while(next_permutation(data.begin(), data.end()));

		std::vector<std::vector<std::vector<std::shared_ptr<Site>>>> index_clusters;
		for( const auto& site : site_vec ){
			std::vector<std::vector<std::shared_ptr<Site>>> site_clusters;
			for(const auto& d_b4 : all_d_b4){
				for( const auto& i : site->getLinkedSiteviaDisncace(d_b4[0]) ){
					for( const auto& j : i->getLinkedSiteviaDisncace(d_b4[1]) ){
						for( const auto& k : j->getLinkedSiteviaDisncace(d_b4[2]) ){
							if( site->getLinkedSite(d_b4[3], j->getSiteNum()) and  site->getLinkedSite(d_b4[4], k->getSiteNum()) and i->getLinkedSite(d_b4[5], k->getSiteNum()) ) {
								std::vector<std::shared_ptr<Site>> tmp = {i ,j, k};
								sort(tmp.begin(), tmp.end());
								if( find(site_clusters.begin(), site_clusters.end(), tmp) == site_clusters.end() ){
									// !! 3*3*3 bug
									// 4 5.65685 5.65685
									// if( std::abs(d_b3[0]-4)<0.00001 and std::abs(d_b3[1]-5.65685 )<0.00001 and std::abs(d_b3[2]-5.65685)<0.00001 ){
									// 	std::cout << site->getSiteNum() << " " << tmp[0]->getSiteNum() <<" "<< tmp[1]->getSiteNum() << std::endl;
									// 	exit(1);
									// }
									site_clusters.push_back( tmp );
								}
							}
						}
					}
				}
			}
			index_clusters.push_back(site_clusters);
		}
		clusters_4body[representative_d_b4] = index_clusters;
		break;
	}

	for(const auto& rep_d_b4 : d_b4_index){
		std::cout << rep_d_b4[0] << " " << rep_d_b4[1] << " " << rep_d_b4[2] << " " << rep_d_b4[3] << " " << rep_d_b4[4] << " " << rep_d_b4[5] << " " << clusters_4body[rep_d_b4][0].size() << std::endl;
		int i = 0;
		for(const auto& site_clusters : clusters_4body[rep_d_b4]){
			assert( clusters_4body[rep_d_b4][0].size() == site_clusters.size() );
			for(const auto& site_vec : site_clusters){
				clusters_out << i << " ";
				for(const auto& site : site_vec){
					clusters_out << site->getSiteNum() << " ";
				}
			}
			++i;
		}
		clusters_out << std::endl;
	}

	clusters_out.close();


		// 		Eigen::Vector3d valid_ditance = validDistance(site3->getCoordinate(), site_vec[0]->getCoordinate());
		// 		double distance = (lattice_ex * valid_ditance).norm();
		// 		distance = site3->getDistance(distance);
		// 		if( site3->getLinkedSite(distance, 0) ){
		// 			std::vector<double> d_vec = {linked_site.first, linked_site2.first, distance};
		// 			sort(d_vec.begin(), d_vec.end());
		// 			auto it2 = find( d_b3_index.begin(), d_b3_index.end(), d_vec );
		// 			if( it2 == d_b3_index.end() ){
		// 				d_b3_index.push_back(d_vec);
		// 			}
		// 		}
		// //
		// for(const auto& linked_site : site_vec[0]->getLinkedSite() ){
		// 	auto it = find(linked_site.second.begin(), linked_site.second.end(),
	// 		}
	// 	}
	// }


	//
	//
	// for(const auto& linked_site : site_vec[0]->getLinkedSite() ){
	// 	for(const auto& site2 : linked_site.second ){
	// 		for(const auto& linked_site2 : site2->getLinkedSite() ){
	// 			for(const auto& site3 : linked_site2.second ){
	// 				for(const auto& linked_site3 : site3->getLinkedSite() ){
	// 					for(const auto& site4 : linked_site3.second ){
	// 						for(const auto& linked_site4 : site4->getLinkedSite() ){
	// 							for(const auto& site5 : linked_site4.second ){
	// 								for(const auto& linked_site5 : site5->getLinkedSite() ){
	// 									for(const auto& site6 : linked_site5.second ){
	// 										Eigen::Vector3d valid_ditance = validDistance(site6->getCoordinate(), site_vec[0]->getCoordinate());
	// 										double distance = (lattice_ex * valid_ditance).norm();
	// 										if( site4->getLinkedSite(distance, 0) ){
	// 											std::vector<double> d_vec = {linked_site.first, linked_site2.first, linked_site3.first, linked_site4.first, linked_site5.first, distance};
	// 											sort(d_vec.begin(), d_vec.end());
	// 											auto it2 = find( d_b4_index.begin(), d_b4_index.end(), d_vec );
	// 											if( it2 == d_b4_index.end() ){
	// 												std::cout << d_vec[0] << " " << d_vec[1] << "  " << d_vec[2] << " " << d_vec[3] << " " << d_vec[4] << " " << d_vec[5] << std::endl;
	// 												d_b4_index.push_back(d_vec);
	// 												// is_searched = true;
	// 											}
	// 										}
	// 									}
	// 								}
	// 							}
	// 						}
	// 					}
	// 				}
	// 			}
	// 		}
	// 	}
	// }


	// std::vector<std::vector<std::shared_ptr<Site>>> clusters_4body;		/* site -> other 2 sites*/


}
