// read poscar.in
// command line argument -2=10 -3=11

// requiremnt Eigen boost
#include <unordered_map>
#include <memory>
#include <iterator>
#include <regex>
#include <Eigen/Core>
#include <Eigen/LU>
#include "./parsePoscar.hpp"
#include "./site.hpp"
#include "./conf2corr.hpp"

using allclusters = std::vector<std::vector<std::vector<std::vector<int>>>>;

class hash_vecd {
public:
  double operator()(const std::vector<double> &x) const {
    const int C = 997;
    double result = 0;
    for (int i=0; i<x.size(); ++i) {
        result = result * C + x[i];
    }
    return result;
  }
};

class hash_vec2d {
public:
  double operator()(const std::vector<std::vector<double>> &x) const {
    const int C = 997;
    double result = 0;
		for (int i=0; i<x.size(); ++i) {
			for (int j=0; j<x[i].size(); ++j) {
        result = result * C + x[i][j];
			}
    }
    return result;
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

int main(int argc, char* argv[]){
	double d2 = -1;
	double prec = 0.00001;

	std::regex re("(-d|-d2|-d3|-d4)\\=(\\d+(\\.\\d+)?)?");
	std::smatch match;
	for(int i=1; i<argc; i++) {
		std::string str(argv[i]);
		if(regex_match(str, match, re)){
			if( match[1] == "-d" ) d2=std::stod(match[2]);
		} else {
			std::cerr << " ERROR : invalid commandline argument [" << str << "]" << std::endl;
			exit(1);
		}
		std::cout << std::endl;
	}

	if( argc == 1 or d2==-1 ){
		std::cerr << " ERROR : commandline argument [-d=*] is required." << std::endl;
		exit(1);
	}

	ParsePoscar poscar("poscar.in");
	const int N = poscar.getAtoms().size();
	std::vector<double> spins_poscar = {-1, 1};

	double lattice[3][3];
	double lattice_unit[3][3];
	Eigen::Map<Eigen::Matrix3d>(&(lattice[0][0]), 3, 3)      = poscar.getLatticeBasis();
	Eigen::Map<Eigen::Matrix3d>(&(lattice_unit[0][0]), 3, 3) = poscar.getLatticeBasis();

	double position[N][3];
	double position_unit[N][3];
	const auto atoms = poscar.getAtoms();
	for( int i=0; i<N; ++i ){
		for( int j=0; j<3; ++j ){
			position[i][j]      = atoms[i].second[j];
			position_unit[i][j] = position[i][j];
		}
 	}

	std::vector<double> spins;
	for(int i=0; i<poscar.getAtomTypes().size(); ++i) {
		for(int j=0; j<poscar.getAtomTypes()[i]; ++j){
			spins.push_back(spins_poscar[i]);
		}
	}

	int N_unit = N;
	int expand_x = 1;
	int expand_y = 1;
	int expand_z = 1;
	Eigen::Matrix3d lattice_unit_matrix;
	lattice_unit_matrix = Eigen::Map<Eigen::Matrix3d>(&lattice[0][0]);

	std::vector<double> spins_ex;
	std::vector<Eigen::Vector3d> position_ex;
	Eigen::Matrix3d lattice_ex;
	double lattice_ex_arr[3][3];
	double position_ex_arr[N_unit*expand_x*expand_y*expand_z][3];
	{
		Eigen::Matrix3d div_mat;
		div_mat << 1./double(expand_x),0,0,
							 0,1./double(expand_y),0,
							 0,0,1./double(expand_z);
		/* !! Eigen::Map<Eigen::MatrixXd>(&(position_unit[0][0]), N_unit,3) does not work */
		std::vector<Eigen::Vector3d> position_unit_vec ;
		for(int i=0; i<N_unit; ++i){
			position_unit_vec.push_back( div_mat * Eigen::Map<Eigen::Vector3d>(&(position_unit[i][0])));
		}
		for(int x=0; x<expand_x; ++x){
			for(int y=0; y<expand_y; ++y){
				for(int z=0; z<expand_z; ++z){
					Eigen::Vector3d block_vec;
					block_vec << double(x)/double(expand_x), double(y)/double(expand_y), double(z)/double(expand_z);
					for(int i=0; i<N_unit; ++i){
						for(int ii=0; ii<3; ++ii){
							position_ex_arr[position_ex.size()][ii] = block_vec[ii] + position_unit_vec[i](ii);
						}
						position_ex.push_back( block_vec + position_unit_vec[i]);
					}
				}
			}
		}
		Eigen::Matrix3d max_ex_mat;
		max_ex_mat << expand_x,0,0, 0,expand_y,0, 0,0,expand_z;
		lattice_ex = max_ex_mat * lattice_unit_matrix;
		for(int i=0; i<3; ++i){
			for(int j=0; j<3; ++j){
				lattice_ex_arr[i][j] = lattice_ex(i,j);
			}
		}
	}
	outputPoscar(lattice_ex, position_ex, position_ex.size(), "expand");

	std::vector<std::shared_ptr<Site>> site_vec;
	for(int i=0; i<position_ex.size(); ++i) {
		site_vec.push_back( std::shared_ptr<Site>( new Site(i, position_ex[i])) );
	}
	for( const auto& site : site_vec){
		site->setRelativeSite(site_vec);
		assert( site->getSiteRelative().size() == site_vec.size() );
	}

	std::ofstream clusters_out( "clusters.out", std::ios::out );
	std::ofstream multiplicity_out( "multiplicity.out", std::ios::out );

	/*  set nbody = 1 */
	std::cout << " -- point cluster" << std::endl;
	std::cout << "0  0" << std::endl;
	std::vector<std::vector<std::vector<int>>> point_clusters;
	for( int i=0; i<position_ex.size(); ++i){
		clusters_out << i << " ";
		point_clusters.push_back(std::vector<std::vector<int>>());
	}
	clusters_out << std::endl;
	multiplicity_out << "1 " << position_ex.size() << std::endl;

	/*  set nbody = 2 */
	std::cout << " -- 2 body cluster" << std::endl;
	std::unordered_map<double, Eigen::Vector3d> distance_atom;
	for( int i=1; i<position_ex.size(); ++i ){
		Eigen::Vector3d d_vec = validCoordinate(position_ex[i], position_ex[0]);
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

	std::unordered_map<double, std::vector<std::vector<int>>> distance_site_to_sites;
	{
		for( const auto& d_a : distance_atom ){
			std::vector<std::vector<int>> a_d_a;
			for(int i=0; i<position_ex.size(); ++i){
			// for(int i=0; i<1; ++i){
				std::vector<int> vec_a;
				for(int j=0; j<position_ex.size(); ++j){
					Eigen::Vector3d d_vec = validCoordinate(position_ex[i], position_ex[j]);
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
		multiplicity_out << "2 " << i.second[0].size()*position_ex.size()/2  << std::endl;
		for( int j=0; j<i.second.size(); ++j) {
			std::vector<std::vector<int>> site_clusters;
			assert( i.second[j].size() == i.second[0].size() );
			for( const auto& k : i.second[j] ){
				clusters_out << j << " " << k << " ";
			}
		}
		clusters_out << std::endl;
	}


	/*  set nbody = 3 */
	std::cout << " -- 3 body cluster" << std::endl;
	/*
	 *	triplet_cluster;
	 *	key: All lines of the 3body cluster
	 *	val: all relative sites for site_vec[0]
	 */
	std::unordered_map<std::vector<double>, std::vector<std::vector<Eigen::Vector3d>>, hash_vecd> triplet_cluster;
	std::vector<std::vector<double>> d_b3_index;
	for(const auto& linked_site : site_vec[0]->getLinkedSite() ){
		for(const auto& site2 : linked_site.second ){
			for(const auto& linked_site2 : site2->getLinkedSite() ){
				for(const auto& site3 : linked_site2.second ){

					/* reject line cluster */
					// Eigen::Matrix3d triplet;
					// auto vector_02 = validCoordinate(site_vec[0]->getCoordinate(), site2->getCoordinate());
					// auto vector_03 = validCoordinate(site_vec[0]->getCoordinate(), site3->getCoordinate());
					// triplet << vector_02[0], vector_02[1], vector_02[2],
					// 	 vector_03[0], vector_03[1], vector_03[2],
					// 		1, 1, 1;
					// double triplet_area  = triplet.determinant();
					// if( triplet_area < prec ) {
					// 	continue;
					// }

					Eigen::Vector3d valid_ditance = validCoordinate(site3->getCoordinate(), site_vec[0]->getCoordinate());
					double distance = (lattice_ex * valid_ditance).norm();
					distance = site3->getDistance(distance);
					if( distance > 0 and site3->getLinkedSite(distance, 0) ){
						std::vector<double> d_vec = {linked_site.first, linked_site2.first, distance};
						std::vector<Eigen::Vector3d> relative_coords = {
							validCoordinate(site2->getCoordinate(), site_vec[0]->getCoordinate()),
							validCoordinate(site3->getCoordinate(), site_vec[0]->getCoordinate())
						};
						for(int i=2; i>=0; --i){
							stable_sort(relative_coords.begin(), relative_coords.end(), [i](const Eigen::Vector3d& lhs, const Eigen::Vector3d& rhs){ return lhs[i] < rhs[i]; });
						}
						sort(d_vec.begin(), d_vec.end());
						auto it = find( d_b3_index.begin(), d_b3_index.end(), d_vec );
						if( it == d_b3_index.end() ){
							d_b3_index.push_back(d_vec);
							triplet_cluster[d_vec].push_back(relative_coords);
						} else {
							auto it2 = find_if( triplet_cluster[d_vec].begin(), triplet_cluster[d_vec].end(),
							[relative_coords, prec](const std::vector<Eigen::Vector3d>& obj){
								for(int i=0; i<2; ++i){
									if( (relative_coords[i]-obj[i]).norm() > prec ) return false;
								}
								return true;;
							});
							if( it2 == triplet_cluster[d_vec].end() ){
								triplet_cluster[d_vec].push_back(relative_coords);
							}
						}
					}
				}
			}
		}
	}

	for(int i=2; i>=0; --i){
		stable_sort(d_b3_index.begin(), d_b3_index.end(), [i](const std::vector<double>& lhs, const std::vector<double>& rhs){ return lhs[i] < rhs[i]; });
	}

	for(const auto& index : d_b3_index){
		std::cout << "multiplicity = " << triplet_cluster[index].size() << std::endl;
		multiplicity_out << "3 " << triplet_cluster[index].size()*position_ex.size()/3 << std::endl;
		std::cout << index[0] << " " << index[1] << " " << index[2] << std::endl;
		std::cout << " -- relative_coords -- "<< std::endl;
		std::cout << site_vec[0]->getCoordinate().transpose() << std::endl;
		std::cout << triplet_cluster[index][0][0].transpose() << std::endl;
		std::cout << triplet_cluster[index][0][1].transpose() << std::endl;
		std::cout << std::endl;
	}

	for(const auto& index : d_b3_index){
		for(const auto& site : site_vec ){
			std::vector<std::vector<int>> site_clusters;
			for(const auto& two_relative_coord : triplet_cluster[index]){
				assert( two_relative_coord.size() == 2 );
				assert( site->getSiteRelative().size() == position_ex.size() );
				assert( site->getSiteRelative().size() == site_vec.size() );
				auto site2 = site->getRelativeSite(two_relative_coord[0]);
				auto site3 = site->getRelativeSite(two_relative_coord[1]);
				if( site2 and site3 ) {
					clusters_out << site->getSiteNum() << " ";
					clusters_out << site2->getSiteNum() << " ";
					clusters_out << site3->getSiteNum() << " ";
				}
			}
		}
		clusters_out << std::endl;
	}

	// /*  set nbody = 4 */
	std::cout << " -- 4 body cluster" << std::endl;
	std::vector<std::vector<std::vector<double>>> d_b4_index;
	std::unordered_map<std::vector<std::vector<double>>, std::vector<std::vector<Eigen::Vector3d>>, hash_vec2d> quadlet_cluster;
	for( const auto& index : d_b3_index ){
		for( int i=0; i<triplet_cluster[index].size(); ++i){
			for( const auto& linked_sites : site_vec[0]->getLinkedSite() ){
				for( const auto& linked_site : linked_sites.second ){

					if(  validCoordinate(linked_site->getCoordinate(), triplet_cluster[index][i][0]).norm() < prec
						or validCoordinate(linked_site->getCoordinate(), triplet_cluster[index][i][1]).norm() < prec ){
							continue;
					}

					auto all_triplets = getAllTriplets( std::vector<Eigen::Vector3d>{
						site_vec[0]->getCoordinate(),
						triplet_cluster[index][i][0],
						triplet_cluster[index][i][1],
						linked_site->getCoordinate()
					}, lattice_ex , d2);
					if( all_triplets.size() != 4 ) continue;

					std::vector<Eigen::Vector3d> relative_coords = {
						validCoordinate(triplet_cluster[index][i][0], site_vec[0]->getCoordinate()),
						validCoordinate(triplet_cluster[index][i][1], site_vec[0]->getCoordinate()),
						validCoordinate(linked_site->getCoordinate(), site_vec[0]->getCoordinate())
					};

					for(int i=2; i>=0; --i){
						stable_sort(relative_coords.begin(), relative_coords.end(), [i](const Eigen::Vector3d& lhs, const Eigen::Vector3d& rhs){ return lhs[i] < rhs[i]; });
					}

					// auto it = find(d_b4_index.begin(), d_b4_index.end(), all_triplets);
					auto it = find_if( d_b4_index.begin(), d_b4_index.end(),
						[all_triplets, prec](const std::vector<std::vector<double>>& obj){
							for(int i=0; i<obj.size(); ++i){
								double norm = 0;
								for(int j=0; j<obj[i].size(); ++j){
									norm += std::pow(all_triplets[i][j]-obj[i][j], 2.);
								}
								if( std::sqrt(norm) > 0.001 ){
									return false;
								}
							}
							return true;
						});
					if( it==d_b4_index.end() ) {
						d_b4_index.push_back(all_triplets);
						quadlet_cluster[all_triplets].push_back(relative_coords);
					} else {
						auto it2 = find_if( quadlet_cluster[(*it)].begin(), quadlet_cluster[(*it)].end(),
						[relative_coords, prec](const std::vector<Eigen::Vector3d>& obj){
							for(int i=0; i<3; ++i){
								if( (relative_coords[i]-obj[i]).norm() > prec ) return false;
							}
							return true;
						});
						if( it2 == quadlet_cluster[(*it)].end() ){
							quadlet_cluster[(*it)].push_back(relative_coords);
						}
					}
				}
			}
		}
	}
	//
	for(int i=3; i>=0; --i){
		for(int j=2; j>=0; --j){
		stable_sort(d_b4_index.begin(), d_b4_index.end(),
			[i,j](const std::vector<std::vector<double>>& lhs, const std::vector<std::vector<double>>& rhs){
				return lhs[i][j] < rhs[i][j];
			});
		}
	}

	for(const auto& index : d_b4_index){
		std::cout << "multiplicity = " << quadlet_cluster[index].size() << std::endl;
		multiplicity_out << "4 " << quadlet_cluster[index].size()*position_ex.size()/4 << std::endl;
		std::cout << " -- relative_coords -- "<< std::endl;
		std::cout << site_vec[0]->getCoordinate().transpose() << std::endl;
		std::cout << quadlet_cluster[index][0][0].transpose() << std::endl;
		std::cout << quadlet_cluster[index][0][1].transpose() << std::endl;
		std::cout << quadlet_cluster[index][0][2].transpose() << std::endl;
		std::cout << std::endl;

		std::cout << " -- all triangles  -- "<< std::endl;
		for(const auto& i : index){
			std::cout << "[";
			for(const auto& j : i){
				std::cout << j << ",";
			}
			std::cout << "]" << std::endl;
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}

	for(const auto& index : d_b4_index){
		for(const auto& site : site_vec ){
			for(const auto& three_relative_coord : quadlet_cluster[index]){
				auto site2 = site->getRelativeSite(three_relative_coord[0]);
				auto site3 = site->getRelativeSite(three_relative_coord[1]);
				auto site4 = site->getRelativeSite(three_relative_coord[2]);
				if( site2 and site3 and site4 ) {
					clusters_out << site->getSiteNum() << " ";
					clusters_out << site2->getSiteNum() << " ";
					clusters_out << site3->getSiteNum() << " ";
					clusters_out << site4->getSiteNum() << " ";
				}
			}
		}
		clusters_out << std::endl;
	}

	clusters_out.close();
	multiplicity_out.close();

	// Conf2corr("poscar.in", std::vector<double>{-1, 1}, std::vector<double>{-1, 1}, pall_clusters);

}
