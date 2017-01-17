#include <unordered_map>
#include <memory>
#include <iterator>
#include <Eigen/Core>
#include <Eigen/LU>
// #include <boost/program_options.hpp>
#include "./parser.hpp"
#include "./lattice.hpp"
#include "./site.hpp"
#include "./search_equiv_sites.hpp"
#include "./myhash.hpp"
#include "./myindex.hpp"
#include "./input.hpp"

using allclusters = std::vector<std::vector<std::vector<std::vector<int>>>>;

void dispInput(std::shared_ptr<Input> in){
	//  [INPUT]
	std::cout << "MAXBODY     : " << in->getDataByString("MAXBODY") << std::endl;
	std::cout << "TRUNCATION  : " << in->getDataByString("TRUNCATION") << std::endl;

	//  [OPTION]
	std::cout << "PAIRTRUNC   : " << in->getDataByString("PAIRTRUNC")  << std::endl;
	std::cout << "TRITRUNC    : " << in->getDataByString("TRITRUNC") << std::endl;
	std::cout << "QUADTRUNC   : " << in->getDataByString("QUADTRUNC") << std::endl;
}

int main(int argc, char* argv[]){
	std::string filename_poscar_in = "poscar.in";

	int num_max_body = 0;
	double truncation = 0;
	double precision  = 0.00001;
	double pair_truncation = 0;
	double triplet_truncation = 0;
	double quadruplet_truncation = 0;
	std::shared_ptr<Input> in(new Input("cluster.ini"));

	in->setData("MAXBODY",      num_max_body, true);
	in->setData("TRUNCATION",   truncation, true);

	in->setData("PRECISION", precision);
	in->setData("PAIRTRUNC", pair_truncation);
	in->setData("TRITRUNC",  triplet_truncation);
	in->setData("QUADTRUNC", quadruplet_truncation);

	if( pair_truncation == 0 ) pair_truncation = truncation;
	if( triplet_truncation == 0 ) triplet_truncation = truncation;
	if( quadruplet_truncation == 0 ) quadruplet_truncation = truncation;

	if( pair_truncation < triplet_truncation or triplet_truncation < quadruplet_truncation ) {
		std::cerr << "ERROR : TRITRUNC <= PAIRTRUNC(orTRUNATION) and QUADTRUNC <= TRITRUNC(orTRUNATION)" << std::endl;
		exit(1);
	}

	ParsePoscar poscar(filename_poscar_in.c_str());
	const auto atoms_unit = poscar.getAtoms();

	Eigen::Vector3d unit_x, unit_y, unit_z;
	unit_x << 1, 0, 0;
	unit_y << 0, 1, 0;
	unit_z << 0, 0, 1;

	double unit_length_x = (poscar.getLatticeBasis() * unit_x).norm();
	double unit_length_y = (poscar.getLatticeBasis() * unit_y).norm();
	double unit_length_z = (poscar.getLatticeBasis() * unit_z).norm();

	int N_unit = poscar.getAtoms().size();
	int expand_x = std::ceil((pair_truncation*2.)/unit_length_x)+2;
	int expand_y = std::ceil((pair_truncation*2.)/unit_length_y)+2;
	int expand_z = std::ceil((pair_truncation*2.)/unit_length_z)+2;

	poscar.expandPoscar(expand_x, expand_y, expand_z);

	const auto atoms = poscar.getAtoms();
	const int N = atoms.size();

	/*  setting lattice and site */
	std::vector<std::shared_ptr<Lattice>> vec_lattices;
	std::vector<std::shared_ptr<Site>> vec_sites;
	setSiteAndLattice(vec_lattices, vec_sites, poscar);

	std::ofstream cluster_info("cluster.info", std::ios::out );
	std::ofstream cluster_out( "cluster.out",  std::ios::out );

	for( int i=0; i<atoms_unit.size(); ++i){
		cluster_out << atoms_unit[i].first << " " <<  atoms_unit[i].second.transpose() << std::endl;
	}
	cluster_out << "--" << std::endl;

	/*  set nbody = 1 */
	cluster_info << " -- point cluster" << std::endl;
	int num_index = 1;
	for( const auto& lattice : vec_lattices ){
		cluster_info << num_index << " : "  << "[" << lattice->getLatticeNum() << "] : 0 : 1 : 0 0 0 " << std::endl;
		cluster_out << num_index << " ";
		cluster_out << "1" << " ";
		cluster_out << N_unit << " ";
		cluster_out << N_unit << " ";
		for( int i=0; i<N_unit; ++i){
			cluster_out << vec_sites[i]->getSiteNum() << " ";
		}
		cluster_out << std::endl;
		++num_index;
	}

	/*  set nbody = 2 */
	cluster_info << " -- 2 body cluster" << std::endl;
	std::unordered_map<PairIndex , Eigen::Vector3d, hash_pairindex> index_pair_clusters;
	if( pair_truncation == 0 ) pair_truncation = truncation;
	setSiteReferencesAndPairIndex(vec_lattices, vec_sites, index_pair_clusters, pair_truncation, precision);

	std::vector<PairIndex> sorted_vec_pairindex;
	for(const auto& index_pair : index_pair_clusters){
		sorted_vec_pairindex.push_back(index_pair.first);
	}
	sort(sorted_vec_pairindex.begin(), sorted_vec_pairindex.end());

	cluster_info << " # [index] : distance : multiplicity : coordination " << std::endl;
	for( const auto& pairindex : sorted_vec_pairindex ){
		cluster_info << num_index << " : ";
		cluster_info << "[" << pairindex.lattice_index[0] << "," << pairindex.lattice_index[1]  << "] : ";
		const auto pair_itr = vec_sites[0]->getLinkedSiteIterator(pairindex.distance);
		cluster_info << pairindex.distance << " : " << std::distance(pair_itr.first, pair_itr.second) << " : ";
		cluster_info << index_pair_clusters[pairindex].transpose() << std::endl;

		for( const auto& site : vec_sites ) {
			const auto test_pair_itr = site->getLinkedSiteIterator(pairindex.distance);
			assert( std::distance(pair_itr.first, pair_itr.second) == std::distance(test_pair_itr.first, test_pair_itr.second) );
		}

		cluster_out << num_index << " ";
		cluster_out << "2" << " ";
		cluster_out << N_unit << " ";
		cluster_out << N_unit * std::distance(pair_itr.first, pair_itr.second) << " ";
		for(int i=0; i<N_unit; ++i){
			auto pair_itr = vec_sites[i]->getLinkedSiteIterator(pairindex.distance);
			while( pair_itr.first != pair_itr.second ){
				cluster_out << vec_sites[i]->getSiteNum() << " " << (*(pair_itr.first)).site->getSiteNum() << " ";
				++pair_itr.first;
			}
		}
		cluster_out << std::endl;
		++num_index;
	}
	if( num_max_body == 2 ){
		cluster_out.close();
		cluster_info.close();
		return 1;
	}

	/*  set nbody = 3 */
	cluster_info<< " -- 3 body cluster" << std::endl;
	std::vector<TripletIndex> keys_index_triplet_cluster;
	std::unordered_map<TripletIndex, std::vector<std::vector<Eigen::Vector3d>>, hash_tripletindex> index_triplet_cluster;
	setTripletIndex(vec_lattices, vec_sites, sorted_vec_pairindex, keys_index_triplet_cluster, index_triplet_cluster, triplet_truncation, precision);

	for(int i=2; i>=0; --i){
		stable_sort(keys_index_triplet_cluster.begin(), keys_index_triplet_cluster.end(), [i](const TripletIndex& lhs, const TripletIndex& rhs){
				return lhs.vec_pairindex[i] < rhs.vec_pairindex[i];
		});
	}

	for(const auto& tripletindex : keys_index_triplet_cluster){

		cluster_info << num_index << " : ";
		cluster_info << "[] : ";
		cluster_info << " [ ";
		cluster_info << tripletindex.vec_pairindex[0].distance << ",";
		cluster_info << tripletindex.vec_pairindex[1].distance << ",";
		cluster_info << tripletindex.vec_pairindex[2].distance;
		cluster_info << " ] : ";
		cluster_info << index_triplet_cluster[tripletindex].size() <<  " : ";

		cluster_info << index_triplet_cluster[tripletindex][0][0].transpose() << ",";
		cluster_info << index_triplet_cluster[tripletindex][0][1].transpose() << std::endl;

		cluster_out << num_index << " ";
		cluster_out << "3" << " ";
		cluster_out << N_unit << " ";
		cluster_out << N_unit*index_triplet_cluster[tripletindex].size() << " ";
		for(int i=0; i<N_unit; ++i){
			for(int j=0; j<index_triplet_cluster[tripletindex].size(); ++j){
				cluster_out  << vec_sites[i]->getSiteNum() << " ";
				cluster_out  << vec_sites[i]->getLinkedSite(index_triplet_cluster[tripletindex][j][0])->getSiteNum() << " ";
				cluster_out  << vec_sites[i]->getLinkedSite(index_triplet_cluster[tripletindex][j][1])->getSiteNum() << " ";
			}
		}
		cluster_out << std::endl;
		++num_index;
	}
	if( num_max_body == 3 ){
		cluster_out.close();
		cluster_info.close();
		return 1;
	}

	// /*  set nbody = 4 */
	cluster_info << " -- 4 body cluster" << std::endl;
	std::vector<QuadrupletIndex> keys_index_quadruplet_cluster;
	std::unordered_map<QuadrupletIndex, std::vector<std::vector<Eigen::Vector3d>>, hash_quadrupletindex> index_quadruplet_cluster;
	setQuadrupletIndex(vec_lattices, vec_sites, sorted_vec_pairindex,  index_triplet_cluster, keys_index_quadruplet_cluster, index_quadruplet_cluster, quadruplet_truncation, precision);

	for(int i=3; i>=0; --i){
		for(int j=2; j>=0; --j){
		stable_sort(keys_index_quadruplet_cluster.begin(), keys_index_quadruplet_cluster.end(),
			[i,j](const QuadrupletIndex& lhs, const QuadrupletIndex& rhs){
				return lhs.vec_tripletindex[i].vec_pairindex[j] < rhs.vec_tripletindex[i].vec_pairindex[j];
			});
		}
	}

	for(const auto& key_index_quadruplet_cluster : keys_index_quadruplet_cluster){

		cluster_info << num_index << " : ";
		cluster_info << "[] : ";
		cluster_info << " [ ";
		for( const auto distance : key_index_quadruplet_cluster.getVecDistance() ){
			cluster_info << distance << ",";
		}
		cluster_info << " ] : ";
		cluster_info << index_quadruplet_cluster[key_index_quadruplet_cluster].size() <<  " : ";

		cluster_info << index_quadruplet_cluster[key_index_quadruplet_cluster][0][0].transpose() << ",";
		cluster_info << index_quadruplet_cluster[key_index_quadruplet_cluster][0][1].transpose() << ",";
		cluster_info << index_quadruplet_cluster[key_index_quadruplet_cluster][0][2].transpose() << std::endl;

		cluster_out << num_index << " ";
		cluster_out << "4" << " ";
		cluster_out << N_unit << " ";
		cluster_out << N_unit*index_quadruplet_cluster[key_index_quadruplet_cluster].size() << " ";
		for(int i=0; i<N_unit; ++i){
			for(int j=0; j<index_quadruplet_cluster[key_index_quadruplet_cluster].size(); ++j){
				cluster_out  << vec_sites[i]->getSiteNum() << " ";
				cluster_out  << vec_sites[i]->getLinkedSite(index_quadruplet_cluster[key_index_quadruplet_cluster][j][0])->getSiteNum() << " ";
				cluster_out  << vec_sites[i]->getLinkedSite(index_quadruplet_cluster[key_index_quadruplet_cluster][j][1])->getSiteNum() << " ";
				cluster_out  << vec_sites[i]->getLinkedSite(index_quadruplet_cluster[key_index_quadruplet_cluster][j][2])->getSiteNum() << " ";
			}
		}
		cluster_out << std::endl;
		++num_index;
	}
	cluster_out.close();
	cluster_info.close();

}
