#ifndef __SEARCH_EQUIV_SITES_HPP
#define __SEARCH_EQUIV_SITES_HPP

#include <unordered_map>
#include "./myhash.hpp"
#include "./myindex.hpp"
#include "./site.hpp"
#include "./lattice.hpp"
#include "./parser.hpp"

void setSiteAndLattice(std::vector<std::shared_ptr<Lattice>>&, std::vector<std::shared_ptr<Site>>&, const ParsePoscar& poscar);
void setSiteReferencesAndPairIndex(
		const std::vector<std::shared_ptr<Lattice>>& vec_lattices,
		std::vector<std::shared_ptr<Site>>& vec_sites,
		std::unordered_map<PairIndex , Eigen::Vector3d, hash_pairindex>& index_pair_clusters,
		double pair_truncation,
		double precision
	);
void setTripletIndex(
		const std::vector<std::shared_ptr<Lattice>>& vec_lattices,
		const std::vector<std::shared_ptr<Site>>& vec_sites,
		const std::vector<PairIndex>& sorted_vec_pairindex,
		std::vector<TripletIndex>& keys_index_triplet_cluster,
		std::unordered_map<TripletIndex, std::vector<std::vector<Eigen::Vector3d>>, hash_tripletindex>& index_triplet_cluster,
		double triplet_truncation,
		double precision
	);
void setQuadrupletIndex(
		const std::vector<std::shared_ptr<Lattice>>& vec_lattices,
		const std::vector<std::shared_ptr<Site>>& vec_sites,
		const std::vector<PairIndex>& keys_index_pair_cluster,
		const std::unordered_map<TripletIndex, std::vector<std::vector<Eigen::Vector3d>>, hash_tripletindex>& index_triplet_cluster,
		std::vector<QuadrupletIndex>& keys_index_quadruplet_cluster,
		std::unordered_map<QuadrupletIndex, std::vector<std::vector<Eigen::Vector3d>>, hash_quadrupletindex>& index_quadruplet_cluster,
		double quadruplet_truncation,
		double precision
	);

#endif
