
#include "./search_equiv_sites.hpp"

void setSiteAndLattice(std::vector<std::shared_ptr<Lattice>>& vec_lattices, std::vector<std::shared_ptr<Site>>& vec_sites, const ParsePoscar& poscar){
	const auto atoms = poscar.getAtoms();
	const std::vector<int> vec_atom_types = poscar.getAtomTypes();
	for(int i=0, j=0; i<vec_atom_types.size(); ++i) {
		auto lattice = std::shared_ptr<Lattice>( new Lattice(i) );
		lattice->setBasis(poscar.getLatticeBasis());
		vec_lattices.push_back(lattice);
		for(int k=0; k<vec_atom_types[i]; ++j, ++k){
			auto site = std::shared_ptr<Site>( new Site(atoms[j].first, atoms[j].second, lattice) );
			vec_sites.push_back(site);
			lattice->setSite(site);
		}
	}
}


void setSiteReferencesAndPairIndex(
		const std::vector<std::shared_ptr<Lattice>>& vec_lattices,
		std::vector<std::shared_ptr<Site>>& vec_sites,
		std::unordered_map<PairIndex , Eigen::Vector3d, hash_pairindex>& index_pair_clusters,
		double pair_truncation,
		double precision
	){

	std::vector<double> vec_pair_distance;  /*  hash function cannot correctly handle float value, so we should cast a specific float value  */

	for(const auto lattice : vec_lattices){
		for(const auto base_site : lattice->getSites()){
			const auto lattice_basis = lattice->getLatticeBasis();
			for( const auto site : vec_sites ){

				const Eigen::Vector3d relative_coord = validCoordinate(site->getCoordinate(), base_site->getCoordinate());

				double distance = (lattice_basis * relative_coord).norm();
				if ( precision < distance and distance <= pair_truncation ) {

					auto itr_distance = std::find_if( vec_pair_distance.begin(), vec_pair_distance.end(), [precision, distance](const double x){
						return std::fabs(x-distance) < precision;
					});
					if( itr_distance != vec_pair_distance.end() ) distance = (*itr_distance);
					else vec_pair_distance.push_back(distance);

					PairIndex pair_index(base_site, site, distance);

					auto it = index_pair_clusters.find(pair_index);
					if( it == index_pair_clusters.end() ){
						index_pair_clusters[pair_index] = relative_coord;
					}

					base_site->setLinkedSite(site, distance, relative_coord);

				}

			}
		}
	}

}

void setTripletIndex(
		const std::vector<std::shared_ptr<Lattice>>& vec_lattices,
		const std::vector<std::shared_ptr<Site>>& vec_sites,
		const std::vector<PairIndex>& sorted_vec_pairindex,
		std::vector<TripletIndex>& keys_index_triplet_cluster,
		std::unordered_map<TripletIndex, std::vector<std::vector<Eigen::Vector3d>>, hash_tripletindex>& index_triplet_cluster,
		double triplet_truncation,
		double precision
	){

	for(const auto& pairindex : sorted_vec_pairindex ){
		auto pair_itr = vec_sites[0]->getLinkedSiteIterator(pairindex.distance);
		while( pair_itr.first != pair_itr.second ){
			for(const auto& pairindex2 : sorted_vec_pairindex ){
			auto pair_itr2 = (*(pair_itr.first)).site->getLinkedSiteIterator(pairindex2.distance);
			while( pair_itr2.first != pair_itr2.second ){

				auto triplet_index = getAllTripletIndex( std::vector<std::shared_ptr<Site>>{
					vec_sites[0],
					(*(pair_itr.first)).site,
					(*(pair_itr2.first)).site
				}, vec_lattices[0]->getLatticeBasis(), triplet_truncation, sorted_vec_pairindex);

				if( triplet_index.size() != 1 ) {
					++pair_itr2.first;
					continue;
				}

				std::vector<Eigen::Vector3d> relative_coords = {
					validCoordinate((*(pair_itr.first)).site->getCoordinate(), vec_sites[0]->getCoordinate()),
					validCoordinate((*(pair_itr2.first)).site->getCoordinate(), vec_sites[0]->getCoordinate())
				};

				for(int i=2; i>=0; --i){
					stable_sort(relative_coords.begin(), relative_coords.end(), [i](const Eigen::Vector3d& lhs, const Eigen::Vector3d& rhs){ return lhs[i] < rhs[i]; });
				}

				auto it = find( keys_index_triplet_cluster.begin(), keys_index_triplet_cluster.end(), triplet_index[0] );
				if( it == keys_index_triplet_cluster.end() ){
					keys_index_triplet_cluster.push_back(triplet_index[0]);
					index_triplet_cluster[triplet_index[0]].push_back(relative_coords);
				} else {
					auto it2 = find_if( index_triplet_cluster[triplet_index[0]].begin(), index_triplet_cluster[triplet_index[0]].end(),
					[&relative_coords, precision](const std::vector<Eigen::Vector3d>& obj){
						for(int i=0; i<2; ++i){
							if( (relative_coords[i]-obj[i]).norm() > precision ) return false;
						}
						return true;
					});
					if( it2 == index_triplet_cluster[triplet_index[0]].end() ){
						index_triplet_cluster[triplet_index[0]].push_back(relative_coords);
					}
				}

				++pair_itr2.first;
				}
			}
			++pair_itr.first;
		}
	}

}


void setQuadrupletIndex(
		const std::vector<std::shared_ptr<Lattice>>& vec_lattices,
		const std::vector<std::shared_ptr<Site>>& vec_sites,
		const std::vector<PairIndex>& keys_index_pair_cluster,
		const std::unordered_map<TripletIndex, std::vector<std::vector<Eigen::Vector3d>>, hash_tripletindex>& map_index_triplet_cluster,
		std::vector<QuadrupletIndex>& keys_index_quadruplet_cluster,
		std::unordered_map<QuadrupletIndex, std::vector<std::vector<Eigen::Vector3d>>, hash_quadrupletindex>& index_quadruplet_cluster,
		double quadruplet_truncation,
		double precision
	){

	for( const auto& pair_index_triplet_cluster :  map_index_triplet_cluster ){
		for( const auto pairindex : keys_index_pair_cluster){
			if( pairindex.distance > quadruplet_truncation ) continue;

			for( const auto& vec_index_triplet_cluster :  pair_index_triplet_cluster.second ){

				auto pair_itr = vec_sites[0]->getLinkedSiteIterator(pairindex.distance);
				while( pair_itr.first != pair_itr.second ){

					// const auto b = index_triplet_cluster.second[0];
					if(  validCoordinate((*(pair_itr.first)).site->getCoordinate(), vec_index_triplet_cluster[0]).norm() < precision
						or validCoordinate((*(pair_itr.first)).site->getCoordinate(), vec_index_triplet_cluster[1]).norm() < precision ){
							++pair_itr.first;
							continue;
					}

					auto all_triplets = getAllTripletIndex( std::vector<std::shared_ptr<Site>>{
						vec_sites[0],
						vec_sites[0]->getLinkedSite(vec_index_triplet_cluster[0]),
						vec_sites[0]->getLinkedSite(vec_index_triplet_cluster[1]),
						(*(pair_itr.first)).site
					}, vec_lattices[0]->getLatticeBasis(), quadruplet_truncation, keys_index_pair_cluster);
					if( all_triplets.size() != 4 ) {
						++pair_itr.first;
						continue;
					}

					for(int i=2; i>=0; --i){
						stable_sort(all_triplets.begin(), all_triplets.end(),
							[i](const TripletIndex& lhs, const TripletIndex& rhs){ return lhs.vec_pairindex[i] < rhs.vec_pairindex[i]; });
					}

					QuadrupletIndex quadrupletindex = QuadrupletIndex(all_triplets);

					std::vector<Eigen::Vector3d> relative_coords = {
						vec_index_triplet_cluster[0],
						vec_index_triplet_cluster[1],
						validCoordinate((*(pair_itr.first)).site->getCoordinate(), vec_sites[0]->getCoordinate())
					};

					for(int i=2; i>=0; --i){
						stable_sort(relative_coords.begin(), relative_coords.end(), [i](const Eigen::Vector3d& lhs, const Eigen::Vector3d& rhs){ return lhs[i] < rhs[i]; });
					}

					auto it = find( keys_index_quadruplet_cluster.begin(), keys_index_quadruplet_cluster.end(), quadrupletindex );
					if( it == keys_index_quadruplet_cluster.end() ){
						keys_index_quadruplet_cluster.push_back(quadrupletindex);
						index_quadruplet_cluster[quadrupletindex].push_back(relative_coords);
					} else {
						auto it2 = find_if( index_quadruplet_cluster[quadrupletindex].begin(), index_quadruplet_cluster[quadrupletindex].end(),
						[&relative_coords, precision](const std::vector<Eigen::Vector3d>& obj){
							for(int i=0; i<3; ++i){
								if( (relative_coords[i]-obj[i]).norm() > precision ) return false;
							}
							return true;
						});
						if( it2 == index_quadruplet_cluster[quadrupletindex].end() ){
							index_quadruplet_cluster[quadrupletindex].push_back(relative_coords);
						}
					}
					++pair_itr.first;
				}
			}
		}
	}

}
