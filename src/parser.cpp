#include <algorithm>
#include <tuple>
#include "./parser.hpp"

using allclusters = std::vector<std::vector<std::vector<std::vector<int>>>>;

ParseClusterOut::ParseClusterOut (const char* f, std::vector<int> vec_index_read):filename(f), pall_clusters(new allclusters()), plabels(new labels()){
	ifs.open(filename);
	std::string buf;
	if(!ifs){
		std::cout << "ERROR : flle [" << filename << "] does not exist." << std::endl;
		exit(1);
	}

	/* reading labels in cluster.out */
	int num_label;
	double x,y,z;
	while( ifs >> buf ){
		if( buf == "--" ) break;

		num_label = std::stoi(buf);
		ifs >> x >> y >> z;

		Eigen::Vector3d coord;
		coord << x, y, z;
		plabels->push_back( std::pair<int, Eigen::Vector3d>(num_label, coord) );
	}

	/* reading clusters in cluster.out */
	int index;
	int nbody;
	int num_in_cluster;
	int num_of_cluster_per_site;
	while( ifs >> index >> num_in_cluster >> nbody >> num_of_cluster_per_site ) {

		auto it = std::find(vec_index_read.begin(), vec_index_read.end(), index);
		if( vec_index_read.size()>0 and it==vec_index_read.end() ) {
			std::getline(ifs, buf);
			continue;
		}

		num_of_cluster_per_site /= nbody;
		// std::cout << index << ":" << num_in_cluster << ":" << num_of_cluster_per_site << std::endl;
		std::vector<std::vector<std::vector<int>>> vec_clusters_for_index(nbody);
		for(int i=0; i<nbody; ++i){
			std::vector<std::vector<int>> vec_clusters_for_a_site(num_of_cluster_per_site);
			for(int j=0; j<num_of_cluster_per_site; ++j){
				std::vector<int> vec_sites_of_clusters_for_a_site(num_in_cluster);
				for(int k=0; k<num_in_cluster; ++k){
					ifs >> vec_sites_of_clusters_for_a_site[k];
				}
				vec_sites_of_clusters_for_a_site.erase(vec_sites_of_clusters_for_a_site.begin());
				vec_clusters_for_a_site[j] = std::move(vec_sites_of_clusters_for_a_site);
			}
			vec_clusters_for_index[i] = std::move(vec_clusters_for_a_site);
		}
		pall_clusters->push_back(vec_clusters_for_index);
	}

}

ParseEcicar::ParseEcicar (const char* filename) {
	ifs.open(filename);
	std::string buf;
	if(!ifs){
		std::cout << "ERROR : flle [" << filename << "] does not exist." << std::endl;
		exit(1);
	}

	while( getline(ifs, buf) ){
		std::vector<double> vec_eci;

		std::stringstream ss(buf);
		int i;
		ss >> i;

		double eci;
		while( ss >> eci ) vec_eci.push_back(eci);
		ecicar[i] = std::move(vec_eci);

		index.push_back(i);
	}

}

ParsePoscar::ParsePoscar(const char *filename){
  std::ifstream ifs(filename);

	if ( !ifs ){
		std::cout << "ERROR : flle [" << filename << "] does not exist." << std::endl;
		exit(1);
	}

	std::string tmp_str;
	std::getline(ifs, this->comment);
	std::getline(ifs, tmp_str);
	this->unit = std::stod(tmp_str);

	for (int i=0; i<3; ++i){
		double a,b,c;
		ifs >> a >> b >> c;
		lattice_basis.col(i) << a, b, c;
		std::getline(ifs, tmp_str);
	}
	lattice_basis *= this->unit;

	std::string buf_num_atom_type;
	std::vector<std::string> vec_num_atom_type;
	std::getline( ifs, buf_num_atom_type);
	boost::algorithm::trim(buf_num_atom_type);
	boost::algorithm::split(vec_num_atom_type, buf_num_atom_type, boost::is_space());

	for (int i=0; i<vec_num_atom_type.size(); ++i){
		int num_atom_type;
		if( vec_num_atom_type[i].size() == 0 ) continue; /* through space */
		try {
			num_atom_type = std::stoi(vec_num_atom_type[i]);
		} catch (std::invalid_argument e) {
			break;
		}
		atom_types.push_back(num_atom_type);
	}

	std::getline(ifs, coordinate_type);

	int num_all_atoms = std::accumulate(atom_types.begin(), atom_types.end(), 0);

	for( int i=0; i<num_all_atoms; ++i ){
		double a,b,c;
		ifs >> a >> b >> c;
		std::pair<int, Eigen::Vector3d> atom(i, Eigen::Vector3d(a, b, c));
		atoms.push_back(atom);
		std::getline(ifs, tmp_str);
	}

}

Eigen::Matrix3d ParsePoscar::getLatticeBasis() const{
	return lattice_basis;
};

std::vector<int> ParsePoscar::getAtomTypes() const{
	return atom_types;
};

std::vector<std::pair<int, Eigen::Vector3d>> ParsePoscar::getAtoms() const{
	return atoms;
};

void ParsePoscar::expandPoscar(int expand_x, int expand_y, int expand_z){

	std::vector<double> spins;
	for(int i=0; i<this->atom_types.size(); ++i) {
		for(int j=0; j<this->atom_types[i]; ++j){
			spins.push_back(i);
		}
	}

	int N_unit = this->atoms.size();
	// std::vector<std::pair<int, Eigen::Vector3d>> new_positions;
	// std::vector<std::pair<int, Eigen::Vector3d>> new_positions_unit;

	using Position = std::tuple<int, int, Eigen::Vector3d>;

	/* site_num, atom_type, position */
	std::vector<Position> new_positions;
	std::vector<Position> new_positions_unit;

	Eigen::Matrix3d div_mat;
	div_mat << 1./double(expand_x),0,0,
						 0,1./double(expand_y),0,
						 0,0,1./double(expand_z);
	for( int i=0; i<atoms.size(); ++i ){
		//  std::pair<int, Eigen::Vector3d> -> spin_type and position
		new_positions_unit.push_back( std::make_tuple( i, spins[i], div_mat * atoms[i].second ));
	}

	int count = 0;
	for(int x=0; x<expand_x; ++x){
		for(int y=0; y<expand_y; ++y){
			for(int z=0; z<expand_z; ++z){
				Eigen::Vector3d block_vec;
				block_vec << double(x)/double(expand_x), double(y)/double(expand_y), double(z)/double(expand_z);
				for(int i=0; i<N_unit; ++i){
					auto position = block_vec + std::get<2>(new_positions_unit[i]);
					new_positions.push_back( std::make_tuple( i, spins[i], position ) );
					++count;
				}
			}
		}
	}

	std::sort( std::begin(new_positions), std::end(new_positions), [](const Position& x, const Position& y) -> int { return std::get<1>(x) > std::get<1>(y);});

	std::vector<int> new_atom_types;
	for(int i=0; i<this->atom_types.size(); ++i){
		int n = std::count_if(new_positions.begin(), new_positions.end(),  [i](const Position& position)->bool{return (std::get<1>(position) == i);} );
		new_atom_types.push_back(n);
	}

	std::vector<std::pair<int, Eigen::Vector3d>> new_atoms;
	for( const auto& position : new_positions ){
		new_atoms.push_back( std::pair<int, Eigen::Vector3d> ( std::get<0>(position), std::get<2>(position) ));
	}

	this->atoms = new_atoms;
	this->atom_types = new_atom_types;
	Eigen::Matrix3d max_ex_mat;
	max_ex_mat << expand_x,0,0, 0,expand_y,0, 0,0,expand_z;
	this->lattice_basis = max_ex_mat * this->lattice_basis;

}
