#include <algorithm>
#include <tuple>
#include "./parser.hpp"

ParseClusterIn_::ParseClusterIn_ (const char* f):filename(f) {
	ifs.open(filename);
	std::string buf;
	if(!ifs){
		std::cout << "ERROR : flle [" << filename << "] does not exist." << std::endl;
		exit(1);
	}

	std::string str;
	ifs >> str;
	while( !ifs.eof() ) {
		ifs >> str;

	}

}


Parser::Parser (const char* f):filename(f) {
	ifs.open(filename);
	std::string buf;
	if(!ifs){
		std::cout << "ERROR : flle [" << filename << "] does not exist." << std::endl;
		exit(1);
	}
	while(ifs && getline(ifs, buf)){
		if(!buf.size()) continue;
		std::vector<std::string> line;
		boost::algorithm::trim(buf);
		boost::algorithm::split(line, buf, boost::is_space());
		std::vector<double> tmp(line.size());
		for(int i=0, imax=tmp.size(); i<imax; ++i){
			tmp[i] = atof(line[i].c_str());
		}
		content.push_back(tmp);
	}
}

Parser::Parser (const char* f, const std::vector<int>& vindex) : filename(f) {
	ifs.open(filename);
	if(!ifs){
		std::cerr << "ERROR : flle [" << filename << "] does not exist." << std::endl;
		exit(1);
	}
	std::string buf;
	int index_count = 0;
	while(ifs && getline(ifs, buf)){
		++index_count;
		auto it = find(vindex.begin(), vindex.end(), index_count);
		if( it == vindex.end() ) continue;

		std::vector<std::string> line;
		boost::algorithm::trim(buf);
		boost::algorithm::split(line, buf, boost::is_space());
		std::vector<double> tmp(line.size());
		for(int i=0, imax=tmp.size(); i<imax; ++i){
			tmp[i] = atof(line[i].c_str());
		}
		content.push_back(tmp);
	}
}


ParseEcicar::ParseEcicar (const char* filename):Parser(filename) {
	for(int i=0, imax=this->getContent().size(); i<imax; ++i){
		std::vector<double> tmp;
		for(const auto& j: this->getContent(i)){
			tmp.push_back(j);
		}
		tmp.erase(tmp.begin()); /* delete index */
		ecicar[this->getContent(i,0)] = tmp;
		index.push_back(getContent(i,0));
	}
}

ParseMultiplicityIn::ParseMultiplicityIn (const char* filename):Parser(filename) {
	multiplicity[0] = std::pair<int, int>(1, this->getContent(0,1));
	for(int i=0, imax=this->getContent().size(); i<imax; ++i){
		multiplicity[(i+1)] = std::pair<int, int>(this->getContent(i,0), this->getContent(i,1));
	}
}

ParseMultiplicityIn::ParseMultiplicityIn (const char* filename, const std::vector<int>& index):Parser(filename) {
	multiplicity[0] = std::pair<int, int>(1, this->getContent(0,1));
	for(int i=1; i<index.size(); ++i){
		if( (index[i]-1) > this->getContent().size() ){
			std::cerr << "ERROR : indices does not correspond to those of [" << filename  << "]" << std::endl;
			std::cerr << "        index " << index[i]-1 << " in ecicar does not exit." << std::endl;
			exit(1);
		}
		multiplicity[index[i]] = std::pair<int, int>(this->getContent(index[i]-1,0), this->getContent(index[i]-1,1));
	}
}

ParseClusterIn::ParseClusterIn (const char* filename, const std::map<int , std::pair<int, int> >& multiplicity )
	: Parser(filename), pclusters(new std::vector<std::vector<std::vector<std::vector<int>>>>())
 {
	for(const auto i : multiplicity) index.push_back(i.first);
	int j = 1;  /* for empty cluster */
	for(int i=0, imax=this->getContent().size(); i<imax; ++i){
		std::vector<std::vector<std::vector<int> > > c;
		std::vector<std::vector<int> > b;
		std::vector<int> a;
		int num_in_cluster = multiplicity.at(index[j]).first;
		/* ここのnum_of_clusterだけ配位数を意味している */
		/* NOTE :: multiplicity.at(0).secondを原子数として代用している */
		int num_of_cluster = multiplicity.at(index[j]).second / multiplicity.at(0).second * num_in_cluster;

		/* NOTE: kmaxも含めないと最後まで含まれない*/
		for(int k=1, kmax=this->getContent(i).size(); k<=kmax; ++k){
			if((k % (num_of_cluster * num_in_cluster)) == 0){
				b.push_back(a);
				c.push_back(b);
				std::vector<int> ().swap(a);
				std::vector<std::vector<int> > ().swap(b);
				continue;
			}else if((k % num_in_cluster) == 0){
				b.push_back(a);
				std::vector<int> ().swap(a);
				continue;
			}
			a.push_back(this->getContent(i, k));
		}
		pclusters->push_back(c);
		std::vector<std::vector<std::vector<int> > > ().swap(c);
		j++;
		if((index.size()) == j)
			break;
	}
	std::vector<std::vector<std::vector<std::vector<int>>>> (*pclusters).swap(*pclusters);

	this->clearContent();
}

ParseClusterIn::ParseClusterIn (
	const char* filename,
	const std::vector<int>& _index,
	const std::map<int , std::pair<int, int> >& multiplicity)
	: Parser(filename, _index),
		index(_index),
	  pclusters(new std::vector<std::vector<std::vector<std::vector<int>>>>())
 {

	int j = 1;  /* for empty cluster */
	for(int i=0, imax=this->getContent().size(); i<imax; ++i){
		std::vector<std::vector<std::vector<int> > > c;
		std::vector<std::vector<int> > b;
		std::vector<int> a;
		int num_in_cluster = multiplicity.at(index[j]).first;
		/* ここのnum_of_clusterだけ配位数を意味している */
		int num_of_cluster = multiplicity.at(index[j]).second / multiplicity.at(0).second * num_in_cluster;

		/* NOTE: kmaxも含めないと最後まで含まれない*/
		for(int k=1, kmax=this->getContent(i).size(); k<=kmax; ++k){
			if((k % (num_of_cluster * num_in_cluster)) == 0){
				b.push_back(a);
				c.push_back(b);
				std::vector<int> ().swap(a);
				std::vector<std::vector<int> > ().swap(b);
				continue;
			}else if((k % num_in_cluster) == 0){
				b.push_back(a);
				std::vector<int> ().swap(a);
				continue;
			}
			a.push_back(this->getContent(i, k));
		}
		pclusters->push_back(c);
		std::vector<std::vector<std::vector<int> > > ().swap(c);
		j++;
		if((index.size()) == j)
			break;
	}
	std::vector<std::vector<std::vector<std::vector<int>>>> (*pclusters).swap(*pclusters);

	this->clearContent();
}

void ParseClusterIn::checkClusterIn(){
	for(int i=0; i<pclusters->size(); ++i){
		if((*pclusters)[i].size() == 0) continue;
		for(int j=0; j<(*pclusters)[i].size(); ++j){
			std::cout << j << " " << (*pclusters)[i][j].size() << std::endl;
		}
	}
}

void outputPoscar(double lattice[3][3], double position[][3], int N, std::string prefix){
	std::ofstream ofs("poscar."+prefix);
	ofs << "POSCAR" << std::endl;
	ofs << "1.0" << std::endl;
	ofs.setf(std::ios::right, std::ios::adjustfield);
	ofs.setf(std::ios_base::fixed, std::ios_base::floatfield);
	ofs << std::setprecision(15);
	ofs << std::setw(20);

	for(int i=0; i<3; ++i){
		for(int j=0; j<3; ++j){
			ofs << lattice[i][j] << " ";
		}
		ofs << std::endl;
	}
	ofs << N << std::endl;
	ofs << "Direct" << std::endl;
	for(int i=0; i<N; ++i){
		for(int j=0; j<3; ++j){
			ofs << position[i][j] << " ";
		}
		ofs << std::endl;
	}
};


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

ParseLabels::ParseLabels(const char* filename) : plabels(new labels())
{
	std::ifstream ifs(filename);
	std::string buf;

	if ( !ifs ){
		std::cout << "ERROR : flle [" << filename << "] does not exist." << std::endl;
		exit(1);
	}

	do{
		int label;
		ifs >> label;
		double a,b,c;
		ifs >> a >> b >> c;
		// std::cout << label << " " << a << " " << " " << b << " " << c << std::endl;
		std::pair<int, Eigen::Vector3d> tmp(label, Eigen::Vector3d(a, b, c));
		this->plabels->push_back(tmp);
	}while(ifs && getline(ifs, buf));
};
