#include "./parser.hpp"

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
	int j = 1;  /* for empty cluster */
	for(int i=0, imax=this->getContent().size(); i<imax; ++i){
		if( index[j] > imax ){
			std::cerr << "ERROR : indices does not correspond to those of [" << filename  << "]" << std::endl;
			std::cerr << "        index " << index[j] << " in ecicar does not exit." << std::endl;
			exit(1);
		}
		if( (index[j] -1 )!=i ) continue;
		multiplicity[index[j]] = std::pair<int, int>(this->getContent(i,0), this->getContent(i,1));
		j++;
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
