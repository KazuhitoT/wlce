#include "./Parser.hpp"

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
	std::string buf;
	while(ifs && getline(ifs, buf)){
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
		ecicar[this->getContent(i,0)] = this->getContent(i,1);
		index.push_back(getContent(i,0));
	}
}

ParseMultiplicityIn::ParseMultiplicityIn (const char* filename, const std::vector<int>& index):Parser(filename) {
	multiplicity[0] = std::pair<int, int>(1, this->getContent(0,1));
	int j = 1;  /* for empty cluster */
	for(int i=0, imax=this->getContent().size(); i<imax; ++i){
		if( (index[j] -1 )!=i ) continue;
		multiplicity[index[j]] = std::pair<int, int>(this->getContent(i,0), this->getContent(i,1));
		#ifdef DEBUG
		std::cout << index[j] << " " << multiplicity[index[j]].first << " " << multiplicity[index[j]].second << std::endl;
		#endif
		j++;
	}
}

ParseClusterIn::ParseClusterIn (const char* filename, const std::vector<int>& index, const std::map<int , std::pair<int, int> >& multiplicity):Parser(filename) {
	/* indexを合わせるために，empty clusterのベクターを挿入しておく */
	std::vector<std::vector<std::vector<int> > > dummy;
	clusters.push_back(dummy);
	/* jはindex  multiplicityにあるやつ読み込む */
	int j = 1;  /* for empty cluster */
	for(int i=0, imax=this->getContent().size(); i<imax; ++i){
		if( (index[j] -1 )!=i ) {
			/* indexを合わせるために，空のベクターを挿入しておく */
			clusters.push_back(dummy);
			continue;
		}
		int lpoint = 0;
		std::vector<double> tmp = this->getContent(i);
		std::vector<std::vector<std::vector<int> > > c;
		std::vector<std::vector<int> > b;
		std::vector<int> a;
		int num_in_cluster = multiplicity.at(index[j]).first;
		/* ここのnum_of_clusterだけ配位数を意味している */
		int num_of_cluster = multiplicity.at(index[j]).second / multiplicity.at(0).second * num_in_cluster;
		/* NOTE: kmaxも含めないと最後まで含まれない*/
		for(int k=1, kmax=tmp.size(); k<=kmax; ++k){
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
			a.push_back(tmp[k]);
		}
		clusters.push_back(c);
		std::vector<std::vector<std::vector<int> > > ().swap(c);
		j++;
		if((index.size()) == j)
			break;
	}
	std::vector<std::vector<std::vector<std::vector<int> > > > (clusters).swap(clusters);
	this->clearContent();
}

void ParseClusterIn::checkClusterIn(){
	for(int i=0; i<clusters.size(); ++i){
		if(clusters[i].size() == 0) continue;
		for(int j=0; j<clusters[i].size(); ++j){
			std::cout << j << " " << clusters[i][j].size() << std::endl;
		}
	}
}
