#include "./conf2corr.hpp"

Conf2corr::Conf2corr(char* filename,
	 std::vector<double> _spinposcar,
	 std::vector<double> _spince,
	 std::shared_ptr<labels> _plabels,
	 std::shared_ptr<allclusters> _ptr,
	 std::shared_ptr<basisfunc>   _pbasis_functions,
	 std::shared_ptr<indexorders> _pindex_orders
 ) : pbasis_functions(new basisfunc()), pindex_orders(new indexorders())
 {
	spinposcar = _spinposcar;
	spince = _spince;
	pall_clusters = _ptr;

	this->setSpins(_plabels);
	this->setCompositions();

	#ifndef SEED
	std::random_device rnd;
	mt.seed(rnd());
	#endif

	rnd_real  = std::uniform_real_distribution<double>(0.0, 1.0);
	rnd_int_N = std::uniform_int_distribution<int>(0,255);
	rnd_int_spince_index = std::uniform_int_distribution<int>(0,1);

	if( _pindex_orders == nullptr )	setIndexOrders();
	else this->pindex_orders = _pindex_orders;
	if( _pbasis_functions == nullptr ) setBasisCoefficient();
	else this->pbasis_functions = _pbasis_functions;
	this->setBasisFunction();

	setInitialCorrelationFunction();
}

Conf2corr::Conf2corr(char* filename,
	std::shared_ptr<Input> _in,
	std::shared_ptr<labels> _plabels,
	std::shared_ptr<allclusters> _pall_clusters,
	std::shared_ptr<basisfunc>   _pbasis_functions,
	std::shared_ptr<indexorders> _pindex_orders
) : pbasis_functions(new basisfunc()), pindex_orders(new indexorders()), pposcar_spin(new ParsePoscar(filename)){
	_in->setData("SPINCE", spince, true);
	_in->setData("SPINPOSCAR", spinposcar, true);
	_in->setData("CHECKINSIDESPIN", spin_for_check_inside);

	if( spin_for_check_inside != -10000 ){
		auto it = std::find(spince.begin(), spince.end(), spin_for_check_inside);
		if( it != spince.end() ){
			spin_for_check_inside = *it;
		} else {
			std::cerr << " ERROR : spin in CHECKINSIDESPIN does not exist.";
			exit(1);
		}
	}

	_in->setData("CHECKINSIDEINDEX", index_for_check_inside);
	pall_clusters = _pall_clusters;

	this->setSpins(_plabels);
	this->setCompositions();

	#ifndef SEED
	std::random_device rnd;
	mt.seed(rnd());
	#endif

	rnd_real  = std::uniform_real_distribution<double>(0.0, 1.0);
	rnd_int_N = std::uniform_int_distribution<int>(0,255);
	rnd_int_spince_index = std::uniform_int_distribution<int>(0,1);

	if( _pindex_orders == nullptr )	setIndexOrders();
	else this->pindex_orders = _pindex_orders;
	if( _pbasis_functions == nullptr ) setBasisCoefficient();
	else this->pbasis_functions = _pbasis_functions;
	this->setBasisFunction();

	this->setInitialCorrelationFunction();
};

void Conf2corr::setSpins(std::shared_ptr<labels> plabels){
	const std::vector<std::pair<int, Eigen::Vector3d>> atoms_spin = pposcar_spin->getAtoms();
	std::vector<int> atom_types = pposcar_spin->getAtomTypes();

	this->spins.resize(atoms_spin.size());
	for(int i=0, imax=atoms_spin.size(); i<imax; ++i){
		auto it = find_if(plabels->begin(), plabels->end(), [i, &atoms_spin](const std::pair<int, Eigen::Vector3d>& obj){
			if( (obj.second - atoms_spin[i].second).norm() < 0.00001 ) return true;
			else return false;
		});
		int count=-1;
		for(int index=i; index>=0; index-=atom_types[count]){
			count++;
		}
		this->spins[it->first] = count;
	}

	if( atoms_spin.size() != spins.size() ){
		std::cerr << "ERROR : the number of atoms in [pposcar_spin->in]  != total spins in [pposcar_spin->in]" << std::endl;
		std::cerr << atoms_spin.size() << " " << spins.size() << std::endl;
		exit(1);
	}

}

void Conf2corr::setCompositions(){
	for(int i=0; i<spince.size(); ++i ){
		 double composition = std::count(spins.begin(), spins.end(), spince[i]) / double(spins.size());
		 compositions.push_back(composition);
	}
	compositions_before = compositions;
}

void Conf2corr::setIndexOrders(){

	int num_composition = spince.size();

	unsigned int max_num_cluster=0;
	for(const auto& i: *pall_clusters ){
		if( i.size()>0 and i[0].size()>0 and max_num_cluster < (i[0][0].size()+1) )
			max_num_cluster = i[0][0].size()+1;
	}

	// std::cout << "max_num_cluster  " << max_num_cluster << std::endl;

	std::vector<std::vector<std::vector<int>>> tmp;
	for(int n=1; n<=max_num_cluster; ++n){
		std::vector<std::vector<int>> permutation;
		for(int i=0; i<std::pow(double(num_composition), double(n)); ++i){
			int ten_ans = i;
			int n_ans = 0;
			for (int j=0; ten_ans>0 ; j++)
			{
					n_ans = n_ans+(ten_ans%num_composition)*std::pow(10.0, double(j));
					ten_ans = ten_ans/num_composition;
			}
			std::ostringstream sout;
			sout << std::setfill('0') << std::setw(n) << n_ans;
			std::string s = sout.str();
			std::vector<int> tmp;
			for(int j=0; j<s.size(); ++j){
				tmp.push_back((int)(s[j]-'0'));
			}
			std::sort(tmp.begin(), tmp.end());
			if( tmp[0] != 0 ) permutation.push_back(tmp);
		}
		std::sort(permutation.begin(), permutation.end());
		permutation.erase(std::unique(permutation.begin(), permutation.end()), permutation.end());

		auto base_combination = permutation;
		for(auto& i : base_combination) std::sort(i.begin(), i.end());
		std::sort(base_combination.begin(), base_combination.end(), [num_composition](const std::vector<int>& lhs, const std::vector<int>& rhs){
			double lsum = 0;
			for(int i=0; i<lhs.size(); ++i) lsum += lhs[i] * std::pow(num_composition, i);
			double rsum = 0;
			for(int i=0; i<rhs.size(); ++i) rsum += rhs[i] * std::pow(num_composition, i);
			return lsum < rsum;
		});

		base_combination.erase(std::unique(base_combination.begin(), base_combination.end()), base_combination.end());

		std::vector<std::vector<std::vector<int>>> combination;
		for( auto data : base_combination){
			std::vector<std::vector<int>> tmp;
		  do{
		    // std::cout << "[ " << data[0];
		    // for(unsigned int i=1; i<data.size(); ++i){
		    //   std::cout << ", " << data[i];
		    // }
		    // std::cout << " ]" << std::endl;
				tmp.push_back(data);
		  }while(next_permutation(data.begin(), data.end()));
			combination.push_back(tmp);
		}
		pindex_orders->push_back(combination);
	}
}

void Conf2corr::setBasisFunction(){

	all_calculated_basis_functions.push_back(std::vector<double>());

	for(int order=1; order<this->compositions.size(); ++order){
		std::vector<double> vec_basis_functions; /* this index corresponds to the index of spince*/

		for(const auto& spin_spince : spince){
			double result = 0;
			/* vec[basis][degree] vec[basis]->polynomials */
			for(int i=0; i<(*pbasis_functions)[order].size(); ++i){
				result += (*pbasis_functions)[order][i] * pow(spin_spince, i);
			}
			vec_basis_functions.push_back(result);
		}

		this->all_calculated_basis_functions.push_back(vec_basis_functions);
	}
}


void Conf2corr::setInitialCorrelationFunction(){
	std::vector<std::vector<double>> corr;
	corr.push_back(std::vector<double>{1});	/* for empty cluster */
	for(int i=0, imax=(*pall_clusters).size(); i<imax; ++i){  // i == cluster index
		std::vector<double> tmp_corr;
		int num_of_cluster = 1;
		int num_in_cluster = 1;
		if( (*pall_clusters)[i].size()>0 ) {
			num_of_cluster = (*pall_clusters)[i][0].size();
			if( num_of_cluster>0 ){
				num_in_cluster = (*pall_clusters)[i][0][0].size()+1;
			}
		} else {
			num_of_cluster = 1;
			num_in_cluster = 1;
		}

		if(num_of_cluster==0 and num_in_cluster==1){
			num_of_cluster = 1;
			num_in_cluster = 1;
		}

		for(const auto& index_order : (*pindex_orders)[num_in_cluster-1] ){
			double corr_averaged_same_index = 0;
			// index_order == [[1112], [1121], [1211],...]
			//  orders == [1112]
			for(const auto& orders : index_order ){
				double average_spin_prod = 0;
				for(int site=0; site<(*pall_clusters)[i].size(); ++site) {
					if( (*pall_clusters)[i][site].size() == 0 ) { /*  point cluster */
						average_spin_prod += all_calculated_basis_functions[orders[0]][spins[site]];
					} else {
						for(const auto& site_clusters : (*pall_clusters)[i][site] ) {
							double spin_prod = all_calculated_basis_functions[orders[0]][spins[site]];
							for(int k=0; k<site_clusters.size(); ++k){
								spin_prod *= all_calculated_basis_functions[orders[k+1]][spins[site_clusters[k]]];
								// std::cout << orders[k+1] <<" " << spins[site_clusters[k]] << std::endl;
							}
							average_spin_prod += spin_prod;
						}
					}
				}
				corr_averaged_same_index += average_spin_prod;
			}
			corr_averaged_same_index /= (double)index_order.size() * (double)num_of_cluster * (double)spins.size() ;
			tmp_corr.push_back(corr_averaged_same_index);
		}
		corr.push_back(tmp_corr);
	}
	this->correlation_functions = corr;
}


//  !! note : point cluster is troublesome
void Conf2corr::setCorrelationFunction_flip(int lattice_point, int after_spin){

	int before_spin;

	while(1){
		if( lattice_point >= 0 ){
			before_spin = this->spins[lattice_point];
			this->spins[lattice_point] = after_spin;
		} else {
			lattice_point = this->rnd_int_N(mt);
			before_spin = this->spins[lattice_point];
			after_spin  = before_spin;
			while( before_spin == after_spin ){
				after_spin = this->rnd_int_spince_index(mt);
			}
			this->spins[lattice_point] = after_spin;
		}

		if( index_for_check_inside.size() > 0 and isInNthNearestNeighborPair(lattice_point) ) {
			this->spins[lattice_point] = before_spin;
			lattice_point = -1;
		} else {
			break;
		}
	}

	this->compositions[before_spin] -= 1/(double)this->spins.size();
	this->compositions[after_spin]  += 1/(double)this->spins.size();

	this->vec_changed_spins.emplace_back(std::array<int,3>{lattice_point, before_spin, after_spin});

	for(int i=0, imax=(*pall_clusters).size(); i<imax; ++i){  // i == cluster index
		int num_of_cluster = 1;
		int num_in_cluster = 1;
		if( (*pall_clusters)[i].size()>0 ) {
			num_of_cluster = (*pall_clusters)[i][0].size();
			if( num_of_cluster>0 ){
				num_in_cluster = (*pall_clusters)[i][0][0].size()+1;
			}
		} else {
			num_of_cluster = 1;
			num_in_cluster = 1;
		}

		if(num_of_cluster==0 and num_in_cluster==1){
			num_of_cluster = 1;
			num_in_cluster = 1;
		}

		int count_index_order = 0;
		for(const auto& index_order : (*pindex_orders)[num_in_cluster-1] ){
			double delta_corr_averaged_same_index = 0;
			// index_order == [[1112], [1121], [1211],...]
			//  orders == [1112]
			for(const auto& orders : index_order ){
				double delta_average_spin_prod = 0;
				if( (*pall_clusters)[i][lattice_point].size() == 0 ) { /*  point cluster */
					delta_average_spin_prod += this->all_calculated_basis_functions[orders[0]][after_spin] - this->all_calculated_basis_functions[orders[0]][ before_spin];
				} else {
					for(const auto& site_clusters : (*pall_clusters)[i][lattice_point] ) {
						// std::cout << (*pall_clusters)[i][lattice_point].size() << " " << num_of_cluster << std::endl;
						// assert( (*pall_clusters)[i][lattice_point].size() == num_of_cluster );
						double delta_spin_after  = this->all_calculated_basis_functions[orders[0]][after_spin];
						double delta_spin_before = this->all_calculated_basis_functions[orders[0]][before_spin];
						double num_equiv_site = 1;
						for(int k=0; k<site_clusters.size(); ++k){
							// std::cout << site_clusters[k] << " " ;
							double delta_spin_prod = this->all_calculated_basis_functions[orders[k+1]][spins[site_clusters[k]]];
							delta_spin_after *= delta_spin_prod;
							// std::cout << lattice_point << " : " <<  site_clusters[k] << std::endl;
							if( lattice_point == site_clusters[k] ){
								delta_spin_before *= this->all_calculated_basis_functions[orders[k+1]][before_spin];
								++num_equiv_site;
							}
							else delta_spin_before *= delta_spin_prod;
						}
						// std::cout << std::endl;
						delta_average_spin_prod += (delta_spin_after - delta_spin_before)/num_equiv_site;
						assert( (delta_spin_after-delta_spin_before)<=2 );
					}
				}
				// std::cout << delta_average_spin_prod << std::endl;
				delta_corr_averaged_same_index += delta_average_spin_prod;
			}
			delta_corr_averaged_same_index /= (double)index_order.size() * (double)spins.size() * (double)num_of_cluster / (double)(num_in_cluster);
			// std::cout << " -- " << this->correlation_functions[i+1][count_index_order] << " " << num_in_cluster << " " << num_of_cluster << std::endl;
			this->correlation_functions[i+1][count_index_order] += delta_corr_averaged_same_index;
			// std::cout << "--" << delta_corr_averaged_same_index << " -- " << this->correlation_functions[i+1][count_index_order] << " -- "<< i << std::endl;
			// assert( this->correlation_functions[i+1][count_index_order] <= 1 and this->correlation_functions[i+1][count_index_order] >= -1 );
			++count_index_order;
		}
	}
}

void Conf2corr::setCorrelationFunction_exchange(){
	int exchanged_spins[2];

	exchanged_spins[0] = rnd_int_N(mt);
	exchanged_spins[1] = rnd_int_N(mt);
	while( spins[exchanged_spins[0]] == spins[exchanged_spins[1]] )
		exchanged_spins[1] = rnd_int_N(mt);

	this->vec_changed_spins.emplace_back(std::array<int,3>{exchanged_spins[0], spins[exchanged_spins[0]], spins[exchanged_spins[1]]});
	this->vec_changed_spins.emplace_back(std::array<int,3>{exchanged_spins[1], spins[exchanged_spins[1]], spins[exchanged_spins[0]]});

	double tmp_spin =  spins[exchanged_spins[0]];
	this->setCorrelationFunction_flip(exchanged_spins[0], spins[exchanged_spins[1]]);
	this->setCorrelationFunction_flip(exchanged_spins[1], tmp_spin);
}

double Conf2corr::calcCorrelationFunctionNorm(double p){
	double result = 0.0;
	for(const auto& corr_index : this->correlation_functions ){
		for(const auto& corr : corr_index){
			result += std::pow(std::fabs(corr), p);
		}
	}
	return std::pow(result, 1.0/p);
}

/**
 * setBasisCoefficient through gramm-schmitt
 *
 * @author
 * @version
 * @param unknown $factory
 * @return
 *

	if three component with -1 0 1,
 * 1 0 0
 * 0 1.22474 0
 * -1.41421 0 2.12132
 */
void Conf2corr::setBasisCoefficient(){
	std::function<double(const std::vector<double>&, const std::vector<double>&)> trace;
	auto tmp_spince = this->spince;
	trace = [&tmp_spince](const std::vector<double>& lhs, const std::vector<double>& rhs){
		const unsigned int r = tmp_spince.size();
		double result = 0;
		for(int i=0; i<r; ++i){
			for(int j=0; j<r; ++j){
				for(int k=0; k<r; ++k){
					result += lhs[j] * rhs[k] * pow(tmp_spince[i], (double)(j+k));
				}
			}
		}
		return result/(double)r;
	};

	(*pbasis_functions).resize(this->spince.size());
	const unsigned int r = tmp_spince.size();
	for(int i=0; i<r; ++i){
		std::vector<double> tmp(r, 0);
		if( i==0 ){
			tmp[0] = 1;
			(*pbasis_functions)[0] = tmp;
			continue;
		}
		std::vector<double> pow_s_m(r, 0);
		pow_s_m[i] = 1;
		for(unsigned int m=0; m<i; ++m){
			double tmp_coef =trace((*pbasis_functions)[m],	pow_s_m);
			for(unsigned int j=0; j<tmp.size(); ++j){
				tmp[j] -= tmp_coef * (*pbasis_functions)[m][j];
			}
		}
		tmp[i] = 1;

		double normalize = sqrt(trace(tmp, tmp));
		std::for_each(tmp.begin(), tmp.end(), [normalize](double &x){ x /= normalize; });
		(*pbasis_functions)[i] = tmp;
	}
};


bool Conf2corr::isInNthNearestNeighborPair(int lattice_point){
	if( this->spins[lattice_point] != spin_for_check_inside ) return false;

	for( auto index : index_for_check_inside ){
		--index; // for empty cluster
		for(const auto& site_clusters : (*(this->pall_clusters))[index][lattice_point] ) {
			for(int k=0; k<site_clusters.size(); ++k){
				if( this->spins[site_clusters[k]] == spin_for_check_inside ) return true;
			}
		}
	}
	return false;
}


void Conf2corr::setMemento(){
	this->vec_changed_spins.clear();
}

void Conf2corr::Memento(){
	for(const auto& changed_spin : this->vec_changed_spins){
		this->setCorrelationFunction_flip(changed_spin[0], changed_spin[1]);
	}
};

void Conf2corr::dispCorr(){
	for(int i=0, imax=correlation_functions.size(); i<imax; ++i){
		std::cout << i << "  ";
		for( const auto& j : correlation_functions[i] )
			std::cout << j << " ";
		std::cout << std::endl;
	}
}

void Conf2corr::outputPoscar(std::string prefix){
	std::ofstream ofs("poscar."+prefix);
	ofs << "POSCAR" << std::endl;
	ofs << "1.0" << std::endl;
	ofs.setf(std::ios::right, std::ios::adjustfield);
	ofs.setf(std::ios_base::fixed, std::ios_base::floatfield);
	ofs << std::setprecision(15);
	ofs << std::setw(20);

	const auto lattice = pposcar_spin->getLatticeBasis();
	for(int i=0; i<3; ++i){
		for(int j=0; j<3; ++j){
			ofs << lattice(i,j) << " ";
		}
		ofs << std::endl;
	}

	for( int i=0; i<this->spince.size(); ++i ){
		int count = 0;
		for( int j=0; j<this->spins.size(); ++j ){
			if( spince[i] == spins[j] ) ++count;
		}
		ofs << count << " ";
	}
	ofs << std::endl;

	ofs << "Direct" << std::endl;
	std::vector<std::vector<Eigen::Vector3d>> atom_types_coords(pposcar_spin->getAtomTypes().size());

	const auto& atoms = pposcar_spin->getAtoms();
	for( int i=0; i<atoms.size(); ++i ){
		for( int j=0; j<this->spince.size(); ++j ){
			if( this->spins[i] == this->spince[j] ){
				atom_types_coords[j].push_back(atoms[i].second);
				break;
			}
		}
	}

	for(const auto& atoms_coord: atom_types_coords){
		for(const auto& atom_coord : atoms_coord){
			ofs << atom_coord[0] << " " << atom_coord[1] << " " << atom_coord[2] << std::endl;
		}
	}

	return;
}
