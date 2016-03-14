#include "./Conf2corr.hpp"

/* for corrdump */
Conf2corr::Conf2corr(std::vector<double> _spins,
	 std::vector<double> _spinposcar,
	 std::vector<double> _spince,
	 std::shared_ptr<allclusters> _ptr,
	 std::shared_ptr<basisfunc>   _pbasis_functions,
	 std::shared_ptr<indexorders> _pindex_orders
 ) : pbasis_functions(new basisfunc()), pindex_orders(new indexorders())
 {
	spins = _spins;
	spins_before = spins;
	spinposcar = _spinposcar;
	spince = _spince;
	pall_clusters = _ptr;

	rnd_int_N = std::bind( std::uniform_int_distribution<int>(0, spins.size()-1), mt);
	rnd_real  = std::bind( std::uniform_real_distribution<double>(0.0, 1.0), mt);

	assert( spinposcar.size()>0 );
	assert( spince.size()>0 );

	if( !_pindex_orders )	setIndexOrders();
	else this->pindex_orders = _pindex_orders;
	if( !_pbasis_functions ) setBasisCoefficient();
	else this->pbasis_functions = _pbasis_functions;

	setInitialCorrelationFunction();
}

void Conf2corr::setIndexOrders(){

	int num_composition = spince.size();

	unsigned int max_num_cluster=0;
	for(const auto& i: *pall_clusters ){
		if( i.size()>0 and i[0].size()>0 and max_num_cluster < i[0][0].size() )
			max_num_cluster = i[0][0].size()+1;
	}

	std::cout << "max_num_cluster  " << max_num_cluster << std::endl;

	std::vector<std::vector<std::vector<int>>> tmp;
	for(int n=1; n<=max_num_cluster; ++n){
		std::vector<std::vector<int>> permutation;
		std::vector<int> tmp(n, 1);
		std::function<void(int, int, std::vector<int>&)> get_permutation;
		get_permutation = [&get_permutation, tmp, num_composition, &permutation](int max, int i, std::vector<int>& tmp) -> void
		{
			if( i==tmp.size() ){
				permutation.push_back(tmp);
			} else {
				for(int j = 1; j < num_composition; ++j){
					tmp[i] = j;
					get_permutation(max, i+1, tmp);
				}
			}
 		};
		get_permutation(n, 1, tmp);
		permutation.push_back( std::vector<int>(n, num_composition-1) );

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
		    std::cout << "[ " << data[0];
		    for(unsigned int i=1; i<data.size(); ++i){
		      std::cout << ", " << data[i];
		    }
		    std::cout << " ]" << std::endl;
				tmp.push_back(data);
		  }while(next_permutation(data.begin(), data.end()));
			combination.push_back(tmp);
		}
		pindex_orders->push_back(combination);
	}
}

//  !! note : point cluster is troublesome
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
		}
		if( i==0 ) num_of_cluster = 1;

		for(const auto& index_order : (*pindex_orders)[num_in_cluster-1] ){
			double corr_averaged_same_index = 0;
			// index_order == [[1112], [1121], [1211],...]
			//  orders == [1112]
			for(const auto& orders : index_order ){
				double average_spin_prod = 0;
				for(int site=0; site<(*pall_clusters)[i].size(); ++site) {
					if( (*pall_clusters)[i][site].size() == 0 ) { /*  point cluster */
						average_spin_prod += getBasisFunction(orders[0], spins[site]);
					} else {
						for(const auto& site_clusters : (*pall_clusters)[i][site] ) {
							double spin_prod = getBasisFunction(orders[0], spins[site]);
							for(int k=0; k<site_clusters.size(); ++k){
								spin_prod *= getBasisFunction(orders[k+1], spins[site_clusters[k]]);
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
	this->correlation_functions_before = corr;
}

//  !! note : point cluster is troublesome
void Conf2corr::setCorrelationFunction_flip(){
	const int lattice_point = this->rnd_int_N();
	auto before_spin = this->spins[lattice_point];
	auto after_spin  = before_spin;
	auto tmp_spins = this->spince;
	while( before_spin == after_spin ){
		after_spin = tmp_spins[0];
		tmp_spins.erase(tmp_spins.begin());
	}
	this->spins[lattice_point] = after_spin;

	for(int i=0, imax=(*pall_clusters).size(); i<imax; ++i){  // i == cluster index
		std::vector<double> tmp_corr;
		int num_of_cluster = 1;
		int num_in_cluster = 1;
		if( (*pall_clusters)[i].size()>0 ) {
			num_of_cluster = (*pall_clusters)[i][0].size();
			if( num_of_cluster>0 ){
				num_in_cluster = (*pall_clusters)[i][0][0].size()+1;
			}
		}
		if( i==0 ) num_of_cluster = 1;

		int count_index_order = 0;
		for(const auto& index_order : (*pindex_orders)[num_in_cluster-1] ){
			double delta_corr_averaged_same_index = 0;
			// index_order == [[1112], [1121], [1211],...]
			//  orders == [1112]
			for(const auto& orders : index_order ){
				double delta_average_spin_prod = 0;
				// for(const auto& site : (*pall_clusters)[i] ) {
				if( (*pall_clusters)[i][lattice_point].size() == 0 ) { /*  point cluster */
					delta_average_spin_prod += getBasisFunction(orders[0], after_spin) - getBasisFunction(orders[0], before_spin);
				} else {
					for(const auto& site_clusters : (*pall_clusters)[i][lattice_point] ) {
						double delta_spin_prod = getBasisFunction(orders[0], after_spin) - getBasisFunction(orders[0], before_spin);
						for(int k=0; k<site_clusters.size(); ++k){
							delta_spin_prod *= getBasisFunction(orders[k+1], spins[site_clusters[k]]);
						}
						delta_average_spin_prod += delta_spin_prod;
					}
				}
				delta_corr_averaged_same_index += delta_average_spin_prod;
			}
			delta_corr_averaged_same_index /= (double)index_order.size() * (double)spins.size() * (double)num_of_cluster / (double)num_in_cluster;
			this->correlation_functions[i+1][count_index_order] += delta_corr_averaged_same_index;
			++count_index_order;
		}
	}

}

/* vec[basis][degree] vec[basis]->polynomials */
double Conf2corr::getBasisFunction(/* order */  int order, /* spince_num */ int spin){
	double result = 0;
	for(int i=0; i<(*pbasis_functions)[order].size(); ++i){
		result += (*pbasis_functions)[order][i] * pow(spin, i);
	}
	return result;
};

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

void Conf2corr::dispCorr(){
	for(int i=0, imax=correlation_functions.size(); i<imax; ++i){
		std::cout << i << "  ";
		for( const auto& j : correlation_functions[i] )
			std::cout << j << " ";
		std::cout << std::endl;
	}
}
