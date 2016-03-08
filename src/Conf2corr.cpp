#include "./Conf2corr.hpp"

Conf2corr::Conf2corr(char* filename,
	 vector<double> _spinposcar,
	 vector<double> _spince,
	 Ensemble _ensemble,
	 vector<double> _chemical_potential,
	 int _setrandom,
	 std::string _input_spin_filename
 ){
	spinposcar = _spinposcar;
	spince = _spince;
	ensemble = _ensemble;
	chemical_potential = _chemical_potential;

	assert( spinposcar.size()>0 );
	assert( spince.size()>0 );

	setSpins(filename);

	if( _setrandom>0 and _input_spin_filename.size() ){
		cout << "ERROR : Both SETRANDOM and SPININPUT are available." << endl;
		cout << "Plese make available only one configuration." << endl;
		exit(1);
	} else if(_input_spin_filename.size()){
		this->setSpinsFromDat(_input_spin_filename, -DBL_MAX, DBL_MAX);
		this->setCorrelationFunctionFromClucar();
	} else if(_setrandom){
		this->setSpinsRandom();
	}

	setInitialCondition();
}

void Conf2corr::setInitialCondition(){
	setEci();
	Conf2corr::setIndexOrders();
	Conf2corr::setBasisCoefficient();
	setCorrelationFunctionFromClucar();
	setTotalEnergy();
	Conf2corr::setInitializeCov(this->getCorrelationFunctions().size());
	mt = std::mt19937(rd());
	rnd_int_N = std::bind( std::uniform_int_distribution<int>(0, spins.size()-1), mt);
	rnd_real = std::bind( std::uniform_real_distribution<double>(0.0, 1.0), mt);
}

Conf2corr::Conf2corr(const Conf2corr &obj){
	spins = obj.getSpins();
	spins_before = obj.getBeforeSpins();
	eci = obj.getEci();
	correlation_functions = obj.getCorrelationFunctions();
	correlation_functions_before = obj.getBeforeCorrelationFunctions();
	totalEnergy = obj.getTotalEnergy();
	mahalanobis_distance = obj.getMahalanobisDistance();
}

void Conf2corr::setIndexOrders(){
	if( index_orders.size()>0 ){
		// cerr << "error : Conf2corr::setIndexOrders" << endl;
		// exit(1);
		// return ;
		index_orders = vector<vector<vector<vector<int>>>> {};
	}

	vector<int> list;
	for(int i=1; i<spince.size(); ++i){
		list.push_back(i);
	}

	unsigned int max_num_cluster=0;
	for(const auto& i: dimcar.getDimcar()){
		if( max_num_cluster<i.second.first )
			max_num_cluster = i.second.first;
	}

	for(int i=1; i<spince.size() ;++i){
		list.push_back(i);
	}

	// index_orders.push_back(vector<vector<double>>{0});
	for(int n=1; n<=max_num_cluster; ++n){
		vector<vector<vector<int>>> tmp;
		vector<vector<int>> combination;
		RepeatedCombinationGenerator<int> g(list.begin(), list.end(), n);
		do{
			combination.push_back(g.data());
		}while(g.next());

		for(int i=0;i<combination.size();++i){
			vector<vector<int>> permutation;
			do{
				permutation.push_back(combination[i]);
			}while( next_permutation( combination[i].begin(),combination[i].end() ) );
			tmp.push_back(permutation);
		}
		// cout << combination.size()<< endl;
		index_orders.push_back(tmp);
	}
}


/* setNewConfとは独立させる 	*/
void Conf2corr::setSpinsRandom(){
	random_shuffle(spins.begin(), spins.end());
	// this->setCorrelationFunctionFromClucar();
	// this->setTotalEnergy();
}

bool Conf2corr::setSpinsFromDat(const string filename, double emin, double emax, bool isSkipNoSpins){
	ifstream ifs(filename.c_str());
	// cout << filename << endl;exit(1);
	if(!ifs){
		ifs.close();
		cout << "ERROR : flle [" << filename << "] does not exist." << endl;
		exit(1);
	}

	bool isSpinSearched = false;
	string buf;
	vector<string> v;
	while(ifs && getline(ifs, buf)){
		if(buf.size() == 0) continue;
		boost::algorithm::trim(buf);
		boost::algorithm::split(v, buf, boost::is_space());
		double e = stod(v.at(1));
		v.erase(v.begin());
		v.erase(v.begin());

		if(v.size() != this->getSpins().size()){
			cout << "ERROR : input spin size in [" << filename << "] differs from [POSCAR.spin]" << endl;
			exit(1);
		}
		if(emin <= e and e <= emax){
			#ifdef DEBUG_SPIN_DISP
			cout << "INPUT :" << endl;
			cout << buf << endl;
			#endif
			for(int i=0; i<v.size(); ++i){
				this->spins.at(i) = stod(v.at(i));
			}
			isSpinSearched = true;
			break;
		}
	}
	ifs.close();

	if(isSkipNoSpins){

	} else if(!isSpinSearched){
		cout << "ERROR : no spin configuration satisfies the condition in [wang-landau.in]" << endl;
		exit(1);
	}
	spins_before = spins;

	if(isSpinSearched)
		return true;
	else
		return false;
}

void Conf2corr::setEci(){
	#ifdef DEBUG
	cout << "setEci() start" << endl;
	#endif
	if(eci.size() > 0)
		return ;
	vector< vector<double> > eci_tmp = ecicar.getContent();
	for(int i=0, imax=eci_tmp.size(); i<imax; ++i){
		const int index = eci_tmp[i][0];
		eci_tmp[i].erase(eci_tmp[i].begin());
		eci.push_back(pair<int, vector<double>>(index, eci_tmp[i]));
	}
}

void Conf2corr::setSpins(char* filename){
	#ifdef DEBUG
	cout << "setSpins(char* filename) start" << endl;
	#endif
	ParsePoscar poscarSpin(filename);
	vector<std::pair<int, Eigen::Vector3d>> pos_cluster = poscarCluster.getAtoms();
	vector<std::pair<int, Eigen::Vector3d>> pos_spin    = poscarSpin.getAtoms();

	vector<int> num_atoms = poscarSpin.getNumAtoms();

	// NOTE : 計算時間が最大NC2
	for(int i=0, imax=pos_cluster.size(); i<imax; ++i){
		vector<pair<int, Eigen::Vector3d>>::iterator it;
		for(it=pos_spin.begin(); it!=pos_spin.end(); it++){
			if( (pos_cluster[i].second - it->second).norm() < 0.00000001 ) {
				int count=-1;
				for(int index=(*it).first; index>=0; index-=num_atoms[count]){
					count++;
				}
				this->spins.push_back(spinposcar[count]);
				break;
			}
		}
	}
	if(this->spins.size() != pos_spin.size()){
		cerr << "ERROR atom numbers in POSCAR.spin is incorrect." << endl;
		exit(1);
	}
}

void Conf2corr::setGrowthFactor(){
	double integer, decimal;
	decimal = modf( sqrt(correlation_functions.size()-1) , &integer);
	if( fabs(decimal) >= 0.00001 ){
		cout << "error : Conf2corr::setGrowthFactor " << endl;
		exit(-1);
	}

	// auto corr_func = this->correlation_functions;
	// for(int i=0; i<100; ++i){
	// 	std::random_shuffle ( corr_func.begin(), corr_func.end() );
	// 	MatrixXd A=Map<MatrixXd>(&corr_func[1], integer, integer);
	// 	PartialPivLU<MatrixXd> lu(A);
	// 	MatrixXd U=lu.matrixLU().triangularView<Upper>();
	// 	cout << this->correlation_functions[1] << " " << U.cwiseAbs().maxCoeff() / A.cwiseAbs().maxCoeff() << " "
	// 				<< U.cwiseAbs().maxCoeff() <<" " << A.cwiseAbs().maxCoeff() <<  endl;
	// }

	// this->growth_factor = result / corr_func.size();
	// 離散的な信号 -> 規則化傾向なら
	// それぞれのcorr_funcをexpするというのは？
	// ノイズ
}

void Conf2corr::setRadialDistributions(){
	vector<vector<double>> comp_spins(spinposcar.size());
	for(int i=0; i<this->spins.size(); ++i){
		double spin = spins[i];  //  for lamda expression
		auto fi = std::find_if(this->spinposcar.begin(), this->spinposcar.end(),
				[spin](double &sp){return((spin == sp));});
		// if(fi != this->spinposcar.end()){
		// 	std::cout << "std::find_if" << std::endl;
		// }
		size_t index = std::distance(this->spinposcar.begin(), fi);
		comp_spins[index].push_back(i);
	}

	auto all_atoms = Conf2corr::poscarCluster.getAtoms();

	vector<vector<Eigen::Vector3d>> comp_atoms(comp_spins.size());
	for( int i=0; i<comp_spins.size(); ++i ){
		for( auto c : comp_spins[i] ){
			auto fi = std::find_if(all_atoms.begin(), all_atoms.end(),
					[c](pair<int, Vector3d> &atom){return((c == atom.first));});
			comp_atoms[i].push_back(fi->second);
		}
	}

	auto axis = Conf2corr::poscarCluster.getAxis();
	Eigen::Vector3d center(0.5, 0.5, 0.5);
	double max_distance = sqrt( (axis * center).transpose() * (axis * center) );

	// the width of bins = max_distance / N
	const int N = this->spins.size();
	const double bin_width = max_distance / (double)N;
	this->radial_distributions.resize(comp_atoms.size());
	for(int i=0; i<this->radial_distributions.size(); ++i){
		this->radial_distributions[i].resize(N+1);
		for(auto atom : comp_atoms[i]){
			for(auto target : comp_atoms[i]){
				Vector3d delta = target - atom;
				for(int tmp=0; tmp<delta.size(); ++tmp){
					if( abs(delta[tmp]) < 0.5 ){
						delta[tmp] = abs(delta[tmp]);
					} else{
						delta[tmp] = 1 - abs(delta[tmp]);
					}
					// cout << delta[tmp] << endl;
				}
				double distance = sqrt( (axis * delta).transpose() * (axis * delta) );
				int index = (max_distance - distance)/bin_width;
				// cout << index << endl;
				++this->radial_distributions[i].at(index);
			}
		}
	}

	for(auto ll : this->radial_distributions[0]){
		cout << ll << endl;
	}
		// int num_spins = 0;
	// for(int i=0; i<this->radial_distributions.size(); ++i){
	// 	for(int j=0; j<this->spins.size(); ++j){
	// 		if( spins[num_spins] spinposcar[i] ){
	//
	// 		}
	// 	}
	// }
};

/* vec[basis][degree] vec[basis]->polynomials */
double Conf2corr::getBasisFunction(/* order */  int order, /* spince_num */ int spin){
	double result = 0;
	for(int i=0; i<basis_functions[order].size(); ++i){
		result += basis_functions[order][i] * pow(spin, i);
	}
	return result;
};


//  clucarからcorrelaton_functionを計算（重い）
void Conf2corr::setCorrelationFunctionFromClucar(){
	#ifdef DEBUG
	cout << "setCorrelationFunctionFromClucar() start" << endl;
	#endif

	vector<int> index = ecicar.getIndex();
	vector<vector<double>> corr;
	corr.push_back(vector<double>{1});	/* for empty cluster */
	for(int i=1, imax=index.size(); i<imax; ++i){  // iはdimcarの行数
		vector<double> tmp_corr;
		vector<int> clucar_tmp = clucar.getClucar(index[i]);
		int num_in_cluster = dimcar.getDimcar(index[i]).first;
		int num_of_cluster = dimcar.getDimcar(index[i]).second;
		for(const auto& index_order : index_orders[num_in_cluster-1] ){
			double corr_averaged_same_index = 0;
			// index_order == [[1112], [1121], [1211],...]
			//  orders == [1112]
			for(const auto& orders : index_order ){
				double spin_prod = 1;
				double average_spin_prod = 0;
				for(int j=0, jmax=clucar_tmp.size(); j<jmax; ++j) {
					if( j!=0 and (j%num_in_cluster) == 0){		/* j=0の時に+1されるから，あとで引く*/
						average_spin_prod += spin_prod;
						spin_prod = 1;
					}
					spin_prod *= getBasisFunction(orders[j%num_in_cluster], spins[clucar_tmp[j]]);
				}
				average_spin_prod += spin_prod;
				corr_averaged_same_index += average_spin_prod / (double)num_of_cluster;
			}
			corr_averaged_same_index /= (double)index_order.size();
			tmp_corr.push_back(corr_averaged_same_index);
		}
		corr.push_back(tmp_corr);
	}
	this->correlation_functions = corr;
	this->correlation_functions_before = corr;
}

void Conf2corr::setCorrelationFunction(){
	this->setMemento();	/* save spins and corr  */
	if( this->ensemble == Ensemble::ALLOY )
		this->setCorrelationFunction_exchange();
	else if( this->ensemble == Ensemble::SPIN )
		this->setCorrelationFunction_flip(rnd_int_N());
}

void Conf2corr::setCorrelationFunction_flip(const int lattice_point){
	auto before_spin = this->spins[lattice_point];
	auto after_spin  = before_spin;
	auto tmp_spins = this->spince;
	while( before_spin == after_spin ){
		after_spin = tmp_spins[0];
		tmp_spins.erase(tmp_spins.begin());
	}
	this->spins[lattice_point] = after_spin;
	for(int i=1; i<eci.size(); ++i){		/* eciのfirstはindex */
		int num_in_cluster = dimcar.getDimcar(eci[i].first).first;
		int num_of_cluster = dimcar.getDimcar(eci[i].first).second;
		vector<vector<int> > tmp = mccar.getMccar(eci[i].first, lattice_point);
		int count_index_order = 0;
		for(const auto& index_order : index_orders[num_in_cluster-1] ){
			double corr_averaged_same_index_delta = 0;
			/*
			 * index_order == [[1112], [1121], [1211],...]
			 *  orders == [1112]
			 */
			 int tmp_count = 0;
			for(const auto& orders : index_order ){
				double sum_delta = 0;
				for(int k=0, kmax=tmp.size(); k<kmax; ++k){
					double spin_prod_delta = getBasisFunction(orders[0], spins[lattice_point]) - getBasisFunction(orders[0], spins_before[lattice_point]);
					/*
					 * if lattice is binary fcc and now considering 1nn pair,
					 * k is index, i.e., first row in ecicar
					 * tmp[k] == [64, 73, 100, 109 , ... ] and its size is 12
					 */
					assert( tmp[k].size() == (orders.size()-1) );
					for(int l=0, lmax=tmp[k].size(); l<lmax; ++l){
						assert( tmp[k][l] != lattice_point );
						spin_prod_delta *= getBasisFunction(orders[l+1], spins[tmp[k][l]]);
					}
					sum_delta += spin_prod_delta;
				}
				corr_averaged_same_index_delta += sum_delta;
			}
			correlation_functions[i][count_index_order] += corr_averaged_same_index_delta / double(num_of_cluster) / double(index_order.size());
			++count_index_order;
		}
	}
}

void Conf2corr::setCorrelationFunction_exchange(){
	const vector<int> swaped_spins = getSwapSpinsIndex();
	#ifdef DEBUG
	cout << "swaped_spins" << endl;
	for(int i=0; i<swaped_spins.size(); ++i){
		cout << swaped_spins[i] << " ";
	}
	cout << endl;
	#endif

	#ifdef DEBUG
	cout << "---------- before corr ----------" << endl;
	dispCorr();
  	#endif

	for(int i=1; i<eci.size(); ++i){		/* eciのfirstはindex */
		int num_in_cluster = dimcar.getDimcar(eci[i].first).first;
		int num_of_cluster = dimcar.getDimcar(eci[i].first).second;
		for(int j=0; j<swaped_spins.size(); ++j){
			vector<vector<int> > tmp = mccar.getMccar(eci[i].first, swaped_spins[j]);
			vector<double> corr_tmp;
			int count_index_order = 0;
			for(const auto& index_order : index_orders[num_in_cluster-1] ){
				double corr_averaged_same_index_before = 0;
				double corr_averaged_same_index_after  = 0;
				// index_order == [[1112], [1121], [1211],...]
				// orders == [1112]
				for(const auto& orders : index_order ){
					double sum_before = 0;
					double sum_after  = 0;
					for(int k=0, kmax=tmp.size(); k<kmax; ++k){
						double spin_prod_before = getBasisFunction(orders[0], spins_before[swaped_spins[j]]);
						double spin_prod_after  = getBasisFunction(orders[0], spins[swaped_spins[j]]);
						if( j==1 ){ // if another spin is included in mccar, continue
							vector<int>::iterator it = find(tmp[k].begin(), tmp[k].end(), swaped_spins[0]);
							if( it != tmp[k].end() ){
								continue;
							}
						}
						for(int l=0, lmax=tmp[k].size(); l<lmax; ++l){
							spin_prod_before *= getBasisFunction(orders[l+1], spins_before[tmp[k][l]]);
							spin_prod_after  *= getBasisFunction(orders[l+1], spins[tmp[k][l]]);
						}
						sum_before += spin_prod_before;
						sum_after  += spin_prod_after;
					}
					corr_averaged_same_index_before += sum_before;
					corr_averaged_same_index_after  += sum_after;
				}
				correlation_functions[i][count_index_order] += (corr_averaged_same_index_after - corr_averaged_same_index_before) / double(num_of_cluster) / double(index_order.size());
				++count_index_order;
			}
		}
	}

	#ifdef DEBUG
	cout << "---------- after  corr ----------" << endl;
	dispCorr();
	cout << "---------- confirm corr ----------" << endl;
	setCorrelationFunctionFromClucar();
	dispCorr();
  #endif
}

void Conf2corr::setCorrelationFunction_withMultiAtomSwap(){
    // std::random_device rnd;
    // std::mt19937 mt(rnd());
    // std::uniform_int_distribution<> rand(1, in->max_swap);
    // int imax = rand(mt);
		int imax = 1;
	for(int i=0; i<imax; ++i)
		this->setCorrelationFunction();
}

void Conf2corr::setTotalEnergy(){
	totalEnergy = 0;
	for(int i=0, imax=eci.size(); i<imax; ++i){
		for(int j=0, jmax=eci[i].second.size(); j<jmax; ++j){
			totalEnergy += correlation_functions[i][j] * eci[i].second[j];
		}
	}
	/*  setCorrに組み込んだほうが早い  */
	if(chemical_potential.size()){
		vector<double> compositions;
		for(const auto& spin : this->spince){
			 double composition = std::count(spins.begin(), spins.end(), spin) / double(spins.size());
			 compositions.push_back(composition);
		}
		assert( compositions.size() == chemical_potential.size() );
		for(int i=0; i<compositions.size(); ++i){
			totalEnergy -= chemical_potential[i] * compositions[i];
		}
	}
}

void Conf2corr::updateCovarianceMatrix(){
	if( Conf2corr::num_allconf<=2 ){
		cout << "ERROR updateCovarianceMatrix needs at least 2 configuration for initial covariance matrix" << endl;
		exit(1);
	}
	// Conf2corr::cov_matrix_before = Conf2corr::cov_matrix;
	// Conf2corr::ave_vector_before = Conf2corr::ave_vector;
	// for(int i=0; i<Conf2corr::cov_matrix.cols(); ++i){
	// 	for(int j=i; j<Conf2corr::cov_matrix.rows(); ++j){
	// 		double delta_ave_i = Conf2corr::ave_vector(i) - Conf2corr::ave_vector_before(i);
	// 		double delta_ave_j = Conf2corr::ave_vector(j) - Conf2corr::ave_vector_before(j);
	// 		double tmp1 = Conf2corr::cov_matrix_before(i,j) * (double)(num_allconf-2) / (double)(num_allconf-1) ;
	// 		double tmp2 = (this->getCorrelationFunctions().at(i) - Conf2corr::ave_vector_before(i)) * (this->getCorrelationFunctions().at(j) - Conf2corr::ave_vector_before(j)) / (double)(num_allconf-1);
	// 		double tmp3 = delta_ave_j * (this->getCorrelationFunctions().at(i) - Conf2corr::ave_vector_before(i)) / (double)(num_allconf-1);
	// 		double tmp4 = delta_ave_i * (this->getCorrelationFunctions().at(j) - Conf2corr::ave_vector_before(j)) / (double)(num_allconf-1);
	// 		double tmp5 = delta_ave_i * delta_ave_j * num_allconf / (double)(num_allconf-1);
	// 		Conf2corr::cov_matrix(i,j) = tmp1 + tmp2 - tmp3 - tmp4 + tmp5;
	// 		Conf2corr::cov_matrix(j,i) = Conf2corr::cov_matrix(i,j);
	// 	}
	// }
	// ++Conf2corr::num_allconf;
};

void Conf2corr::setAverageVector(){
	cout << "setAverageVector" << endl;
	exit(1);
	// Conf2corr::ave_vector_before = Conf2corr::ave_vector;
	// for(int i=1; i<Conf2corr::ave_vector.size(); ++i){
	// 	Conf2corr::ave_vector[i] = (Conf2corr::ave_vector_before[i] * (num_allconf-1) + this->getCorrelationFunctions().at(i) ) / (double)(num_allconf);
	// }
};

/* NOTE :: cov_matrixにempty_clusterを入れると要素が全部ゼロになる行，列が存在して逆行列が発散する．  */
/*         以上の理由により，最初の行と列を抜いてから逆行列を計算する                                 */
void Conf2corr::setMahalanobisDistance(){
	MatrixXd block  = Conf2corr::cov_matrix.block(1, 1, Conf2corr::cov_matrix.cols()-1, Conf2corr::cov_matrix.rows()-1);
	vector<double> tmp;
	for(const auto& i: this->correlation_functions )
		tmp.insert(tmp.end(), i.begin(), i.end());
	VectorXd v = Map<VectorXd>(&tmp[1], tmp.size()-1);
	for(int i=1; i<Conf2corr::ave_vector.size(); ++i){
		v(i-1) -= Conf2corr::ave_vector(i);
	}
	this->mahalanobis_distance = sqrt(v.transpose() * block.inverse() * v);
};

vector<int> Conf2corr::getSwapSpinsIndex(){
	vector<int> swaped_spins(2,0);

	swaped_spins[0] = rnd_int_N();
	swaped_spins[1] = rnd_int_N();

	while( spins[swaped_spins[0]] == spins[swaped_spins[1]] )
		swaped_spins[1] = rnd_int_N();
	swap(spins[swaped_spins[0]], spins[swaped_spins[1]]);
	return swaped_spins;
}

Conf2corr &Conf2corr::operator=(const Conf2corr &ob)
{
    if (this == &ob)
    	return *this;
  	this->spins = ob.getSpins();
		this->spins_before = ob.getBeforeSpins();
    this->correlation_functions = ob.getCorrelationFunctions();
		this->correlation_functions_before = ob.getBeforeCorrelationFunctions();
		this->mahalanobis_distance = ob.getMahalanobisDistance();
		this->totalEnergy = ob.getTotalEnergy();
    return *this;
};

void Conf2corr::Memento()
{
	this->spins = this->spins_before;
	this->correlation_functions = this->correlation_functions_before;
	this->totalEnergy = this->totalEnergy_before;
};

void Conf2corr::setMemento()
{
	this->spins_before = this->spins;
	this->correlation_functions_before = this->correlation_functions;
	this->totalEnergy_before = this->totalEnergy;
};

//
// double Conf2corr::getProperty(int corrDosIndex = -1){
// 	if(corrDosIndex > 0){
// 		return this->getCorrelationFunctions().at(corrDosIndex);
// 	} else {
// 		return this->getTotalEnergy();
// 	}
// };

void Conf2corr::outputCorr(string filename){
	ofstream ofs(filename, std::ios::app);
	for( auto i : correlation_functions )
		for( auto j : i )
			ofs << j << " ";
		ofs << endl;
	ofs.close();
}

void Conf2corr::outputCorrSpin(int num, int mode, string filename){
	ofstream ofs(filename, std::ios::app);
	ofs << num << " ";
	ofs.precision(10);
	ofs << this->getTotalEnergy() << " ";
	if( mode==OUTPUT_ENERGY_AND_SPIN ){
		for( auto j : this->getSpins() ){
			ofs << j << " ";
		}
	}
	ofs << endl;
	ofs.close();
}

void Conf2corr::outputPoscar(string filename){
	vector<vector<int>> sorted_spins(spinposcar.size());
	for(int i=0; i<spins.size(); ++i){
		for(int j=0; j<spinposcar.size(); ++j){
			if( spins[i] == spinposcar[j] ){
				sorted_spins[j].push_back(i);
				break;
			}
		}
	}

	auto atoms = poscarCluster.getAtoms();
	ofstream ofs(filename);
	ofs.precision(10);
	ofs << poscarCluster.getComment() << endl;
	ofs << "1.0" << endl;
	ofs << poscarCluster.getAxis() << endl;
	for(auto same_spins : sorted_spins )
		ofs << same_spins.size() << " ";
	ofs << endl;
	ofs << poscarCluster.getCoordinateType() << endl;
	for(auto same_spins : sorted_spins ){
		for(auto spin_num : same_spins){
			for(int i=0; i<atoms[spin_num].second.size(); ++i)
				ofs << atoms[spin_num].second[i] << " ";
			ofs << endl;
		}
	}
	ofs.close();
}

void Conf2corr::setCovarianceMatrix(string filename){
	const int num_loop = this->getSpins().size()*this->getSpins().size();
	Conf2corr::setInitializeCov(this->getCorrelationFunctions().size());
	vector<double> all_corrs;
	Conf2corr tmp(*this);
	for(int i=0; i<num_loop; ++i){
		tmp.setCorrelationFunction();
		vector<vector<double>> corrs_tmp = tmp.getCorrelationFunctions();
		for(const auto& tmp : corrs_tmp)
			all_corrs.insert(all_corrs.end(), tmp.begin(), tmp.end());
	}

	// 行ごとにcorrをいれてる
	if(all_corrs.size()==0){
		return;
	}

	MatrixXd A=Map<MatrixXd>(&all_corrs[0], all_corrs.size()/num_loop, num_loop);
	VectorXd V = VectorXd::Zero(this->getCorrelationFunctions().size());
	for(int i=0; i<this->getCorrelationFunctions().size(); ++i){
		V(i) = A.row(i).mean();
		for(int j=0; j<A.cols(); ++j){
			A(i,j) -= V(i);
		}
	}
	auto cov = A * A.transpose();
	Conf2corr::cov_matrix = cov/(double)(num_loop);
	Conf2corr::ave_vector = V;
	Conf2corr::num_allconf = num_loop;

	#ifdef DEBUG_COV_DISP
	cout << " -- average correlation functions (DEBUG) " << endl;
	cout << ave_vector << endl;
	cout << " -- covariance matrix (DEBUG) " << endl;
	cout << cov_matrix << endl;
	#endif
}

/*  Metropolis algo. */
void Conf2corr::setNewConf(){
	this->setCorrelationFunction();
	this->setTotalEnergy();
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
	const unsigned int r = Conf2corr::spince.size();
	Conf2corr::basis_functions.resize(r);
	for(int i=0; i<r; ++i){
		vector<double> tmp(r, 0);
		if( i==0 ){
			tmp[0] = 1;
			Conf2corr::basis_functions[0] = tmp;
			continue;
		}
		vector<double> pow_s_m(r, 0);
		pow_s_m[i] = 1;
		for(unsigned int m=0; m<i; ++m){
			double tmp_coef = Conf2corr::trace(Conf2corr::basis_functions[m],	pow_s_m);
			for(unsigned int j=0; j<tmp.size(); ++j){
				tmp[j] -= tmp_coef * Conf2corr::basis_functions[m][j];
			}
		}
		tmp[i] = 1;

		double normalize = sqrt(Conf2corr::trace(tmp, tmp));
		std::for_each(tmp.begin(), tmp.end(), [normalize](double &x){ x /= normalize; });
		basis_functions[i] = tmp;
	}
};

/**
 * setBasisCoefficient through gramm-schmitt
 *
 * @author
 * @version
 * @param vector<double> basis_functions_coefficient, spin^m
 * @return
 */
double Conf2corr::trace(const vector<double>& lhs, const vector<double>& rhs){
	const unsigned int r = Conf2corr::spince.size();
	double result = 0;
	for(int i=0; i<r; ++i){
		for(int j=0; j<r; ++j){
			for(int k=0; k<r; ++k){
				result += lhs[j] * rhs[k] * pow(Conf2corr::spince[i], (double)(j+k));
			}
		}
	}
	return result/(double)r;
}


void Conf2corr::dispCorr(){
	cout << "index    correlation_functions   (DEBUG)" << endl;
	for(int i=0, imax=correlation_functions.size(); i<imax; ++i){
		cout << eci[i].first << "  ";
		for( const auto& j : correlation_functions[i] )
			cout << j << " ";
		cout << endl;
	}
}

void Conf2corr::dispSpin(){
	for(auto i : this->spins) cout << " "  << spins[i];
	cout << endl;
}

void Conf2corr::dispEci(){
	cout << "    ECI   (DEBUG)" << endl;
	for(const auto& i : eci){
		cout << i.first << " ";
		for(const auto& j : i.second )
			cout << j << " ";
		cout << endl;
	}
}

void Conf2corr::dispCovarianceCoefficient(){
	int size = Conf2corr::cov_matrix.cols();
	MatrixXd coef_matrix = MatrixXd::Zero(size, size);
	coef_matrix(0,0) = 1;

	vector<double> std;
	for(int i=1; i<size; ++i){
		std.push_back(sqrt(Conf2corr::cov_matrix(i,i)));
	}
	for(int i=1; i<size; ++i){
		for(int j=1; j<size; ++j){
			coef_matrix(i,j) = Conf2corr::cov_matrix(i,j) / std.at(i-1) / std.at(j-1);
		}
	}
	cout << coef_matrix << endl;
}

void Conf2corr::dispBasisCoefficient(){
	for(const auto& i : basis_functions){
		for(const auto& j : i){
			cout << j << " ";
		}
		cout << endl;
	}
}

const ParseEcicar Conf2corr::ecicar("./ecicar");
const ParseDimcar Conf2corr::dimcar("./dimcar", ecicar.getIndex());
const ParseMccar  Conf2corr::mccar("./mccar.conf", ecicar.getIndex(), dimcar.getDimcar());
const ParseClucar Conf2corr::clucar("./clucar", ecicar.getIndex());

ParsePoscar Conf2corr::poscarCluster("./POSCAR.cluster");
vector<double> Conf2corr::spinposcar;
vector<double> Conf2corr::spince;
vector<double> Conf2corr::chemical_potential;
vector<vector<double>> Conf2corr::basis_functions;
vector<vector<vector<vector<int>>>> Conf2corr::index_orders;
MatrixXd Conf2corr::cov_matrix;
MatrixXd Conf2corr::cov_matrix_before;
VectorXd Conf2corr::ave_vector;
VectorXd Conf2corr::ave_vector_before;
int Conf2corr::num_allconf;
Ensemble Conf2corr::ensemble;

std::random_device Conf2corr::rd;
std::mt19937 Conf2corr::mt;
std::function<int ()> Conf2corr::rnd_int_N;
std::function<double ()> Conf2corr::rnd_real;
