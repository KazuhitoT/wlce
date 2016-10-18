#include "./site.hpp"
#include "./lattice.hpp"

const int Site::getLatticeNum(){return this->lattice->getLatticeNum();};

void Site::setLinkedSite(std::shared_ptr<Site> site, double distance, const Eigen::Vector3d& relative_coord){
		linked_site.insert( RelativeSite(site, distance, relative_coord) );
};

double myfloor( double dSrc, int iLen )
{
	double dRet;
	dRet = dSrc * pow(10.0, iLen);
	dRet = (double)(int)(dRet);
	return dRet * pow(10.0, -iLen);
}

Eigen::Vector3d validCoordinate(const Eigen::Vector3d& lhs, const Eigen::Vector3d& rhs){
	Eigen::Vector3d d_vec = lhs - rhs;
	for( int i=0; i<3; ++i ){ /* boundary condition */
		if( d_vec(i) > (0.5) )       {d_vec(i) = d_vec(i)-1;}
		else if( d_vec(i) < -(0.5) ) {d_vec(i) = 1+d_vec(i);}
	}
	return d_vec;
};

Eigen::Vector3d validCoordinate(Eigen::Vector3d lhs){
	for( int i=0; i<3; ++i ){ /* boundary condition */
		if( lhs(i) >  (0.5) )      {lhs(i) = lhs(i)-1;}
		else if( lhs(i) < -(0.5) ) {lhs(i) = 1+lhs(i);}
	}
	return lhs;
};

std::vector<std::vector<double>> getAllTriplets(const std::vector<Eigen::Vector3d>& points, const Eigen::Matrix3d& lattice, double maxd, const std::vector<double>& distances){
	std::vector<std::vector<double>> result;
	for(int i=0; i<points.size(); ++i){
		for(int j=(i+1); j<points.size(); ++j){
			for(int k=(j+1); k<points.size(); ++k){
				auto coord_ij = validCoordinate(points[i], points[j]);
				auto coord_jk = validCoordinate(points[j], points[k]);
				auto coord_ki = coord_ij + coord_jk;
				double distance1 = ( lattice * coord_ij ).norm();
				double distance2 = ( lattice * coord_jk ).norm();
				double distance3 = ( lattice * coord_ki ).norm();
				if( distance1 > maxd or distance2 > maxd or distance3 > maxd
						or distance1 < 0.00001 or distance2 < 0.00001  or distance3 < 0.00001) return std::vector<std::vector<double>>{};
				auto it1 = find_if(distances.begin(), distances.end(), [distance1](const double d){
					if( std::abs(distance1-d) < 0.00001 ) return true;
					else return false;
				});
				auto it2 = find_if(distances.begin(), distances.end(), [distance2](const double d){
					if( std::abs(distance2-d) < 0.00001 ) return true;
					else return false;
				});
				auto it3 = find_if(distances.begin(), distances.end(), [distance3](const double d){
					if( std::abs(distance3-d) < 0.00001 ) return true;
					else return false;
				});
				if( it1==distances.end() or it2==distances.end() or it3==distances.end()) return std::vector<std::vector<double>>{};
				std::vector<double> lines = {*it1, *it2, *it3};
				sort(lines.begin(), lines.end());
				result.push_back(lines);
			}
		}
	}

	for(int i=2; i>=0; --i){
		stable_sort(result.begin(), result.end(), [i](const std::vector<double>& lhs, const std::vector<double>& rhs){
			// std::cout <<  lhs[i] << " < " <<  rhs[i] << std::endl;
			if (std::abs(lhs[i] - rhs[i])<0.00001) return false;
			else return lhs[i] < rhs[i];
		});
	}
	//
	// for(const auto& i : result) {
	// 	std::cout  << "[ ";
	// 	for(const auto& j : i){
	// 		std::cout << j << " ";
	// 	}
	// 	std::cout  << " ]";
	// }
	// std::cout << std::endl;
	return result;
}


std::vector<TripletIndex> getAllTripletIndex(
	const std::vector<std::shared_ptr<Site>>& sites,
	const Eigen::Matrix3d& lattice,
	double maxd,
	const std::vector<PairIndex>& vec_pair_index
){

	std::vector<TripletIndex> result;
	for(int i=0; i<sites.size(); ++i){
		for(int j=(i+1); j<sites.size(); ++j){
			for(int k=(j+1); k<sites.size(); ++k){

				auto coord_ij = validCoordinate(sites[i]->getCoordinate(), sites[j]->getCoordinate());
				auto coord_jk = validCoordinate(sites[j]->getCoordinate(), sites[k]->getCoordinate());
				auto coord_ki = coord_ij + coord_jk;
				double distance1 = ( lattice * coord_ij ).norm();
				double distance2 = ( lattice * coord_jk ).norm();
				double distance3 = ( lattice * coord_ki ).norm();
				if( distance1 > maxd or distance2 > maxd or distance3 > maxd
						or distance1 < 0.00001 or distance2 < 0.00001  or distance3 < 0.00001) continue;

				auto pair_index1 = PairIndex(sites[i], sites[j], distance1);
				auto pair_index2 = PairIndex(sites[j], sites[k], distance2);
				auto pair_index3 = PairIndex(sites[k], sites[i], distance3);

				auto index_it1 = std::find(vec_pair_index.begin(), vec_pair_index.end(), pair_index1);
				auto index_it2 = std::find(vec_pair_index.begin(), vec_pair_index.end(), pair_index2);
				auto index_it3 = std::find(vec_pair_index.begin(), vec_pair_index.end(), pair_index3);

				if( index_it1==vec_pair_index.end() or index_it2==vec_pair_index.end() or index_it3==vec_pair_index.end() ){
					continue;
				};

				std::vector<PairIndex> vec_pair_index = {*index_it1, *index_it2, *index_it3};
				sort(vec_pair_index.begin(), vec_pair_index.end());
				result.push_back( TripletIndex(vec_pair_index) );
			}
		}
	}

	for(int i=2; i>=0; --i){
		stable_sort(result.begin(), result.end(), [i](const TripletIndex& lhs, const TripletIndex& rhs){
				return lhs.vec_pairindex[i] < rhs.vec_pairindex[i];
		});
	}

	return result;
}


const std::shared_ptr<Site> Site::not_found;
