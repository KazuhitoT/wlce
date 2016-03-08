
#include <Eigen/Core>
#include <unordered_map>
#include <vector>

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

std::vector<std::vector<double>> getAllTriplets(const std::vector<Eigen::Vector3d>& points, const Eigen::Matrix3d& lattice, double maxd){
	std::vector<std::vector<double>> result;
	for(int i=0; i<points.size(); ++i){
		for(int j=(i+1); j<points.size(); ++j){
			for(int k=(j+1); k<points.size(); ++k){
				double distance1 = ( lattice * validCoordinate(points[i], points[j]) ).norm();
				double distance2 = ( lattice * validCoordinate(points[j], points[k]) ).norm();
				double distance3 = ( lattice * validCoordinate(points[k], points[i]) ).norm();
				if( distance1 > maxd or distance2 > maxd or distance3 > maxd
						or distance1 < 0.00001 or distance2 < 0.00001  or distance3 < 0.00001) continue;
				std::vector<double> lines = {distance1, distance2, distance3};
				sort(lines.begin(), lines.end());
				result.push_back(lines);
			}
		}
	}

	for(int i=(points.size()-1); i>=0; --i){
		stable_sort(result.begin(), result.end(), [i](const std::vector<double>& lhs, const std::vector<double>& rhs){
			return lhs[i] < rhs[i];
		});
	}

	return result;
}

class Site {
private :

	class hash_vec3d {
	public:
	  double operator()(const Eigen::Vector3d &x) const {
	    const int C = 997;      // 素数
	    double t = 0;
	    for (int i = 0; i != x.size(); ++i) {
	        t = t * C + x[i];
	    }
	    return t;
	  }
	};

	const int site_num;
	const Eigen::Vector3d coordinate;
	std::unordered_map<double, std::vector<std::shared_ptr<Site>>> linked_site;

	/*  site-> pair< corrdination, site >  */
	// std::vector<std::pair<Eigen::Vector3d, std::shared_ptr<Site>>> site_relative;
	std::unordered_map<Eigen::Vector3d, std::shared_ptr<Site>, hash_vec3d> site_relative;

public :
	Site(int _n, const Eigen::Vector3d& _c):site_num(_n), coordinate(_c) {};

	bool operator <(const Site& rhs) {
		return this->site_num < rhs.site_num;
	}

	void setRelativeSite(const std::vector<std::shared_ptr<Site>>& all_sites){
		for( const auto& site : all_sites){
			// std::cout << site_num << " " << site->getSiteNum() << " - " <<  validCoordinate(site->getCoordinate(), this->coordinate).transpose() << std::endl;
			// site_relative.push_back( std::pair<Eigen::Vector3d, std::shared_ptr<Site>>(validCoordinate(site->getCoordinate(), this->coordinate), site) );
			site_relative[validCoordinate(site->getCoordinate(), this->coordinate)] = site;
		}
		assert( this->site_relative.size() == all_sites.size() );
	}

	std::shared_ptr<Site> getRelativeSite(const Eigen::Vector3d& coord){
		return this->site_relative[coord];
	}

	void setLinkedSite(double distance, const std::shared_ptr<Site> another_site){
	this->linked_site[distance].push_back(another_site);
	};

	const std::unordered_map<double, std::vector<std::shared_ptr<Site>>>& getLinkedSite() const {
		return linked_site;
	}

	const std::shared_ptr<Site>& getLinkedSite(double distance, const int num) const {
		if( this->linked_site.find(distance) == this->linked_site.end() ){
			bool is_found = false;
			for( const auto& i : this->linked_site ){
				if( std::abs( distance - i.first ) < 0.00001 ) {
					distance = i.first;
					is_found = true;
					break;
				}
			}
			/*  distance is not found  */
			if(!is_found)	return nullptr;
		}

		auto it = find_if( linked_site.at(distance).begin(), linked_site.at(distance).end(),
				[num](const std::shared_ptr<Site> s){ return ( s->getSiteNum() == num );} );
		if( it != linked_site.at(distance).end() )
			return *it;
		else
			return nullptr;
	};

	const std::vector<std::shared_ptr<Site>>& getLinkedSiteviaDisncace(const double distance) const {
		return this->linked_site.at(distance);
	};

	int getSiteNum() const { return this->site_num; };
	const Eigen::Vector3d& getCoordinate() const { return this->coordinate; };
	const double getDistance(double distance) const {
		if( this->linked_site.find(distance) == this->linked_site.end() ){
			for( const auto& i : this->linked_site ){
				if( std::abs( distance - i.first ) < 0.00001 ) {
					return i.first;
				}
			}
			/*  distance is not found  */
		} else {
		return distance;
		}
		return -1;
	}

	bool isCoordinate( const Eigen::Vector3d& obj ){
		for( int i=0; i<3; ++i ){
			if( std::abs( this->coordinate[i]-obj[i] ) > 0.00001 ) {
				return false;
			}
		}
		return true;
	}
};
