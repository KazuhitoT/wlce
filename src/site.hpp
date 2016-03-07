
#include <Eigen/Core>
#include <unordered_map>
#include <vector>

class Site {
private :
	const int site_num;
	const Eigen::Vector3d coordinate;
	std::unordered_map<double, std::vector<std::shared_ptr<Site>>> linked_site;

public :
	Site(int _n, const Eigen::Vector3d& _c):site_num(_n), coordinate(_c) {};

	bool operator <(const Site& rhs) {
		return this->site_num < rhs.site_num;
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
