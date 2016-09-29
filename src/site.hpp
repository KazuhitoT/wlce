#ifndef _SITE_HPP
#define _SITE_HPP

#include <Eigen/Core>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <boost/mem_fn.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/composite_key.hpp>
#include "myindex.hpp"

class Lattice;
class Site;
class PairIndex;
class TripletIndex;

Eigen::Vector3d validCoordinate(const Eigen::Vector3d& lhs, const Eigen::Vector3d& rhs);

Eigen::Vector3d validCoordinate(Eigen::Vector3d lhs);

std::vector<std::vector<double>> getAllTriplets(const std::vector<Eigen::Vector3d>& points, const Eigen::Matrix3d& lattice, double maxd, const std::vector<double>& distances);

std::vector<TripletIndex> getAllTripletIndex(
	const std::vector<std::shared_ptr<Site>>& sites,
	const Eigen::Matrix3d& lattice,
	double maxd,
	const std::vector<PairIndex>& vec_pair_index
);


struct RelativeSite{
	double distance;
	std::shared_ptr<Site> site;
	double x,y,z;

	RelativeSite(std::shared_ptr<Site> _s, double _d, const Eigen::Vector3d& _r)
		: site(_s), distance(_d) {
			x = _r(0);
			y = _r(1);
			z = _r(2);
		}

	bool operator==(const RelativeSite &r) const {
		if( std::fabs(this->x-r.x) < 0.00001
				&& std::fabs(this->y-r.y) < 0.00001
				&& std::fabs(this->z-r.z) < 0.00001
			){
			return true;
		} else {
			return false;
		}
	}

	friend size_t hash_value(const RelativeSite &r){
		const int C = 997;	// prime_number
		double result = 0;
		result = result*C + r.x;
		result = result*C + r.y;
		result = result*C + r.z;
		return result;
	}
};


using namespace std;
using namespace boost::multi_index;

struct tag_distance {};
struct tag_relative_coord {};

typedef composite_key <
    RelativeSite,
		BOOST_MULTI_INDEX_MEMBER(RelativeSite, double, x),
		BOOST_MULTI_INDEX_MEMBER(RelativeSite, double, y),
		BOOST_MULTI_INDEX_MEMBER(RelativeSite, double, z)
> composite_key_site;

typedef composite_key_hash <
	std::hash<double>,
	std::hash<double>,
	std::hash<double>
> composite_key_hash_coord;

typedef composite_key_equal_to <
	std::equal_to<double>,
	std::equal_to<double>,
	std::equal_to<double>
> composite_key_equal_coord;

typedef multi_index_container<
    RelativeSite,
    indexed_by<
        hashed_non_unique<tag<tag_distance>, member<RelativeSite, double, &RelativeSite::distance> >,
				hashed_unique<composite_key_site, composite_key_hash_coord, composite_key_equal_coord>
    >
> LinkedSite;


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
	std::shared_ptr<Lattice> lattice;

	LinkedSite linked_site;

	const static std::shared_ptr<Site> not_found;

public :
	Site(int _n, const Eigen::Vector3d& _c):site_num(_n), coordinate(_c) {};
	Site(int _n, const Eigen::Vector3d& _c, std::shared_ptr<Lattice> _l):site_num(_n), coordinate(_c), lattice(_l){};

	bool operator <(const Site& rhs) {
		return this->site_num < rhs.site_num;
	}

	void setLinkedSite(const std::vector<std::shared_ptr<Site>>& all_sites, double pair_truncation);
	void setLinkedSite(const std::shared_ptr<Site> all_sites, double distance, const Eigen::Vector3d& relative_coord);

	std::pair<LinkedSite::iterator, LinkedSite::iterator> getLinkedSiteIterator(double distance){
		std::pair<LinkedSite::iterator, LinkedSite::iterator> pair_itr  = linked_site.equal_range(distance);
		return pair_itr;
	}

	std::shared_ptr<Site> getLinkedSite(const Eigen::Vector3d& coord){
		LinkedSite::nth_index<1>::type& a = linked_site.get<1>();
		LinkedSite::nth_index<1>::type::iterator itr = a.find(boost::make_tuple(coord(0), coord(1), coord(2)));
		if( itr != a.end() ) return itr->site;
		else {
			for( auto it=a.begin(); it!=a.end(); ++it){
				if( (validCoordinate(it->site->getCoordinate(), this->getCoordinate())-coord).norm() < 0.0001 ) return it->site;
			}
		}
		std::cerr << "ERROR : site not found" << std::endl;
		exit(1);

		return not_found;
	}

	// std::shared_ptr<LinkedSite> getLinkedSite(){return linked_site;}


	const int getSiteNum() const { return this->site_num; };
	const int getLatticeNum();
	const Eigen::Vector3d& getCoordinate() const { return this->coordinate; };

	bool isCoordinate( const Eigen::Vector3d& obj ){
		for( int i=0; i<3; ++i ){
			if( std::abs( this->coordinate[i]-obj[i] ) > 0.00001 ) {
				return false;
			}
		}
		return true;
	}
};


#endif
