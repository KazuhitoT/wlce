#ifndef __MYINDEX_HPP
#define __MYINDEX_HPP

#include <vector>
#include <cmath>
#include <memory>
#include "./site.hpp"

class Site;

struct PairIndex{
	std::vector<int> lattice_index;
	double distance;

	PairIndex(const std::shared_ptr<Site>&, const std::shared_ptr<Site>&);          /* using non-casted distance temporary*/
	PairIndex(const std::shared_ptr<Site>&, const std::shared_ptr<Site>&, double);  /* using casted distance */

	bool operator==(const PairIndex& rhs) const;
	bool operator!=(const PairIndex& rhs) const;
	bool operator<(const PairIndex& rhs) const {
		return this->distance < rhs.distance;
	}
};

inline bool PairIndex::operator==(const PairIndex& rhs) const {
	return this->lattice_index == rhs.lattice_index and std::fabs(this->distance-rhs.distance) < 0.00001;
}

inline bool PairIndex::operator!=(const PairIndex& rhs) const {
    return !(this->operator==(rhs));
}

class hash_pairindex {
public:
  double operator()(const PairIndex &x) const {
    const int C = 997;  // prime_number
    double t = 0;
    for (int i = 0; i<2; ++i) {
        t = t * C + x.lattice_index[i];
    }
		t = t*C + x.distance;
    return t;
  }
};


struct TripletIndex{
	std::vector<PairIndex> vec_pairindex;

	TripletIndex(const std::vector<PairIndex>& _r) : vec_pairindex(_r){};

	std::vector<double> getVecDistance() const;
	bool operator==(const TripletIndex& rhs) const;
	bool operator!=(const TripletIndex& rhs) const;
};

inline std::vector<double> TripletIndex::getVecDistance() const {
	std::vector<PairIndex> vec_pairindex;
	for( const auto& pairindex : vec_pairindex ){
		vec_pairindex.push_back(pairindex);
	}
	std::sort(vec_pairindex.begin(), vec_pairindex.end());
	vec_pairindex.erase(std::unique(vec_pairindex.begin(), vec_pairindex.end()), vec_pairindex.end());
	std::vector<double> vec_distance;
	for(const auto& pairindex : vec_pairindex ){
		vec_distance.push_back(pairindex.distance);
	}
	return vec_distance;
}

inline bool TripletIndex::operator==(const TripletIndex& rhs) const {
    return this->vec_pairindex[0] == rhs.vec_pairindex[0]
				and this->vec_pairindex[1] == rhs.vec_pairindex[1]
				and this->vec_pairindex[2] == rhs.vec_pairindex[2];
}

inline bool TripletIndex::operator!=(const TripletIndex& rhs) const {
    return !(this->operator==(rhs));
}

class hash_tripletindex {
public:
  double operator()(const TripletIndex &x) const {
    const int C = 997;  // prime_number
    double t = 0;
		for( const auto pairindex : x.vec_pairindex ){
			hash_pairindex hashp;
			t = t*C + hashp(pairindex);
		}
    return t;
  }
};

struct QuadrupletIndex{
	std::vector<TripletIndex> vec_tripletindex;

	QuadrupletIndex(const std::vector<TripletIndex>& _r) : vec_tripletindex(_r){};

	std::vector<double> getVecDistance() const;
	bool operator==(const QuadrupletIndex& rhs) const;
	bool operator!=(const QuadrupletIndex& rhs) const;
};

inline bool QuadrupletIndex::operator==(const QuadrupletIndex& rhs) const {
    return this->vec_tripletindex[0] == rhs.vec_tripletindex[0]
				and this->vec_tripletindex[1] == rhs.vec_tripletindex[1]
				and this->vec_tripletindex[2] == rhs.vec_tripletindex[2]
				and this->vec_tripletindex[3] == rhs.vec_tripletindex[3];
}

inline bool QuadrupletIndex::operator!=(const QuadrupletIndex& rhs) const {
    return !(this->operator==(rhs));
}

inline std::vector<double> QuadrupletIndex::getVecDistance() const {
	std::vector<PairIndex> vec_pairindex;
	for( const auto& tripletindex : vec_tripletindex ){
		for( const auto& pairindex : tripletindex.vec_pairindex ){
			vec_pairindex.push_back(pairindex);
		}
	}
	std::sort(vec_pairindex.begin(), vec_pairindex.end());
	vec_pairindex.erase(std::unique(vec_pairindex.begin(), vec_pairindex.end()), vec_pairindex.end());
	std::vector<double> vec_distance;
	for(const auto& pairindex : vec_pairindex ){
		vec_distance.push_back(pairindex.distance);
	}
	return vec_distance;
}

class hash_quadrupletindex {
public:
  double operator()(const QuadrupletIndex &x) const {
    const int C = 997;  // prime_number
    double t = 0;
		for( const auto tripletindex : x.vec_tripletindex ){
			hash_tripletindex hasht;
      t = t * C + hasht(tripletindex);
		}
    return t;
  }
};


#endif
