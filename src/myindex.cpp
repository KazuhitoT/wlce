#include "./myindex.hpp"

PairIndex::PairIndex(const std::shared_ptr<Site>& lhs, const std::shared_ptr<Site>& rhs, double _distance){
	this->lattice_index = {lhs->getLatticeNum(), rhs->getLatticeNum()};
	this->distance = _distance;
}
