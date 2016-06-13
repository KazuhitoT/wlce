#include <vector>

class hash_vecd {
public:
  double operator()(const std::vector<double> &x) const {
    const int C = 997;
    double result = 0;
    for (int i=0; i<x.size(); ++i) {
        result = result * C + x[i];
    }
    return result;
  }
};

class hash_vec2d {
public:
  double operator()(const std::vector<std::vector<double>> &x) const {
    const int C = 997;
    double result = 0;
		for (int i=0; i<x.size(); ++i) {
			for (int j=0; j<x[i].size(); ++j) {
        result = result * C + x[i][j];
			}
    }
    return result;
  }
};
