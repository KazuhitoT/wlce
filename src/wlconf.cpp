#include "./wlconf.hpp"

WLconf::WLconf(char* filename,
	std::shared_ptr<Input> _in,
	std::shared_ptr<allclusters> _pall_clusters,
	std::map<int /*index*/ , double /*eci*/> ecicar,
	std::shared_ptr<basisfunc>   _pbasis_functions,
	std::shared_ptr<indexorders> _pindex_orders,
	bool is_exchange
):Conf2corr(filename, _in, _pall_clusters, _pbasis_functions, _pindex_orders){

	_in->setData("EMIN", emin, true);
	_in->setData("EMAX", emax, true);
	_in->setData("EMAX", bin, true);
	_in->setData("CHEMIPOT", chemical_potential);

}
