#include <iomanip>
#include <iterator>
#include <numeric>
#include <random>
#include <chrono>
#include <mpi.h>

#include "../src/parser.hpp"
#include "../src/input.hpp"
#include "../src/metroconf.hpp"

constexpr double kb = 0.0000861733; // [eV/K]

double metropolis_like_step(Metroconf& conf, double p){
	double norm_before = conf.calcCorrelationFunctionNorm(p);

	conf.setNewConf();
	double norm_after  = conf.calcCorrelationFunctionNorm(p);

	// std::cout << norm_after << " " << norm_before << std::endl;

	double b = exp( -(norm_after-norm_before) );
	if(b>=1.0 or b>conf.RandReal()){
	} else {
		conf.Memento();
	}
	return conf.getTotalEnergy();
}


int main(int argc, char* argv[]){

	MPI::Init(argc, argv);

	int rank = MPI::COMM_WORLD.Get_rank();
	int size = MPI::COMM_WORLD.Get_size();

	std::vector<double> vec_p;
	int mcstep, samplestep, exstep;
	std::vector<double> p_input;

	std::function<int ()> rnd_rank = std::bind( std::uniform_int_distribution<int>(0, size-1), std::mt19937());
	std::function<int ()> rnd_real = std::bind( std::uniform_real_distribution<double>(0.0, 1.0), std::mt19937());

	std::shared_ptr<Input> in(new Input("remc-gs-search.ini"));

	in->setData("MCSTEP",      mcstep, true);
	in->setData("SAMPLESTEP",  samplestep, true);
	in->setData("EXSTEP",      exstep, true);
	in->setData("P", p_input, true);


	if ( p_input.size() == 2 ) {
		double tdelta = fabs(p_input[1] - p_input[0])/(double)(size-1);
		if( p_input[0] > p_input[1] ){ // 900 100 100
			for( int i=0; i<size; ++i){
				vec_p.push_back( p_input[1]+tdelta*i  );
			}
		} else if( p_input[0] < p_input[1] ){  // 100 100 900
			for( int i=0; i<size; ++i){
				vec_p.push_back( p_input[0]+tdelta*i  );
			}
		}
	} else {
		std::cerr << " ERROR : TEMPERATURE in remc.ini must be 2 arguments." << std::endl;
		exit(1);
	}

	const ParseLabels label("./labels.in");
	const ParseEcicar ecicar("./ecicar");
	const ParseMultiplicityIn multiplicity_in("./multiplicity.in", ecicar.getIndex());
	const ParseClusterIn  cluster_in("./clusters.in", ecicar.getIndex(), multiplicity_in.getMultiplicityIn());

	Metroconf PoscarSpin("./poscar.spin", in, label.getLabels(), cluster_in.getCluster(), ecicar.getEci(), nullptr, nullptr);

	const int N = label.getLabels()->size();

	MPI::COMM_WORLD.Barrier();

	std::ofstream ofs_trace("trace-"+std::to_string(rank)+".out");
	ofs_trace.setf(std::ios_base::fixed, std::ios_base::floatfield);

	for(int i=1; i<=mcstep; ++i){
		for(int j=0; j<N; ++j){
			metropolis_like_step(PoscarSpin, vec_p[rank]);
		}
		if( rank == 0 ) std::cout << i << " step done." << std::endl;

		ofs_trace << i << " ";
		for(const auto& corrs : PoscarSpin.getCorrelationFunctions() )
			ofs_trace << corrs[0] << " ";
		ofs_trace << std::endl;

		if( i%exstep == 0 ){
			int swp_rank1 = rnd_rank();
			int swp_rank2 = swp_rank1;
			if( swp_rank1==0 ) swp_rank2 = 1;
			else if( swp_rank1==size-1 ) swp_rank2 = swp_rank1-1;
			else{
				if( rnd_real()<=0.5 ) swp_rank2 = swp_rank1-1;
				else swp_rank2 = swp_rank1+1;
			}

			// if( rank == 0 ) std::cout << "+++ " << vec_p[swp_rank1] << " and " << vec_p[swp_rank2] << " was exchanged." << std::endl;
			std::swap(vec_p[swp_rank1], vec_p[swp_rank2]);
			MPI_Bcast(&vec_p[0], vec_p.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI::COMM_WORLD.Barrier();
		}

	}

	ofs_trace.close();

	std::vector<double> vec_ene;
	for(int i=1; i<=samplestep; ++i){
		for(int j=0; j<N; ++j){
			double e = metropolis_like_step(PoscarSpin, vec_p[rank]);
			vec_ene.push_back(e);
		}
		if( rank == 0 ) std::cout << "SAMPLE STEP " << i << " step done." << std::endl;
	}

	double average  = std::accumulate(vec_ene.begin(), vec_ene.end(), 0.0 ) / double(vec_ene.size());
	double variance = std::accumulate(vec_ene.begin(), vec_ene.end(), 0.0,
		[average]( double sum, double val ){ return sum+(val-average)*(val-average);} ) / double(vec_ene.size());

	MPI::COMM_WORLD.Barrier();

	if (rank == 0){
		std::ofstream ofs("remc-gs-search.out");
		ofs.setf(std::ios_base::fixed, std::ios_base::floatfield);

		ofs << vec_p[rank] << " " << average << " " << variance << std::endl;

		double recv_ave, recv_var;
		for(int i=1; i<size; ++i){
			MPI::Status status;
			MPI::COMM_WORLD.Recv(&recv_ave, 1, MPI::DOUBLE, i, 0, status);
			MPI::COMM_WORLD.Recv(&recv_var, 1, MPI::DOUBLE, i, 0, status);
			std::cout << vec_p[i] << " " << recv_ave << " " << recv_var << std::endl;
			ofs << vec_p[i] << " " << recv_ave << " " << recv_var << std::endl;
		}

		ofs.close();
	} else {
		MPI::Status status;
		MPI::COMM_WORLD.Send(&average,  1, MPI::DOUBLE, 0, 0);
		MPI::COMM_WORLD.Send(&variance, 1, MPI::DOUBLE, 0, 0);
	}


	MPI::Finalize();
	return 0;
}
