#include "./parsePoscar.hpp"

void outputPoscar(double lattice[3][3], double position[][3], int N, std::string prefix){
	std::ofstream ofs("poscar."+prefix);
	ofs << "POSCAR" << std::endl;
	ofs << "1.0" << std::endl;
	ofs.setf(std::ios::right, std::ios::adjustfield);
	ofs.setf(std::ios_base::fixed, std::ios_base::floatfield);
	ofs << std::setprecision(15);
	ofs << std::setw(20);

	for(int i=0; i<3; ++i){
		for(int j=0; j<3; ++j){
			ofs << lattice[i][j] << " ";
		}
		ofs << std::endl;
	}
	ofs << N << std::endl;
	ofs << "Direct" << std::endl;
	for(int i=0; i<N; ++i){
		for(int j=0; j<3; ++j){
			ofs << position[i][j] << " ";
		}
		ofs << std::endl;
	}
};


ParsePoscar::ParsePoscar(const char *filename){

  std::ifstream input(filename);

	if (input.fail()){
		std::cerr << "Error: Could not open " << filename << "\n";
		exit(1);
	}
	std::getline( input, comment );
	double unit;
	input >> unit;

	double x,y,z;
	for (int i = 0; i < 3; ++i){
		input >> x >> y >> z;
		axis.col(i) << x ,y, z;
	}
	axis *= unit;

	std::string del1;
	std::getline (input, del1);

	for (int i = 0; i < 20; ++i){
		int a = 0;
		input >> a;
		if ( a == 0){
			break;
		}
		else if ( a != 0){
			num_atoms.push_back(a);
		}
	}

	input.close();

	std::ifstream input2(filename);
	for (int i = 0; i < 6; ++i){
		std::string del2;
		std::getline( input2, del2);
	}
	std::getline(input2, coordinate_type);

	int num_all_atoms = std::accumulate(num_atoms.begin(), num_atoms.end(), 0);

	for( int i=0; i<num_all_atoms; ++i ){
		input2 >> x >> y >> z;
		std::pair<int, Eigen::Vector3d> atom(i, Eigen::Vector3d(x, y, z));
		atoms.push_back(atom);
	}

}

Eigen::Matrix3d ParsePoscar::getAxis(){
	return axis;
};

std::vector<int> ParsePoscar::getNumAtoms(){
	return num_atoms;
};

std::vector<std::pair<int, Eigen::Vector3d>> ParsePoscar::getAtoms(){
	return atoms;
};
