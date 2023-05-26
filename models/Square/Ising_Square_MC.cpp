#include <random>
#include <iostream>
#include "Ising_Square.hpp"


#include<fstream>

int main(int argc, char **argv)
{
	std::cout << "# MC simulation (square lattice)" << std::endl;
	if (argc <= 1) {
		std::cerr << "We need a parameter file (see sample_parameters/...)" << std::endl;
		exit(1);
	}
	SquareLatticeParameterBundle_MC parameter(argv[1], true);
	
	/* Square lattice system */
	const int dim = 2;
	std::vector<int> L(dim);
	L[0] = parameter._L0();
	L[1] = parameter._L1();
	const int n_sites = SquareLatticeIsingSystem::eval_n_spins(L);

	/* list of temperatures */
	const double Tmin = parameter._Tmin();
	const double Tmax = parameter._Tmax();
	const int nT = parameter._N_TemperaturePoints();
	const double Tdlt = (Tmax - Tmin) / (nT - 1);
	std::vector<double> beta(nT);
	for (int i = 0; i < nT; i++) beta[i] = 1.0 / (Tmin + i * Tdlt);
	
	/* MC parameters */
	const int n_bins = parameter._n_bins();
	const int mcs_thermalization = parameter._mcs_thermalization();
	const int n_samples_per_bin = parameter._n_samples_per_bin();
	const int mcs_interval_btwn_bins = parameter._mcs_interval_btwn_bins();
	
	/* simulation */
	const int seed = 745099137;
	RandomNumberGenerator mtwist(seed, n_sites);
	IsingSystem::setup_RNGen(mtwist);
	SquareLatticeDataBundle data_bundle(n_sites, beta, n_bins, n_samples_per_bin);
	SquareLatticeIsingSystem system(L, data_bundle, mcs_thermalization, mcs_interval_btwn_bins);
	
	//insert
 	std::ofstream file("Square_MC.txt");  // 打开文件
    std::streambuf* backup = std::cout.rdbuf();  // 备份 cout 缓存
    std::cout.rdbuf(file.rdbuf());  // 将 cout 缓存指向文件

	data_bundle.output_legends_MC(system._system_size());
	system.run_MC();
	
	return 0;
}
