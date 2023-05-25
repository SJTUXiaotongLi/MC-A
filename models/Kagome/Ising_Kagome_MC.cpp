#include <random>
#include <iostream>
#include "Ising_Kagome.hpp"

int main(int argc, char **argv)
{
	std::cout << "# MC simulation (Kagome lattice)" << std::endl;
	if (argc <= 1) {
		std::cerr << "We need a parameter file (see sample_parameters/...)" << std::endl;
		exit(1);
	}
	KagomeLatticeParameterBundle_MC parameter(argv[1], true);
	
	/* Kagome lattice system */
	const int dim = 2;
	std::vector<int> L(dim);
	L[0] = parameter._L0();
	L[1] = parameter._L1();
	const int n_sites = KagomeLatticeIsingSystem::eval_n_spins(L);

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
	KagomeLatticeDataBundle data_bundle(n_sites, beta, n_bins, n_samples_per_bin);
	KagomeLatticeIsingSystem system(L, data_bundle, mcs_thermalization, mcs_interval_btwn_bins);
	
	data_bundle.output_legends_MC(system._system_size());
	system.run_MC();
	
	return 0;
}
