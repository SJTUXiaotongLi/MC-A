#include "Ising_Honeycomb.hpp"

int main(int argc, char **argv) {
	std::cout << "# Exact evaluation (Honeycomb lattice)" << std::endl;
	if (argc <= 1) {
		std::cerr << "We need a parameter file (see sample_parameters/...)" << std::endl;
		exit(1);
	}
	HoneycombLatticeParameterBundle_exact parameter(argv[1], true);

	/* Honeycomb lattice system */
	const int dim = 2;
	std::vector<int> L(dim);
	L[0] = parameter._L0();
	L[1] = parameter._L1();
	const int n_sites = HoneycombLatticeIsingSystem::eval_n_spins(L);
	
	/* list of temperatures */
	const double Tmin = parameter._Tmin();
	const double Tmax = parameter._Tmax();
	const int nT = parameter._N_TemperaturePoints();
	const double Tdlt = (Tmax - Tmin) / (nT - 1);
	std::vector<double> beta(nT);
	for (int i = 0; i < nT; i++) beta[i] = 1.0 / (Tmin + i * Tdlt);
	
	/* simulation */
	HoneycombLatticeDataBundle data_bundle(n_sites, beta);
	HoneycombLatticeIsingSystem system(L, data_bundle);
	system.exact_count();
	
	data_bundle.output_legends_exact(system._system_size());
	data_bundle.output_data_exact_all();
	
	return 0;
}
