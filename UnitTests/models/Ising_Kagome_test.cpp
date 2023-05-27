#ifdef Xcode
#include "../Catch2_amalgamated/catch_amalgamated.hpp"
#else
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#endif

#include "../../models/Kagome/Ising_Kagome.hpp"

int KagomeLattice_N_sites(int L0, int L1) {
	const int dim = 2;
	std::vector<int> L(dim);
	L[0] = L0;
	L[1] = L1;
	KagomeLattice lattice(L);
	return lattice._N_sites();
}

int KagomeLattice_NN(int L0, int L1, int site_idx, int bond_idx) {
	const int dim = 2;
	std::vector<int> L(dim);
	L[0] = L0;
	L[1] = L1;
	KagomeLattice lattice(L);
	return lattice._NN_of_Site(site_idx, bond_idx);
}

TEST_CASE("KagomeLattice_0"){
	REQUIRE(KagomeLattice_N_sites(4, 4) == 48);
	REQUIRE(KagomeLattice_NN(4, 4, 0, 0) == 10);
	REQUIRE(KagomeLattice_NN(4, 4, 0, 1) == 38);
	REQUIRE(KagomeLattice_NN(4, 4, 0, 2) == 1);
	REQUIRE(KagomeLattice_NN(4, 4, 0, 3) == 2);
}

TEST_CASE("KagomeLattice_1"){
	REQUIRE(KagomeLattice_N_sites(4, 4) == 48);
	REQUIRE(KagomeLattice_NN(4, 4, 1, 0) == 41);
	REQUIRE(KagomeLattice_NN(4, 4, 1, 1) == 3);
	REQUIRE(KagomeLattice_NN(4, 4, 1, 2) == 2);
	REQUIRE(KagomeLattice_NN(4, 4, 1, 3) == 0);
}

TEST_CASE("KagomeLattice_2"){
	REQUIRE(KagomeLattice_N_sites(4, 4) == 48);
	REQUIRE(KagomeLattice_NN(4, 4, 5, 0) == 15);
	REQUIRE(KagomeLattice_NN(4, 4, 5, 1) == 13);
	REQUIRE(KagomeLattice_NN(4, 4, 5, 2) == 3);
	REQUIRE(KagomeLattice_NN(4, 4, 5, 3) == 4);
}

TEST_CASE("KagomeLatticeIsingSystem(1)"){
	const int dim = 2;
	std::vector<int>L(dim);
	L[0] = 2;
	L[1] = 2;
	KagomeLatticeIsingSystem test_system(L);

	REQUIRE(test_system._n_spins()==12);
	REQUIRE(test_system.eval_energy()==-24.0);
	REQUIRE(test_system.eval_mz()==12.0);
}

TEST_CASE("KagomeLatticeIsingSystem(2)"){
	const int dim = 2;
	std::vector<int>L(dim);
	L[0] = 2;
	L[1] = 2;
	KagomeLatticeIsingSystem test_system(L);
	test_system.set_spin(2, -1);
	test_system.set_spin(6, -1);
	test_system.set_spin(7, -1);
	test_system.set_spin(9, -1);
	test_system.set_spin(11, -1);
	
	REQUIRE(test_system._n_spins() == 12);
	REQUIRE(test_system.eval_energy() == 0.0);
	REQUIRE(test_system.eval_mz() == 2.0);
}

TEST_CASE("KagomeLatticeIsingSystem Integer Rep. of States (1)") {
	const int n_spins = 12;
	IsingSystem test_system(n_spins);
	test_system.set_state_by_code(2363);
	REQUIRE(test_system._sz(0) == 1);
	REQUIRE(test_system._sz(1) == 1);
	REQUIRE(test_system._sz(2) == -1);
	REQUIRE(test_system._sz(3) == 1);
	REQUIRE(test_system._sz(4) == 1);
	REQUIRE(test_system._sz(5) == 1);
	REQUIRE(test_system._sz(6) == -1);
	REQUIRE(test_system._sz(7) == -1);
	REQUIRE(test_system._sz(8) == 1);
	REQUIRE(test_system._sz(9) == -1);
	REQUIRE(test_system._sz(10) == -1);
	REQUIRE(test_system._sz(11) == 1);
}

TEST_CASE("KagomeLatticeIsingSystem Integer Rep. of States (2)") {
	const int dim = 2;
	std::vector<int> L(dim);
	L[0] = 2;
	L[1] = 1;
	KagomeLatticeIsingSystem test_system(L);
	REQUIRE(test_system._maxrep_state() == 63);
	
	test_system.set_state_by_code(0);
	REQUIRE(test_system.eval_mz() == -6.0);
	REQUIRE(test_system.eval_energy() == -12.0);
	
	test_system.set_state_by_code(3);
	REQUIRE(test_system.eval_mz() == -2.0);
	REQUIRE(test_system.eval_energy() == 0.0);

	test_system.set_state_by_code(2);
	REQUIRE(test_system.eval_mz() == -4.0);
	REQUIRE(test_system.eval_energy() == -4.0);

	test_system.set_state_by_code(63);
	REQUIRE(test_system.eval_mz() == 6.0);
	REQUIRE(test_system.eval_energy() == -12.0);
}

TEST_CASE("KagomeLatticeIsingSystem Integer Rep. of States (3)") {
	const int dim = 2;
	std::vector<int> L(dim);
	L[0] = 2;
	L[1] = 2;
	KagomeLatticeIsingSystem test_system(L);
	REQUIRE(test_system._maxrep_state() == 4095);
	
	test_system.set_state_by_code(0);
	REQUIRE(test_system.eval_mz() == -12.0);
	REQUIRE(test_system.eval_energy() == -24.0);
	
	test_system.set_state_by_code(4095);
	REQUIRE(test_system.eval_mz() == 12.0);
	REQUIRE(test_system.eval_energy() == -24.0);
}


TEST_CASE("KagomeLatticeIsingSystem Observables Brute-Force Counting") {
	const int dim = 2;
	std::vector<int> L(dim);
	L[0] = 1;
	L[1] = 1;
	const int n_spins = KagomeLatticeIsingSystem::eval_n_spins(L);
	const std::vector<double> beta{0, 0.3, 1.1};
	KagomeLatticeDataBundle data_bundle(n_spins, beta);
	KagomeLatticeIsingSystem test_system(L, data_bundle);
	test_system.exact_count();
	
	const double Jabs = std::fabs(test_system._J());
	const double mz = 0.0;
	for ( unsigned int beta_idx = 0; beta_idx < beta.size(); beta_idx++ ) {
		const double b = beta[beta_idx];
		const double mz_sq = (3.0 + 9.0 * std::exp(8.0 * b * Jabs)) / (3.0 + std::exp(8.0 * b * Jabs));
		const double energy = 6.0 * (1.0 - std::exp(8.0 * b * Jabs)) / (3.0 + std::exp(8.0 * b * Jabs));
		const double ene_sq = (12.0 + 36.0 * std::exp(8.0 * b * Jabs)) / (3.0 + std::exp(8.0 * b * Jabs));
		
		REQUIRE_THAT(data_bundle._get_exact_magz(beta_idx), Catch::Matchers::WithinRel(mz, 1e-12));
		REQUIRE_THAT(data_bundle._get_exact_magz_sq(beta_idx), Catch::Matchers::WithinRel(mz_sq, 1e-12));
		REQUIRE_THAT(data_bundle._get_exact_energy(beta_idx), Catch::Matchers::WithinRel(energy, 1e-12));
		REQUIRE_THAT(data_bundle._get_exact_C(beta_idx), Catch::Matchers::WithinRel(b * b * (ene_sq - energy * energy), 1e-12));
	}
}

bool KagomeLatticeIsingSystem_MC_stderr_validation(const int L_spec, const double b_spec, const int n_trials, const double dev_tol = 3.0) {
	const double CI = 0.95;
	const int n_bins = 20;
	const double t_inv_0975_19 = 2.093;
	const double stddev_binomial = sqrt(CI * (1.0 - CI) / n_trials);

	const int n_samples_per_bin = 4000;
	const int mcs_thermalization = 1000;
	const int mcs_interval_btwn_bins = 100;
	
	const int dim = 2;
	std::vector<int> L(dim);
	L[0] = L_spec;
	L[1] = L[0];
	const int n_sites = KagomeLatticeIsingSystem::eval_n_spins(L);
	std::cout << "System size in KagomeLatticeIsingSystem_MC_stderr_validation(...) ... " << L[0] << " x " << L[1] << std::endl;
	
	std::vector<double> beta = { b_spec };
	const int seed = 745099137;
	RandomNumberGenerator mtwist(seed, n_sites);
	IsingSystem::setup_RNGen(mtwist);
	
	int n_in_CI_energy = 0;
	int n_in_CI_C = 0;
	int n_in_CI_magz_sq = 0;
	double exact_energy = 0;
	double exact_C = 0;
	double exact_magz_sq = 0;
	const int beta_idx = 0;
	for (int i_trial = 0; i_trial < n_trials; i_trial++) {
		KagomeLatticeDataBundle data_bundle(n_sites, beta, n_bins, n_samples_per_bin);
		KagomeLatticeIsingSystem test_system(L, data_bundle, mcs_thermalization, mcs_interval_btwn_bins);
		if (i_trial==0) {
			test_system.exact_count();
			exact_energy = data_bundle._get_exact_energy(beta_idx);
			exact_C = data_bundle._get_exact_C(beta_idx);
			exact_magz_sq = data_bundle._get_exact_magz_sq(beta_idx);
		}
		test_system.run_MC(false);
		
		const double dlt_energy = std::fabs( data_bundle._get_obs_energy(beta_idx) - exact_energy );
		if ( dlt_energy <= t_inv_0975_19 * data_bundle._get_stderr_energy(beta_idx) ) n_in_CI_energy++;
		const double dlt_C = std::fabs( data_bundle._get_obs_C(beta_idx) - exact_C );
		if ( dlt_C <= t_inv_0975_19 * data_bundle._get_stderr_C(beta_idx) ) n_in_CI_C++;
		const double dlt_magz_sq = std::fabs( data_bundle._get_obs_magz_sq(beta_idx) - exact_magz_sq );
		if ( dlt_magz_sq <= t_inv_0975_19 * data_bundle._get_stderr_magz_sq(beta_idx) ) n_in_CI_magz_sq++;

	}
	double portion_in_CI_energy = static_cast<double>(n_in_CI_energy) / n_trials;
	std::cout << "T = " << 1.0 / b_spec << " : Energy : Portion within CI (95%) = " << portion_in_CI_energy << " (stddev for the binomial dist. = " << stddev_binomial << "; tol. parameter = " << dev_tol << ")";
	bool check_energy = (std::fabs(portion_in_CI_energy - CI) / stddev_binomial < dev_tol);
	if (check_energy) {
		std::cout << " ... OK" << std::endl;
	} else {
		std::cout << " (*)" << std::endl;
	}
	
	double portion_in_CI_C = static_cast<double>(n_in_CI_C) / n_trials;
	std::cout << "T = " << 1.0 / b_spec << " : C : Portion within CI (95%) = " << portion_in_CI_C << " (stddev for the binomial dist. = " << stddev_binomial << "; tol. parameter = " << dev_tol << ")";
	bool check_C = (std::fabs(portion_in_CI_C - CI) / stddev_binomial < dev_tol);
	if (check_C) {
		std::cout << " ... OK" << std::endl;
	} else {
		std::cout << " (*)" << std::endl;
	}

	double portion_in_CI_magz_sq = static_cast<double>(n_in_CI_magz_sq) / n_trials;
	std::cout << "T = " << 1.0 / b_spec << " : M^2 : Portion within CI (95%) = " << portion_in_CI_magz_sq << " (stddev for the binomial dist. = " << stddev_binomial << "; tol. parameter = " << dev_tol << ")";
	bool check_magz_sq = (std::fabs(portion_in_CI_magz_sq - CI) / stddev_binomial < dev_tol);
	if (check_magz_sq) {
		std::cout << " ... OK" << std::endl;
	} else {
		std::cout << " (*)" << std::endl;
	}

	return check_energy && check_magz_sq && check_C;
};

TEST_CASE("KagomeLatticeIsingSystem : MC validation") {
	const int L = 2;
	const int n_trials = 50;

	std::vector<double> beta = { 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1 };
	for (std::vector<double>::iterator itr = beta.begin(); itr != beta.end(); itr++) {
		REQUIRE(KagomeLatticeIsingSystem_MC_stderr_validation(L, *itr, n_trials) == true );
	}
}