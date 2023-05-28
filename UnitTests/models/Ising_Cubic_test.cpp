#ifdef Xcode
#include "../Catch2_amalgamated/catch_amalgamated.hpp"
#else
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#endif

#include "../../models/Cubic/Ising_Cubic.hpp"

TEST_CASE("CubicLatticeIsingSystem 1") {
	const int dim = 3;
	std::vector<int> L(dim);
	L[0] = 3;
	L[1] = 3;
	L[2] = 3;
	CubicLatticeIsingSystem test_system(L);
	test_system.set_spin(2, -1);
	test_system.set_spin(6, -1);
	test_system.set_spin(7, -1);
	test_system.set_spin(9, -1);
	test_system.set_spin(11, -1);
	
	REQUIRE(test_system._n_spins() == 27);
	REQUIRE(test_system.eval_energy() == -33.0);
	REQUIRE(test_system.eval_mz() == 17.0);
}

TEST_CASE("CubicLatticeIsingSystem 2") {
	const int dim = 3;
	std::vector<int> L(dim);
	L[0] = 3;
	L[1] = 3;
	L[2] = 3;
	CubicLatticeIsingSystem test_system(L);
	test_system.set_state_by_code(1339);
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
	REQUIRE(test_system._sz(10) == 1);
	REQUIRE(test_system._sz(11) == -1);
	REQUIRE(test_system.eval_energy() == -29.0);
	REQUIRE(test_system.eval_mz() == -13.0);
}

// TEST_CASE("CubicLatticeIsingSystem Observables Brute-Force Counting") {
// 	const int dim = 3;
// 	std::vector<int> L(dim);
// 	L[0] = 2;
// 	L[1] = 2;
// 	L[2] = 2;
// 	const int n_sites = CubicLatticeIsingSystem::eval_n_spins(L);
// 	const std::vector<double> beta{ 0, 0.3, 1.1, 10.1, 100.1 };
// 	CubicLatticeDataBundle data_bundle(n_sites, beta);
// 	CubicLatticeIsingSystem test_system(L, data_bundle);
// 	REQUIRE(data_bundle._get_obs_energy(0) == 0);
// 	REQUIRE(data_bundle._get_obs_magz(0) == 0);
// 	test_system.exact_count();
	
// 	const double Jabs = std::fabs(test_system._J());
// 	const double mz = 0.0;
// 	for ( unsigned int beta_idx = 0; beta_idx < beta.size(); beta_idx++ ) {
// 		const double b = beta[beta_idx];
// 		const double e4b=std::exp(4.0*b*Jabs);
// 		const double z=2.0*pow(e4b,-6.0)*pow(1+e4b,4.0)*pow(1+e4b*e4b,2.0)*(pow(1-e4b,4.0)+2.0*e4b*e4b);
// 		const double mz_sq = 64*pow(e4b,-3.0)*(1+3*e4b*e4b+8*e4b*e4b*e4b+3*pow(e4b,4.0)+6*pow(e4b,5.0)+9*pow(e4b,6.0)+2*pow(e4b,9.0))/z;

// 		const double energy = -24.0*(-1.0+(3.0*e4b-5.0*e4b*e4b+3.0*e4b*e4b*e4b)*(1-e4b*e4b*e4b)+pow(e4b,7))/((1+e4b)*(1+e4b*e4b)*(1.0-4.0*e4b+8.0*e4b*e4b-4.0*pow(e4b,3)+pow(e4b,4)));
// 		const double ene_sq = 384.0 * (3*pow(e4b,-6.0)+6*pow(e4b,-2.0)+2/e4b+2*e4b+5*e4b*e4b+6*e4b*e4b*e4b+3*pow(e4b,6.0))/z;
			
// 		REQUIRE_THAT(data_bundle._get_exact_magz(beta_idx), Catch::Matchers::WithinRel(mz, 1e-2));
// 		REQUIRE_THAT(data_bundle._get_exact_magz_sq(beta_idx), Catch::Matchers::WithinRel(mz_sq, 1e-2));
// 		REQUIRE_THAT(data_bundle._get_exact_energy(beta_idx), Catch::Matchers::WithinRel(energy, 1e-2));
// 		REQUIRE_THAT(data_bundle._get_exact_C(beta_idx), Catch::Matchers::WithinRel(b * b * (ene_sq - energy * energy), 1e-2));
// 	}
// }

bool CubicLatticeIsingSystem_MC_stderr_validation(const int L_spec, const double b_spec, const int n_trials, const double dev_tol = 3.0) {
	const double CI = 0.95;
	const int n_bins = 20;
	const double t_inv_0975_19 = 2.093;
	const double stddev_binomial = sqrt(CI * (1.0 - CI) / n_trials);

	const int n_samples_per_bin = 10000;
	const int mcs_thermalization = 1000;
	const int mcs_interval_btwn_bins = 100;
	
	const int dim = 3;
	std::vector<int> L(dim);
	L[0] = L_spec;
	L[1] = L[0];
	L[2] = L[0];
	const int n_sites = CubicLatticeIsingSystem::eval_n_spins(L);
	std::cout << "System size in CubicLatticeIsingSystem_MC_stderr_validation(...) ... " << L[0] << " x " << L[1] << " x " << L[2] <<std::endl;
	
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
		CubicLatticeDataBundle data_bundle(n_sites, beta, n_bins, n_samples_per_bin);
		CubicLatticeIsingSystem test_system(L, data_bundle, mcs_thermalization, mcs_interval_btwn_bins);
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

TEST_CASE("CubicLatticeIsingSystem : MC validation") {
	const int L = 2;
	const int n_trials = 200;

	std::vector<double> beta = { 0.1, 0.2, 0.4, 0.6, 0.8};
	for (std::vector<double>::iterator itr = beta.begin(); itr != beta.end(); itr++) {
		REQUIRE(CubicLatticeIsingSystem_MC_stderr_validation(L, *itr, n_trials) == true );
	}
}
