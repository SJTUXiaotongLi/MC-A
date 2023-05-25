#ifdef Xcode
#include "../Catch2_amalgamated/catch_amalgamated.hpp"
#else
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#endif


#include "../../models/Honeycomb/Ising_Honeycomb.hpp"

TEST_CASE("HoneycombLatticeIsingSystem") {
    const int dim = 2;
    std::vector<int> L(dim);
    L[0] = 3;
    L[1] = 3;
    HoneycombLatticeIsingSystem test_system(L);
    test_system.set_spin(0, -1);  
    
    REQUIRE(test_system._n_spins() == 18);
    REQUIRE(test_system.eval_mz() == 16.0);
   //w
    REQUIRE(test_system.eval_energy() == -21.0);
}

bool HoneycombLatticeIsingSystem_MC_stderr_validation(const int L_spec, const double b_spec, const int n_trials, const double dev_tol = 3.0) {
    const double CI = 0.95;
    const int n_bins = 20;
    const double t_inv_0975_19 = 2.093;
    const double stddev_binomial = sqrt(CI * (1.0 - CI) / n_trials);

    const int n_samples_per_bin = 9000;
    const int mcs_thermalization = 1000;
    const int mcs_interval_btwn_bins = 100;
    
    const int dim = 2;
    std::vector<int> L(dim);
    L[0] = L_spec;
    L[1] = L[0];
    const int n_sites = HoneycombLatticeIsingSystem::eval_n_spins(L);
    std::cout << "System size in HoneycombLatticeIsingSystem_MC_stderr_validation(...) ... " << L[0] << " x " << L[1] << std::endl;
    
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
        HoneycombLatticeDataBundle data_bundle(n_sites, beta, n_bins, n_samples_per_bin);
        HoneycombLatticeIsingSystem test_system(L, data_bundle, mcs_thermalization, mcs_interval_btwn_bins);
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



//TEST_CASE("HoneycombLatticeIsingSystem Weiss field") {
//	const int dim = 2;
//	std::vector<int> L(dim);
//	L[0] = 4;
//	L[1] = 3;
//	HoneycombLatticeIsingSystem test_system(L);
//	test_system.set_state_by_code(1339);
//	/*
//	 + - + -
//	 + + - -
//	 + + - +
//	 */
//	REQUIRE_THAT(test_system.eval_weiss_at(0), Catch::Matchers::WithinAbs(4.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_weiss_at(1), Catch::Matchers::WithinAbs(0.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_weiss_at(2), Catch::Matchers::WithinAbs(2.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_weiss_at(3), Catch::Matchers::WithinAbs(-2.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_weiss_at(4), Catch::Matchers::WithinAbs(2.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_weiss_at(5), Catch::Matchers::WithinAbs(0.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_weiss_at(6), Catch::Matchers::WithinAbs(0.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_weiss_at(7), Catch::Matchers::WithinAbs(0.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_weiss_at(8), Catch::Matchers::WithinAbs(0.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_weiss_at(9), Catch::Matchers::WithinAbs(4.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_weiss_at(10), Catch::Matchers::WithinAbs(-4.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_weiss_at(11), Catch::Matchers::WithinAbs(2.0, 1e-15));
//}

//TEST_CASE("HoneycombLatticeIsingSystem Metropolis probability") {
//	const int dim = 2;
//	std::vector<int> L(dim);
//	L[0] = 4;
//	L[1] = 3;
//	std::vector<double> beta = { 0.5 };
//	DataBundle data_bundle(beta);
//	HoneycombLatticeIsingSystem test_system(L, &data_bundle);
//	test_system.set_state_by_code(1339);
//	/*
//	 + - + -
//	 + + - -
//	 + + - +
//	 */
//	REQUIRE_THAT(test_system.eval_weiss_at(0), Catch::Matchers::WithinAbs(4.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(0), Catch::Matchers::WithinAbs(0.01831563888, 1e-11));
//
//	REQUIRE_THAT(test_system.eval_weiss_at(1), Catch::Matchers::WithinAbs(0.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(1), Catch::Matchers::WithinAbs(1.00000000000, 1e-11));
//
//	REQUIRE_THAT(test_system.eval_weiss_at(2), Catch::Matchers::WithinAbs(2.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(2), Catch::Matchers::WithinAbs(7.38905609893, 1e-11));
//
//	REQUIRE_THAT(test_system.eval_weiss_at(3), Catch::Matchers::WithinAbs(-2.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(3), Catch::Matchers::WithinAbs(7.38905609893, 1e-11));
//
//	REQUIRE_THAT(test_system.eval_weiss_at(4), Catch::Matchers::WithinAbs(2.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(4), Catch::Matchers::WithinAbs(0.13533528323, 1e-11));
//
//	REQUIRE_THAT(test_system.eval_weiss_at(5), Catch::Matchers::WithinAbs(0.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(5), Catch::Matchers::WithinAbs(1.00000000000, 1e-11));
//
//	REQUIRE_THAT(test_system.eval_weiss_at(6), Catch::Matchers::WithinAbs(0.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(6), Catch::Matchers::WithinAbs(1.00000000000, 1e-11));
//
//	REQUIRE_THAT(test_system.eval_weiss_at(7), Catch::Matchers::WithinAbs(0.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(7), Catch::Matchers::WithinAbs(1.00000000000, 1e-11));
//
//	REQUIRE_THAT(test_system.eval_weiss_at(8), Catch::Matchers::WithinAbs(0.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(8), Catch::Matchers::WithinAbs(1.00000000000, 1e-11));
//
//	REQUIRE_THAT(test_system.eval_weiss_at(9), Catch::Matchers::WithinAbs(4.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(9), Catch::Matchers::WithinAbs(54.5981500331, 1e-10));
//
//	REQUIRE_THAT(test_system.eval_weiss_at(10), Catch::Matchers::WithinAbs(-4.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(10), Catch::Matchers::WithinAbs(54.5981500331, 1e-10));
//
//	REQUIRE_THAT(test_system.eval_weiss_at(11), Catch::Matchers::WithinAbs(2.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(11), Catch::Matchers::WithinAbs(7.38905609893, 1e-11));
//}

//TEST_CASE("HoneycombLatticeIsingSystem Metropolis update") {
//	const int dim = 2;
//	std::vector<int> L(dim);
//	L[0] = 4;
//	L[1] = 3;
//	std::vector<double> beta = { 0.5 };
//	DataBundle data_bundle(beta);
//	HoneycombLatticeIsingSystem test_system(L, &data_bundle);
//
//	const int seed = 745099137;
//	RandomNumberGenerator mtwist(seed);
//	test_system.setup_RNGen(&mtwist);
//	/* cf. https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution */
//	// 0.455093 0.137216 0.967235 0.771364 0.244447 0.565477 0.8711 0.441108 0.999349 0.340951 0.907204 0.945926
//
//	test_system.set_state_by_code(1339);
//	/*
//	 + - + -
//	 + + - -
//	 + + - +
//	 */
//	REQUIRE_THAT(test_system.eval_weiss_at(0), Catch::Matchers::WithinAbs(4.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(0), Catch::Matchers::WithinAbs(0.01831563888, 1e-11));
//	REQUIRE(test_system._sz(0) == 1);
//	test_system.spinflip_Metropolis(0);
//	REQUIRE(test_system._sz(0) == 1);
//
//	REQUIRE_THAT(test_system.eval_weiss_at(1), Catch::Matchers::WithinAbs(0.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(1), Catch::Matchers::WithinAbs(1.00000000000, 1e-11));
//	REQUIRE(test_system._sz(1) == 1);
//	test_system.spinflip_Metropolis(1);
//	REQUIRE(test_system._sz(1) == -1);
//	/*
//	 + - + -
//	 + + - -
//	 + - - +
//	 */
//
//	REQUIRE_THAT(test_system.eval_weiss_at(2), Catch::Matchers::WithinAbs(0.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(2), Catch::Matchers::WithinAbs(1.00000000000, 1e-11));
//	REQUIRE(test_system._sz(2) == -1);
//	test_system.spinflip_Metropolis(2);
//	REQUIRE(test_system._sz(2) == 1);
//	/*
//	 + - + -
//	 + + - -
//	 + - + +
//	 */
//
//	REQUIRE_THAT(test_system.eval_weiss_at(3), Catch::Matchers::WithinAbs(0.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(3), Catch::Matchers::WithinAbs(1.00000000000, 1e-11));
//	REQUIRE(test_system._sz(3) == 1);
//	test_system.spinflip_Metropolis(3);
//	REQUIRE(test_system._sz(3) == -1);
//	/*
//	 + - + -
//	 + + - -
//	 + - + -
//	 */
//
//	REQUIRE_THAT(test_system.eval_weiss_at(4), Catch::Matchers::WithinAbs(2.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(4), Catch::Matchers::WithinAbs(0.13533528323, 1e-11));
//	REQUIRE(test_system._sz(4) == 1);
//	test_system.spinflip_Metropolis(4);
//	REQUIRE(test_system._sz(4) == 1);
//
//	REQUIRE_THAT(test_system.eval_weiss_at(5), Catch::Matchers::WithinAbs(-2.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(5), Catch::Matchers::WithinAbs(7.38905609893, 1e-11));
//	REQUIRE(test_system._sz(5) == 1);
//	test_system.spinflip_Metropolis(5);
//	REQUIRE(test_system._sz(5) == -1);
//	/*
//	 + - + -
//	 + - - -
//	 + - + -
//	 */
//
//	REQUIRE_THAT(test_system.eval_weiss_at(6), Catch::Matchers::WithinAbs(0.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(6), Catch::Matchers::WithinAbs(1.00000000000, 1e-11));
//	REQUIRE(test_system._sz(6) == -1);
//	test_system.spinflip_Metropolis(6);
//	REQUIRE(test_system._sz(6) == 1);
//	/*
//	 + - + -
//	 + - + -
//	 + - + -
//	 */
//
//	REQUIRE_THAT(test_system.eval_weiss_at(7), Catch::Matchers::WithinAbs(0.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(7), Catch::Matchers::WithinAbs(1.00000000000, 1e-11));
//	REQUIRE(test_system._sz(7) == -1);
//	test_system.spinflip_Metropolis(7);
//	REQUIRE(test_system._sz(7) == 1);
//	/*
//	 + - + -
//	 + - + +
//	 + - + -
//	 */
//
//	REQUIRE_THAT(test_system.eval_weiss_at(8), Catch::Matchers::WithinAbs(0.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(8), Catch::Matchers::WithinAbs(1.00000000000, 1e-11));
//	REQUIRE(test_system._sz(8) == 1);
//	test_system.spinflip_Metropolis(8);
//	REQUIRE(test_system._sz(8) == -1);
//	/*
//	 - - + -
//	 + - + +
//	 + - + -
//	 */
//
//	REQUIRE_THAT(test_system.eval_weiss_at(9), Catch::Matchers::WithinAbs(-2.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(9), Catch::Matchers::WithinAbs(0.13533528323, 1e-10));
//	REQUIRE(test_system._sz(9) == -1);
//	test_system.spinflip_Metropolis(9);
//	REQUIRE(test_system._sz(9) == -1);
//
//	REQUIRE_THAT(test_system.eval_weiss_at(10), Catch::Matchers::WithinAbs(0.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(10), Catch::Matchers::WithinAbs(1.00000000000, 1e-11));
//	REQUIRE(test_system._sz(10) == 1);
//	test_system.spinflip_Metropolis(10);
//	REQUIRE(test_system._sz(10) == -1);
//	/*
//	 - - - -
//	 + - + +
//	 + - + -
//	 */
//
//	REQUIRE_THAT(test_system.eval_weiss_at(11), Catch::Matchers::WithinAbs(-2.0, 1e-15));
//	REQUIRE_THAT(test_system.eval_prob_spinflip_Metropolis(11), Catch::Matchers::WithinAbs(0.13533528323, 1e-11));
//	REQUIRE(test_system._sz(11) == -1);
//	test_system.spinflip_Metropolis(11);
//	REQUIRE(test_system._sz(11) == -1);
//}

//TEST_CASE("HoneycombLatticeIsingSystem Metropolis update (2)") {
//	const int dim = 2;
//	std::vector<int> L(dim);
//	L[0] = 4;
//	L[1] = 3;
//	std::vector<double> beta = { 0.5 };
//	HoneycombLatticeDataBundle data_bundle(beta);
//	HoneycombLatticeIsingSystem test_system(L, &data_bundle);
//
//	const int seed = 745099137;
//	RandomNumberGenerator mtwist(seed);
//	test_system.setup_RNGen(&mtwist);
//	/* cf. https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution */
//	// 0.455093 0.137216 0.967235 0.771364 0.244447 0.565477 0.8711 0.441108 0.999349 0.340951 0.907204 0.945926
//
//	test_system.set_state_by_code(1339);
//	/*
//	 + - + -
//	 + + - -
//	 + + - +
//	 */
//
//	test_system.update_spins_Metropolis_sequential();
//	/*
//	 - - - -
//	 + - + +
//	 + - + -
//	 */
//
//	REQUIRE(test_system._sz(0) == 1);
//	REQUIRE(test_system._sz(1) == -1);
//	REQUIRE(test_system._sz(2) == 1);
//	REQUIRE(test_system._sz(3) == -1);
//	REQUIRE(test_system._sz(4) == 1);
//	REQUIRE(test_system._sz(5) == -1);
//	REQUIRE(test_system._sz(6) == 1);
//	REQUIRE(test_system._sz(7) == 1);
//	REQUIRE(test_system._sz(8) == -1);
//	REQUIRE(test_system._sz(9) == -1);
//	REQUIRE(test_system._sz(10) == -1);
//	REQUIRE(test_system._sz(11) == -1);
//}

