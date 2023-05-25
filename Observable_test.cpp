#ifdef Xcode
#include "Catch2_amalgamated/catch_amalgamated.hpp"
#else
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#endif

#include "../Utils/Observable.hpp"

TEST_CASE("Observable Usage for a direct summation (demo)") {
	Observable test_obs;
	REQUIRE(test_obs._q() == 0);
	REQUIRE(test_obs._q_abs() == 0);
	REQUIRE(test_obs._q_sq() == 0);
	REQUIRE(test_obs._q_quar() == 0);
	
	test_obs.update_direct_by(10, 1.0);
	REQUIRE(test_obs._q() == 10);
	REQUIRE(test_obs._q_abs() == 10);
	REQUIRE(test_obs._q_sq() == 100);
	REQUIRE(test_obs._q_quar() == 10000);
	
	test_obs.update_direct_by(-10, 1.0);
	REQUIRE(test_obs._q() == 0);
	REQUIRE(test_obs._q_abs() == 20);
	REQUIRE(test_obs._q_sq() == 200);
	REQUIRE(test_obs._q_quar() == 20000);
	
	test_obs.normalize_by_Z();
	REQUIRE(test_obs._q() == 0);
	REQUIRE(test_obs._q_abs() == 10);
	REQUIRE(test_obs._q_sq() == 100);
	REQUIRE(test_obs._q_quar() == 10000);
}

TEST_CASE("MonteCarloObservable Updates & Normalization (demo)") {
	const int n_bin = 2;
	const int n_samples_per_bin = 2;
	MonteCarloObservable test_obs(n_bin, n_samples_per_bin);
	
	test_obs.update_sampling_by(-0.5);
	if (test_obs.check_if_bin_filled()) test_obs.switch_bin();
//	REQUIRE_THAT(test_obs._q_bin(0), Catch::Matchers::WithinRel(-0.5, 1e-10));
//	REQUIRE_THAT(test_obs._q_abs_bin(0), Catch::Matchers::WithinRel(0.5, 1e-10));
//	REQUIRE_THAT(test_obs._q_sq_bin(0), Catch::Matchers::WithinRel(0.25, 1e-10));
//	REQUIRE_THAT(test_obs._q_quar_bin(0), Catch::Matchers::WithinRel(0.0625, 1e-10));
	REQUIRE(test_obs.check_if_complete_sampling() == false);
	
	test_obs.update_sampling_by(0.1);
	if (test_obs.check_if_bin_filled()) test_obs.switch_bin();
//	REQUIRE_THAT(test_obs._q_bin(0), Catch::Matchers::WithinRel(-0.4, 1e-10));
//	REQUIRE_THAT(test_obs._q_abs_bin(0), Catch::Matchers::WithinRel(0.6, 1e-10));
//	REQUIRE_THAT(test_obs._q_sq_bin(0), Catch::Matchers::WithinRel(0.26, 1e-10));
//	REQUIRE_THAT(test_obs._q_quar_bin(0), Catch::Matchers::WithinRel(0.0626, 1e-10));
	REQUIRE(test_obs.check_if_complete_sampling() == false);
	
	test_obs.update_sampling_by(-0.1);
	if (test_obs.check_if_bin_filled()) test_obs.switch_bin();
//	REQUIRE_THAT(test_obs._q_bin(1), Catch::Matchers::WithinRel(-0.1, 1e-10));
//	REQUIRE_THAT(test_obs._q_abs_bin(1), Catch::Matchers::WithinRel(0.1, 1e-10));
//	REQUIRE_THAT(test_obs._q_sq_bin(1), Catch::Matchers::WithinRel(0.01, 1e-10));
//	REQUIRE_THAT(test_obs._q_quar_bin(1), Catch::Matchers::WithinRel(0.0001, 1e-10));
	REQUIRE(test_obs.check_if_complete_sampling() == false);
	
	test_obs.update_sampling_by(-2.0);
	if (test_obs.check_if_bin_filled()) test_obs.switch_bin();
//	REQUIRE_THAT(test_obs._q_bin(1), Catch::Matchers::WithinRel(-2.1, 1e-10));
//	REQUIRE_THAT(test_obs._q_abs_bin(1), Catch::Matchers::WithinRel(2.1, 1e-10));
//	REQUIRE_THAT(test_obs._q_sq_bin(1), Catch::Matchers::WithinRel(4.01, 1e-10));
//	REQUIRE_THAT(test_obs._q_quar_bin(1), Catch::Matchers::WithinRel(16.0001, 1e-10));
	REQUIRE(test_obs.check_if_complete_sampling() == true);
	
	test_obs.MC_normalize();
	REQUIRE_THAT(test_obs._q(), Catch::Matchers::WithinRel(-0.625, 1e-10));
	REQUIRE_THAT(test_obs._q_abs(), Catch::Matchers::WithinRel(0.675, 1e-10));
	REQUIRE_THAT(test_obs._q_sq(), Catch::Matchers::WithinRel(1.0675, 1e-10));
	REQUIRE_THAT(test_obs._q_quar(), Catch::Matchers::WithinRel(4.015675, 1e-10));
	REQUIRE_THAT(test_obs._q_fluc(), Catch::Matchers::WithinRel(0.49625, 1e-10));
	
	REQUIRE_THAT(test_obs._stderr_q(), Catch::Matchers::WithinRel(0.425, 1e-10));
	REQUIRE_THAT(test_obs._stderr_q_abs(), Catch::Matchers::WithinRel(0.375, 1e-10));
	REQUIRE_THAT(test_obs._stderr_q_sq(), Catch::Matchers::WithinRel(0.9375, 1e-10));
	REQUIRE_THAT(test_obs._stderr_q_quar(), Catch::Matchers::WithinRel(3.984375, 1e-10));
	REQUIRE_THAT(test_obs._stderr_q_fluc(), Catch::Matchers::WithinRel(0.40625, 1e-10));
};
