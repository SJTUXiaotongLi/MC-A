#ifdef Xcode
#include "Catch2_amalgamated/catch_amalgamated.hpp"
#else
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#endif

#include "../Utils/Lattice.hpp"
#include "../Utils/IsingSystem.hpp"
#include "../Utils/RandomNumberGenerator.hpp"

TEST_CASE("IsingSpin") {
	IsingSpin test_spin;
	REQUIRE(test_spin._sz() == 1);
	
	test_spin.set_sz(-1);
	REQUIRE(test_spin._sz() == -1);
	
	test_spin.flip();
	REQUIRE(test_spin._sz() == 1);
}

TEST_CASE("IsingSystem") {
	const int n_spins = 4;
	IsingSystem test_system(n_spins);
	REQUIRE(test_system._J() == -1);
	
	REQUIRE(test_system._sz(0) == 1);
	REQUIRE(test_system._sz(1) == 1);
	REQUIRE(test_system._sz(2) == 1);
	REQUIRE(test_system._sz(3) == 1);
	REQUIRE(test_system.eval_mz() == 4.0);
	
	test_system.flip_spin(0);
	REQUIRE(test_system._sz(0) == -1);
	REQUIRE(test_system.eval_mz() == 2.0);
	
	test_system.set_spin(2, -1);
	REQUIRE(test_system._sz(2) == -1);
	REQUIRE(test_system.eval_mz() == 0.0);
	
	test_system.set_spin(3, -1);
	REQUIRE(test_system._sz(3) == -1);
	REQUIRE(test_system.eval_mz() == -2.0);
}

TEST_CASE("IsingSystem Randomization") {
	const int seed = 745099137;
	RandomNumberGenerator mtwist(seed);
	IsingSystem::setup_RNGen(mtwist);

	const int n_spins = 10;
	IsingSystem test_system(n_spins);
	
	/* cf. https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution */
	// 0.455093 0.137216 0.967235 0.771364 0.244447 0.565477 0.8711 0.441108 0.999349 0.340951
	
	test_system.randomize_spins();
	REQUIRE(test_system._sz(0) == -1);
	REQUIRE(test_system._sz(1) == -1);
	REQUIRE(test_system._sz(2) == 1);
	REQUIRE(test_system._sz(3) == 1);
	REQUIRE(test_system._sz(4) == -1);
	REQUIRE(test_system._sz(5) == 1);
	REQUIRE(test_system._sz(6) == 1);
	REQUIRE(test_system._sz(7) == -1);
	REQUIRE(test_system._sz(8) == 1);
	REQUIRE(test_system._sz(9) == -1);
}
