#ifdef Xcode
#include "Catch2_amalgamated/catch_amalgamated.hpp"
#else
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#endif

#include "../Utils/RandomNumberGenerator.hpp"

TEST_CASE("Random Number (64-bit Mersenne Twister in C++ library)") {
	const int seed = 745099137;
	RandomNumberGenerator mtwist(seed);

	/* cf. https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution */
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.455093, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.137216, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.967235, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.771364, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.244447, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.565477, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.871100, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.441108, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.999349, 1e-1));
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(0.340951, 1e-1));
};

TEST_CASE("Save/Load Random Number (64-bit Mersenne Twister in C++ library)") {
	const int seed = 745099137;
	RandomNumberGenerator mtwist(seed);

	/* cf. https://www.reddit.com/r/cpp/comments/n53g4v/standard_library_random_slow_to_loadsave/ */
	for (int n = 0; n < 99; n++) mtwist.gen_rand01();
	std::string save_data = mtwist.save();
	const double rn_sample = mtwist.gen_rand01();

	for (int n = 0; n < 999; n++) mtwist.gen_rand01();
	mtwist.load(save_data);
	REQUIRE_THAT(mtwist.gen_rand01(), Catch::Matchers::WithinRel(rn_sample, 1e-10));
};
