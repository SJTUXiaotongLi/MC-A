#ifdef Xcode
#include "Catch2_amalgamated/catch_amalgamated.hpp"
#else
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#endif

#include "../ParameterBundle.hpp"

/* cf. sample/sample.xml */

TEST_CASE("XML parameter paser test 2") {
#ifdef Xcode
	SomeLatticeParameterBundle Parameters("sample/sample.xml", true);
#else
	SomeLatticeParameterBundle Parameters("../Utils/xml-parser/sample/sample.xml", true);
#endif
	REQUIRE(Parameters._L0() == 12);
	REQUIRE(Parameters._L1() == 34);
	REQUIRE(Parameters._N_TemperaturePoints() == 56);
	REQUIRE_THAT(Parameters._Tmax(), Catch::Matchers::WithinRel(6.5, 1e-10));
	REQUIRE_THAT(Parameters._Tmin(), Catch::Matchers::WithinRel(0.1, 1e-10));
}

