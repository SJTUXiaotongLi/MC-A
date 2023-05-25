cmake_minimum_required(VERSION 3.14)
project(ProjectIsing CXX)
set(CMAKE_CXX_STANDARD 14)

file( GLOB UTILS Utils/*.cpp )
file( GLOB XML Utils/xml-parser/pugixml/pugixml.cpp )

add_subdirectory(UnitTests/Catch2)
#Include(FetchContent)
#FetchContent_Declare(
#  Catch2
#  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
#  GIT_TAG        v3.3.2 # or a later release
#)
#FetchContent_MakeAvailable(Catch2)

##### exact #####
### Target (exact; Honeycomb lattice)
file( GLOB SRCS1_SQ models/Honeycomb/Ising_Honeycomb_exact.cpp ${UTILS} ${XML} *.hpp )
add_executable(Ising_Honeycomb_exact ${SRCS1_SQ})
target_compile_options(Ising_Honeycomb_exact PRIVATE -Wall -Wextra -O3 -static -DNDEBUG)

##### prototype of MC #####
### Target (prototype MC; Honeycomb lattice)
file( GLOB SRCS2_SQ models/Honeycomb/Ising_Honeycomb_Prototype.cpp ${UTILS} ${XML} *.hpp )
add_executable(Ising_Honeycomb_MC_Prototype ${SRCS2_SQ})
target_compile_options(Ising_Honeycomb_MC_Prototype PRIVATE -Wall -Wextra -O3 -static -DNDEBUG)

##### normal MC #####
### Target (normal MC; Honeycomb lattice)
file( GLOB SRCS2_SQ models/Honeycomb/Ising_Honeycomb_MC.cpp ${UTILS} ${XML} *.hpp )
add_executable(Ising_Honeycomb_MC ${SRCS2_SQ})
target_compile_options(Ising_Honeycomb_MC PRIVATE -Wall -Wextra -O3 -static -DNDEBUG)

###### Unit tests
enable_testing()
file( GLOB SRCS_TEST_MISC1 UnitTests/*.hpp UnitTests/*.cpp)
file( GLOB SRCS_TEST_MISC2 Utils/xml-parser/UnitTests/*_test.cpp)
 


### Unit test (misc 2)
add_executable(test_misc2 ${UTILS} ${XML} ${SRCS_TEST_MISC2})
target_compile_options(test_misc2 PRIVATE -Wall -Wextra -O2)
target_link_libraries(test_misc2 PRIVATE Catch2::Catch2WithMain)
add_test(test_misc2 test_misc2)


### Unit test (misc 1)
add_executable(test_misc1 ${UTILS} ${XML} ${SRCS_TEST_MISC1})
target_compile_options(test_misc1 PRIVATE -Wall -Wextra -O2)
target_link_libraries(test_misc1 PRIVATE Catch2::Catch2WithMain)
add_test(test_misc1 test_misc1)

### Unit test (normal MC; Honeycomb lattice)
add_executable(test_Ising_Honeycomb UnitTests/models/Ising_Honeycomb_test.cpp ${UTILS} ${XML})
target_compile_options(test_Ising_Honeycomb PRIVATE -Wall -Wextra -O2)
target_link_libraries(test_Ising_Honeycomb PRIVATE Catch2::Catch2WithMain)
add_test(test_Ising_Honeycomb test_Ising_Honeycomb)