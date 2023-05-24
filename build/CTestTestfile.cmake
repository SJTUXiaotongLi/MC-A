# CMake generated Testfile for 
# Source directory: /home/admini/Downloads/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU CP2022-23
# Build directory: /home/admini/Downloads/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU CP2022-23/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_misc1 "test_misc1")
set_tests_properties(test_misc1 PROPERTIES  _BACKTRACE_TRIPLES "/home/admini/Downloads/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU CP2022-23/CMakeLists.txt;49;add_test;/home/admini/Downloads/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU CP2022-23/CMakeLists.txt;0;")
add_test(test_misc2 "test_misc2")
set_tests_properties(test_misc2 PROPERTIES  _BACKTRACE_TRIPLES "/home/admini/Downloads/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU CP2022-23/CMakeLists.txt;55;add_test;/home/admini/Downloads/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU CP2022-23/CMakeLists.txt;0;")
add_test(test_Ising_Square "test_Ising_Square")
set_tests_properties(test_Ising_Square PROPERTIES  _BACKTRACE_TRIPLES "/home/admini/Downloads/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU CP2022-23/CMakeLists.txt;61;add_test;/home/admini/Downloads/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU CP2022-23/CMakeLists.txt;0;")
add_test(test_Ising_Honey "test_Ising_Honey")
set_tests_properties(test_Ising_Honey PROPERTIES  _BACKTRACE_TRIPLES "/home/admini/Downloads/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU CP2022-23/CMakeLists.txt;67;add_test;/home/admini/Downloads/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU-CP2022-23-Ising-Lecture0424_Square/SJTU CP2022-23/CMakeLists.txt;0;")
subdirs("UnitTests/Catch2")
