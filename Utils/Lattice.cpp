#include "Lattice.hpp"
#include "../models/Honeycomb/Ising_Honeycomb.hpp"

int Site::dim;
int Site::N_SL;
std::vector<int> Site::L;

const int HoneycombLattice::N_SL_SqLatt;
const int HoneycombLattice::z_common;
constexpr int HoneycombLattice::z_SqLatt[N_SL_SqLatt];
