#include "Lattice.hpp"
#include "../models/Square/Ising_Square.hpp"
#include "../models/Square/Ising_Honey.hpp"

int Site::dim;
int Site::N_SL;
std::vector<int> Site::L;

const int SquareLattice::N_SL_SqLatt;
const int SquareLattice::z_common;
const int SquareLattice::z_common_half;
constexpr int SquareLattice::z_SqLatt[N_SL_SqLatt];

const int HoneyLattice::N_SL_HnLatt;
const int HoneyLattice::z_common;
// const int HoneyLattice::z_common_half;
constexpr int HoneyLattice::z_HnLatt[N_SL_HnLatt];
