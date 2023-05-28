#include "Lattice.hpp"
#include "../models/Square/Ising_Square.hpp"
#include "../models/Honeycomb/Ising_Honeycomb.hpp"
#include "../models/Kagome/Ising_Kagome.hpp"
#include "../models/Cubic/Ising_Cubic.hpp"

int Site::dim;
int Site::N_SL;
std::vector<int> Site::L;

const int SquareLattice::N_SL_SqLatt;
const int SquareLattice::z_common;
const int SquareLattice::z_common_half;
constexpr int SquareLattice::z_SqLatt[N_SL_SqLatt];

const int HoneycombLattice::N_SL_HcLatt;
const int HoneycombLattice::z_common;
constexpr int HoneycombLattice::z_HcLatt[N_SL_HcLatt];

const int KagomeLattice::N_SL_KaLatt;
const int KagomeLattice::z_common;
const int KagomeLattice::z_common_half;
constexpr int KagomeLattice::z_KaLatt[N_SL_KaLatt];

const int CubicLattice::N_SL_CuLatt;
const int CubicLattice::z_common;
const int CubicLattice::z_common_half;
constexpr int CubicLattice::z_CuLatt[N_SL_CuLatt];

