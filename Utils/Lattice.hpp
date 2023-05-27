#ifndef Lattice_hpp
#define Lattice_hpp

#include <iostream>
#include <vector>
#include <cassert>

class Site {
public:
	static void set_lattice_info(int N_SL_spec, const std::vector<int>& L_spec) {
		/* set dim, L, N_SL */
		dim = static_cast<int>(L_spec.size());
		L = L_spec;
		N_SL = N_SL_spec;
	};
	
	static int eval_site_index(int sl_index, const std::vector<int>& n) {
		/*
		 s_idx = (n[d-1], n[d-2], ..., n[0], sl_index)
		 */
		int s_idx = n[dim - 1];
		for ( int a = dim - 2; a >= 0; a-- ) {
			assert(n[a] < L[a] && n[a] >= 0);
			s_idx *= L[a];
			s_idx += n[a];
		}
		assert(sl_index < N_SL && sl_index >= 0);
		s_idx *= N_SL;
		s_idx += sl_index;
		
		return s_idx;
	};
	
	static void eval_position_at(int s_idx, int& sl_index, std::vector<int>& n) {
		sl_index = s_idx % N_SL;
		s_idx /= N_SL;
		for (int a = 0; a < dim; a++) {
			n[a] = s_idx % L[a];
			s_idx /= L[a];
		}
	};

private:
	static int dim; /* dimension */
	static int N_SL; /* number of sublattices */
	static std::vector<int> L; /* system size */
	
	int site_index;
	int sublattice_index;
	std::vector<int> coordinate;
	int z; /* coordination number */
	std::vector<int> NN; /* nearest neighbor sites */

public:
	Site() : z(0) {};

	virtual ~Site() {};
	
	int get_site_index() const { return site_index; };
	
	void set_site_index() {
		site_index = eval_site_index(sublattice_index, coordinate);
	};
	
	void get_position(int& sl_index, std::vector<int>& r) const {
		sl_index = sublattice_index;
		r = coordinate;
	};
	
	void set_position(int sl_index, const std::vector<int>& r) {
		sublattice_index = sl_index;
		coordinate = r;
		set_site_index();
	};

	int get_z() const { return z; };
	
	void set_z(int z_spec) {
		z = z_spec;
		NN.resize(z);
	};
	
	int get_NN(int bond_index) const {
		assert(bond_index < z && bond_index >= 0);
		return NN[bond_index];
	};
	
	void set_NN(int bond_index, int i_site) {
		assert(bond_index < z && bond_index >= 0);
		NN[bond_index] = i_site;
	};
	
	void print_info() const {
		std::cout
		<< "Site #" << site_index << ": "
		<< "sublattice " << sublattice_index << "; coordinate ("
		<< coordinate[0];
		for (int a = 1; a < dim; a++) {
			std::cout << ", " << coordinate[a];
		}
		std::cout << "); NN sites (if any) ...";
		for (int mu = 0; mu < z; mu++) {
			std::cout << NN[mu] << ' ';
		}
		std::cout << std::endl;
	};
};

class Lattice {
	/* base class for specific lattices */
protected:
	const int dim;
	const int N_SL; /* number of sublattices */
	const std::vector<int> L; /* system size */
	int N_sites;
	std::vector<int> z; /* coordination number for each sublattice */
	std::vector<Site> LatticeSite;
	
public:
	Lattice( int N_SL_spec,
			  const int z_spec[],
			  const std::vector<int>& L_spec)
	: dim(static_cast<int>(L_spec.size())), N_SL(N_SL_spec), L(L_spec) {
		Site::set_lattice_info(N_SL, L);
		N_sites = N_SL;
		for (int a = 0; a < dim; a++) N_sites *= L[a];
		
		z.resize(N_SL);
		for (int mu = 0; mu < N_SL; mu++) z[mu] = z_spec[mu];
		
		LatticeSite.resize(N_sites);
		for (int i = 0; i < N_sites; i++) {
			std::vector<int> coordinate(dim);
			int sl_index;
			Site::eval_position_at(i, sl_index, coordinate);
			LatticeSite[i].set_position(sl_index, coordinate);
			LatticeSite[i].set_z(z[sl_index]);
		}
	};
	
	virtual ~Lattice() {};

	int _NN_of_Site(int site_index, int bond_index) const {
		assert(bond_index < LatticeSite[site_index].get_z() && bond_index >= 0);
		return LatticeSite[site_index].get_NN(bond_index);
	};
	
	std::vector<int> _L() const { return L; }

	int _N_SL() const { return N_SL; };
	
	int _N_sites() const { return N_sites; };
	
	int _z(int site_index) const { return LatticeSite[site_index].get_z(); };

	void print_info_Site(int site_index) const {
		LatticeSite[site_index].print_info();
	};

protected:
	void setup_NN_network() {
		for (int i = 0; i < N_sites; i++) {
			std::vector<int> coordinate(dim);
			int sl_index;
			LatticeSite[i].get_position(sl_index, coordinate);
			for (int bond_index = 0; bond_index < z[sl_index]; bond_index++) {
				int NN_site_index = eval_NN_site_index(sl_index, coordinate, bond_index);
				LatticeSite[i].set_NN(bond_index, NN_site_index);
			}
		}
	};
	
	virtual int eval_NN_site_index(int sl_index, const std::vector<int>& r, int bond_index) = 0;
	
	void move_forward(std::vector<int>& r, int a) {
		assert(a < dim && a >= 0);
		r[a] += 1;
		r[a] = r[a] % L[a];
	};
	
	void move_backward(std::vector<int>& r, int a) {
		assert(a < dim && a >= 0);
		r[a] += L[a] - 1;
		r[a] = r[a] % L[a];
	};
};

#endif /* Lattice_hpp */
