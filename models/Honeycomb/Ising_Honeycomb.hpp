#ifndef Ising_Honeycomb_HPP
#define Ising_Honeycomb_HPP

#include "../../Utils/IsingSystem.hpp"
#include "../../Utils/xml-parser/ParameterBundle.hpp"

class HoneycombLatticeParameterBundle_exact : public ParameterBundle {
public:
	HoneycombLatticeParameterBundle_exact(std::string input_file, bool print_option_spec = false) : ParameterBundle(print_option_spec) {
		init(input_file);
	};
	
	void read_xml(std::string input_file) {
		pugi::xml_document doc;
		if ( ! doc.load_file(input_file.c_str()) ) {
			std::cerr << input_file << " is not found" << std::endl;
			exit(1);
		}
		read_xml_category(doc, "ID");
		read_xml_category(doc, "System");
		read_xml_category(doc, "Model");
	};

	void add_parameters() {
		ParameterBundle::add_parameters();
		parameters.push_back(Parameter("L0", "int", print_option, ">", 1));
		parameters.push_back(Parameter("L1", "int", print_option, ">", 1));
		parameters.push_back(Parameter("N_TemperaturePoints", "int", print_option, ">", 1));
		parameters.push_back(Parameter("Tmax", "double", print_option, ">", 0));
		parameters.push_back(Parameter("Tmin", "double", print_option, ">", 0));
	};
	
	int _L0() { return _value_int("L0"); };
	int _L1() { return _value_int("L1"); };
	std::vector<int> _L() {
		std::vector<int> L = { _L0(), _L1() };
		return L;
	};
	int _N_TemperaturePoints() { return _value_int("N_TemperaturePoints"); };
	double _Tmin() { return _value_dbl("Tmin"); };
	double _Tmax() { return _value_dbl("Tmax"); };
};

class HoneycombLatticeParameterBundle_MC : public ParameterBundle {
public:
	HoneycombLatticeParameterBundle_MC(std::string input_file, bool print_option_spec = false) : ParameterBundle(print_option_spec) {
		init(input_file);
	};
	
	void read_xml(std::string input_file) {
		pugi::xml_document doc;
		if ( ! doc.load_file(input_file.c_str()) ) {
			std::cerr << input_file << " is not found" << std::endl;
			exit(1);
		}
		read_xml_category(doc, "ID");
		read_xml_category(doc, "System");
		read_xml_category(doc, "Model");
		read_xml_category(doc, "MC");
	};

	void add_parameters() {
		ParameterBundle::add_parameters();
		parameters.push_back(Parameter("L0", "int", print_option, ">", 1));
		parameters.push_back(Parameter("L1", "int", print_option, ">", 1));
		parameters.push_back(Parameter("N_TemperaturePoints", "int", print_option, ">", 1));
		parameters.push_back(Parameter("Tmax", "double", print_option, ">", 0));
		parameters.push_back(Parameter("Tmin", "double", print_option, ">", 0));
		parameters.push_back(Parameter("n_bins", "int", print_option, ">", 1));
		parameters.push_back(Parameter("mcs_thermalization", "int", print_option, ">", 0));
		parameters.push_back(Parameter("n_samples_per_bin", "int", print_option, ">", 0));
		parameters.push_back(Parameter("mcs_interval_btwn_bins", "int", print_option, ">", 0));
	};
	
	int _L0() { return _value_int("L0"); };
	int _L1() { return _value_int("L1"); };
	std::vector<int> _L() {
		std::vector<int> L = { _L0(), _L1() };
		return L;
	};
	int _N_TemperaturePoints() { return _value_int("N_TemperaturePoints"); };
	double _Tmin() { return _value_dbl("Tmin"); };
	double _Tmax() { return _value_dbl("Tmax"); };
	int _n_bins() { return _value_int("n_bins"); };
	int _mcs_thermalization() { return _value_int("mcs_thermalization"); };
	int _n_samples_per_bin() { return _value_int("n_samples_per_bin"); };
	int _mcs_interval_btwn_bins() { return _value_int("mcs_interval_btwn_bins"); };
};

class HoneycombLattice : public Lattice {
private:
	static const int N_SL_HcLatt = 2;
	static const int z_common = 3;
	static constexpr int z_HcLatt[N_SL_HcLatt] = { z_common, z_common };

public:
	HoneycombLattice(const std::vector<int>& L_spec)
	: Lattice(N_SL_HcLatt, z_HcLatt, L_spec) {
		setup_NN_network();
	};
	
	~HoneycombLattice() {};

	int _z_common() const { return z_common; };

	int _L0() const { return L[0]; };
	int _L1() const { return L[1]; };
	std::vector<int> _L() {
		std::vector<int> L = { _L0(), _L1() };
		return L;
	};

private:
	int eval_NN_site_index(int sl_index, const std::vector<int>& r, int bond_index) {
		std::vector<int> NN_r;
		int NN_sl_index;
		/*sublattice index == 0*/
		if (sl_index == 0)
		{
			switch (bond_index)
			{
			case 0:
				NN_r = r;
				NN_sl_index = 1;
				break;
			case 1:
				NN_r = r;
				move_backward(NN_r, 0);
				NN_sl_index = 1;
				break;
			case 2:
				NN_r = r;
				move_backward(NN_r, 1);
				NN_sl_index = 1;
				break;
			default:
				std::cerr << "HoneycombLattice::eval_NN_site_index> error" << std::endl;
				exit(0);
			}
		}
		/*sublattice index == 1*/
		if (sl_index == 1)
		{
			switch (bond_index)
			{
			case 0:
				NN_r = r;
				move_forward(NN_r, 0);
				NN_sl_index = 0;
				break;
			case 1:
				NN_r = r;
				move_forward(NN_r, 1);
				NN_sl_index = 0;
				break;
			case 2:
				NN_r = r;
				NN_sl_index = 0;
				break;
			default:
				std::cerr << "HoneycombLattice::eval_NN_site_index> error" << std::endl;
				exit(0);
			}
		}

		return Site::eval_site_index(NN_sl_index, NN_r);
	}
};

class HoneycombLatticeDataBundle : public DataBundle {
public:
	HoneycombLatticeDataBundle(int n_spins_spec, const std::vector<double>& beta_spec, int n_bins_spec = 0, int n_samples_per_bin_spec = 0) : DataBundle(n_spins_spec, beta_spec, n_bins_spec, n_samples_per_bin_spec) {};
	
	~HoneycombLatticeDataBundle() {};
	
	void MC_normalize(int beta_idx) {
		assert(check_if_complete_sampling_beta_index(beta_idx));
		DataBundle::MC_normalize(beta_idx);
	};
	
	void switch_bin(int beta_idx) {
		DataBundle::switch_bin(beta_idx);
	};
	
	void output_legends_MC(const std::string system_size) {
		DataBundle::output_legends_MC(system_size);
		std::cout << std::endl;
	};
	
	void output_data_MC(int beta_idx) {
		DataBundle::output_data_MC(beta_idx);
		std::cout << std::endl;
	};
	
	void output_legends_exact(const std::string system_size) {
		DataBundle::output_legends_exact(system_size);
		std::cout << std::endl;
	};
	
	void output_data_exact(int beta_idx) {
		DataBundle::output_data_exact(beta_idx);
		std::cout << std::endl;
	};
	
	void output_data_exact_all() {
		for ( unsigned int beta_idx = 0; beta_idx < beta.size(); beta_idx++ ) {
			output_data_exact(beta_idx);
		}
	};
};

class HoneycombLatticeIsingSystem : public LatticeIsingSystem<HoneycombLattice, HoneycombLatticeDataBundle> {
public:
	HoneycombLatticeIsingSystem(const std::vector<int>& L_spec, HoneycombLatticeDataBundle& MCDB_spec, const int mcs_thermalization_spec = 0, const int mcs_interval_btwn_bins_spec = 0)
	: LatticeIsingSystem<HoneycombLattice, HoneycombLatticeDataBundle>(L_spec, eval_n_spins(L_spec), MCDB_spec, mcs_thermalization_spec, mcs_interval_btwn_bins_spec) {};

	HoneycombLatticeIsingSystem(const std::vector<int>& L_spec) : LatticeIsingSystem<HoneycombLattice, HoneycombLatticeDataBundle>(L_spec, eval_n_spins(L_spec)){};
	
	~HoneycombLatticeIsingSystem() {};
	
	double eval_energy() const {
		double energy_tmp = 0;
		for (int site_idx = 0; site_idx < lattice._N_sites(); site_idx += 2) {
			for (int bond_index = 0; bond_index < 3; bond_index++) {
				const int j = lattice._NN_of_Site(site_idx, bond_index);
				energy_tmp += J * spin[site_idx]._sz() * spin[j]._sz();
			}
		}
		return energy_tmp;
	};
	
	static int eval_n_spins(const std::vector<int>& L_spec) {
		return 2 * L_spec.at(0) * L_spec.at(1);
	};
	int _L0() const { return lattice._L0(); };
	int _L1() const { return lattice._L1(); };
	std::vector<int> _L() { return lattice._L(); };
	std::string _system_size() const { return "# System size : " + std::to_string(_L0()) + " x " + std::to_string(_L1()); };
};

#endif /* Ising_Honeycomb_HPP */