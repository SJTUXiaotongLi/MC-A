#ifndef IsingSystem_hpp
#define IsingSystem_hpp

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <cassert>
#include "Lattice.hpp"
#include "Observable.hpp"
#include "RandomNumberGenerator.hpp"

class IsingSpin {
private:
	int sz; /* +/- 1 */
	
public:
	IsingSpin() {	sz = 1; };
	
	~IsingSpin() {};
	
	int _sz() const { return sz; };
	
	void set_sz(int sz_spec) {
		assert(sz_spec == 1 || sz_spec == -1);
		sz = sz_spec;
	};
	
	void flip() {	sz *= -1; };
};

class IsingSystem {
protected:
	static RandomNumberGenerator* ptr_mtwist;

	const double J;
	const int n_spins;
	const long long maxrep_state;
	std::vector<IsingSpin> spin;


public:
	IsingSystem(const int n_spins_spec) : J(-1.0), n_spins(n_spins_spec), maxrep_state(static_cast<long long>(std::pow(2, n_spins)) - 1) {
		spin.resize(n_spins);
	};

	virtual ~IsingSystem() {};
	
	static void setup_RNGen(RandomNumberGenerator& mtwist) {
		ptr_mtwist = &mtwist;
	};

	double _J() const { return J; };
	
	int _n_spins() const { return n_spins;	};
	
	long long _maxrep_state() const { return maxrep_state; };
	
	int _sz(const int site_idx) const { return spin[site_idx]._sz();	}
	
	void set_spin(const int site_idx, int s_spec) {	spin[site_idx].set_sz(s_spec); };
	
	void flip_spin(const int site_idx) { spin[site_idx].flip(); };
	
	double eval_mz() const {
		int mz_tmp = 0;
		for (int site_idx = 0; site_idx < n_spins; site_idx++) {
			mz_tmp += spin[site_idx]._sz();
		}
		return static_cast<double>(mz_tmp);
	};
	
	void set_state_by_code(long long rep_state) {
		assert(rep_state >= 0 && rep_state <= maxrep_state);
		for (int site_idx = 0; site_idx < n_spins; site_idx++) {
			set_spin(site_idx, (rep_state % 2) ? 1 : -1);
			rep_state /= 2;
		}
	};
	
	void randomize_spins() {
		for (int site_idx = 0; site_idx < n_spins; site_idx++) {
			set_spin(site_idx, (ptr_mtwist->gen_rand01() > 0.5) ? 1 : -1);
		}
	};
};

template<class LatticeType, class LatticeDataBundle>
class LatticeIsingSystem : public IsingSystem {
protected:
	LatticeType lattice;
	LatticeDataBundle* ptrMCDB;
	const int beta_list_size;
	const int mcs_thermalization;
	const int mcs_sampling;
	const int mcs_total;
	const int mcs_interval_btwn_bins;
	int beta_index_current;
	int mcs_current;
	int sampling_current_bin;

public:
	LatticeIsingSystem(const std::vector<int>& L_spec, int n_spins_spec, LatticeDataBundle& MCDB_spec, const int mcs_thermalization_spec = 0, const int mcs_interval_btwn_bins_spec = 0)
	: IsingSystem(n_spins_spec), lattice(L_spec), ptrMCDB(&MCDB_spec), beta_list_size(static_cast<int>(MCDB_spec._beta().size())), mcs_thermalization(mcs_thermalization_spec), mcs_sampling(MCDB_spec._mcs_sampling()), mcs_total(mcs_thermalization_spec + MCDB_spec._mcs_sampling()), mcs_interval_btwn_bins(mcs_interval_btwn_bins_spec), beta_index_current(0), mcs_current(0), sampling_current_bin(0) {};

	LatticeIsingSystem(const std::vector<int>& L_spec, int n_spins_spec)
	: IsingSystem(n_spins_spec), lattice(L_spec), ptrMCDB(nullptr), beta_list_size(0), mcs_thermalization(0), mcs_sampling(0), mcs_total(0), mcs_interval_btwn_bins(0), beta_index_current(0), mcs_current(0), sampling_current_bin(0) {};
	
	virtual ~LatticeIsingSystem() { ptrMCDB = nullptr; };
	
	std::vector<int> _L() const { return lattice._L(); };
	
	double _beta(int beta_idx) const {
		assert(ptrMCDB != nullptr);
		assert(beta_idx >= 0 && beta_idx < beta_list_size);
		return ptrMCDB->_beta(beta_idx);
	};
	
	double _T(int beta_idx) const { return 1.0 / _beta(beta_idx); };

	void set_beta_index(int b_idx) { beta_index_current = b_idx; }
	
	void exact_count() {
		assert(ptrMCDB != nullptr);
		double ene, mz, weight;
		long long rep_state = 0;
		double ene_min = 0;
		while (rep_state <= maxrep_state) {
			set_state_by_code(rep_state);
			ene = eval_energy();
			mz = eval_mz();
			
			if (ene < ene_min) { /* weight correction */
				const double dE = ene_min - ene;
				for ( int beta_idx = 0; beta_idx < beta_list_size; beta_idx++ ) {
					const double w_correction = std::exp(-1.0 * _beta(beta_idx) * dE);
					ptrMCDB->weight_correction(beta_idx, w_correction);
				}
				ene_min = ene;
			}
			
			for ( int beta_idx = 0; beta_idx < beta_list_size; beta_idx++ ) {
				weight = std::exp(-1.0 * _beta(beta_idx) * (ene - ene_min));
				ptrMCDB->update_exact_energy(beta_idx, ene, weight);
				ptrMCDB->update_exact_magz(beta_idx, mz, weight);
			}
			
			rep_state++;
		}
		ptrMCDB->normalize_by_Z();
	};
	
	virtual double eval_energy() const = 0;
	
	void reset_mcs() { mcs_current = 0; };
	
	void reset_sampling_current_bin() { sampling_current_bin = 0; };
	
	void run_MC(bool print_output = true) {
		randomize_spins();
		while (! _is_MC_finished_all()) {
			run_MC_given_T();
			MC_finalize(beta_index_current, print_output);
			change_beta_index();
		}
	};

	void run_MC_prototype(int beta_idx) {
		const bool need_exact = (n_spins <= 32);
		if ( need_exact ) exact_count();
		const double ene_exact = ptrMCDB->_get_exact_energy(beta_idx);
		const double magz_exact = ptrMCDB->_get_exact_magz(beta_idx);
		const double mzsq_exact = ptrMCDB->_get_exact_magz_sq(beta_idx);
		
		std::cout << "### MC prototype at T = " << _T(beta_idx) << " : ";
		std::cout << "[1] MCS [2] E [3] M [4] M2 [5] ";
		if ( need_exact ) std::cout << " E(exact) [6] M(exact) [7] M2(exact)";
		std::cout << std::endl;
		
		if ( need_exact ) {
			std::cout << "### MC prototype at T = " << _T(beta_idx) << " : ";
			std::cout << "# E = " << ene_exact << std::endl;
			std::cout << "### MC prototype at T = " << _T(beta_idx) << " : ";
			std::cout << "# Mz = " << magz_exact << std::endl;
			std::cout << "### MC prototype at T = " << _T(beta_idx) << " : ";
			std::cout << "# Mz2 = " << mzsq_exact << std::endl;
		}

		beta_index_current = beta_idx;
		randomize_spins();
		reset_mcs();
		while ( mcs_current < mcs_total) {
			update_spins_Metropolis_random();
			const double ene = eval_energy();
			const double mz = eval_mz();

			std::cout << mcs_current << ' ';
			std::cout << ene << ' ' << mz << ' ' << mz * mz << ' ';
			if ( need_exact ) {
				std::cout << ene_exact << ' ' << magz_exact << ' ' << mzsq_exact << ' ';
			}
			std::cout << std::endl;
			mcs_current++;
		}
	};
	
	void update_spins_Metropolis_sequential() {
		for (int site_idx = 0; site_idx < lattice._N_sites(); site_idx++) {
			spinflip_Metropolis(site_idx);
		}
	};
	
	void update_spins_Metropolis_random() {
		for (int counter = 0; counter < lattice._N_sites(); counter++) {
			spinflip_Metropolis(ptr_mtwist->gen_rand_site());
		}
	};
	
	void update_estimators() {
		assert(ptrMCDB != nullptr);
		const double ene = eval_energy();
		const double mz = eval_mz();
		ptrMCDB->update_sampling_energy(beta_index_current, ene);
		ptrMCDB->update_sampling_magz(beta_index_current, mz);
		sampling_current_bin++;
	};
	
protected:
	void MC_finalize(int b_idx, bool print_output) {
		ptrMCDB->MC_normalize(b_idx);
		if (print_output) ptrMCDB->output_data_MC(b_idx);
	};

private:
	void change_beta_index() {
		if (++beta_index_current == beta_list_size) beta_index_current = -1;
	};
	
	bool _is_MC_finished_all() const { return (beta_index_current == -1); };

	void run_MC_given_T() {
		reset_mcs();
		thermalize_given_T();
		sampling_given_T();
	};
	
	void thermalize_given_T() {
		while ( mcs_current < mcs_thermalization ) {
			update_spins_Metropolis_random();
			mcs_current++;
		}
	};
	
	void sampling_given_T() {
		assert(mcs_current >= mcs_thermalization);
		while ( mcs_current < mcs_total) {
			update_spins_Metropolis_random();
			update_estimators();
			check_bin_status();
			mcs_current++;
		}
	};
	
	void check_bin_status() {
		if ( ptrMCDB->check_if_bin_filled(sampling_current_bin) ) {
			reset_sampling_current_bin();
			ptrMCDB->switch_bin(beta_index_current);
			for (int n = 0; n < mcs_interval_btwn_bins; n++) {
				update_spins_Metropolis_random();
			}
		}
	};
	
	void spinflip_Metropolis(int site_idx) {
		if (ptr_mtwist->gen_rand01() < eval_prob_spinflip_Metropolis(site_idx)) {
			spin[site_idx].flip();
		}
	};
	
	double eval_prob_spinflip_Metropolis(int site_idx) const {
		assert(site_idx >= 0 && site_idx < lattice._N_sites());
		return std::exp(-2.0 * _beta(beta_index_current) * eval_weiss_at(site_idx) * spin[site_idx]._sz());
	};
	
	double eval_weiss_at(int site_idx) const {
		assert(site_idx >= 0 && site_idx < lattice._N_sites());
		double hMF = 0.0;
		const int z_all =  lattice._z_common();
		for (int bond_index = 0; bond_index < z_all; bond_index++) {
			const int j = lattice._NN_of_Site(site_idx, bond_index);
			hMF += spin[j]._sz();
		}
		hMF *= -J;
		return hMF;
	};
};

#endif /* IsingSystem_hpp */
