/*
 * ContactShuffler.h
 *
 *  Created on: Nov 30, 2016
 *      Author: nettam
 */

#ifndef CONTACTSHUFFLER_H_
#define CONTACTSHUFFLER_H_
#include <vector>
#include <map>
#include <unordered_map>
#include <iostream>
using namespace std;

class GenomeGridLog;

// A contact as stored in its grid cell: its index in the contact list, its
// coordinates and its distance bin. Keeping the coordinates in the cell saves
// the pointer chase to the contact list on every partner draw.
struct GridContact {
	int idx;
	int x;
	int y;
	int dist_bin;
};

class ContactShuffler final {
public:
	ContactShuffler(int dist_log_scale,	int dist_resolution,
			int grid_x_resolution, //int grid_dist_resolution,
			int grid_switch_bin_dist, //int grid_switch_x_dist,
			float correction_factor, int decay_smooth, int regularization, int min_dist, int max_dist);
	virtual ~ContactShuffler();
	// contact i is (x[i*stride], y[i*stride]), i < n
	long load_contacts(const int* x, const int* y, int stride, long n, bool symetric);
	int save_contacts(const char* fn, bool symetric, bool with_header);
	int init_obs_decay_from_contacts();
	int	init_exp_decay_from_obs();

	int init_proposal_const();
	int init_proposal_from_area();
	int init_proposal_from_contacts(long proposal_shuffle);
	void correct_proposal_dist();
	int shuffle_contacts(int shuffle_factor, float transition_correction_factor=0.0001,
			float transition_cooling_update=0.01, int debug=0);
	void save_transitions(ostream& out, int id);
	void debug(ostream& out, int id);
	void reset_grid(int grid_x_binsize);
protected:
	void	init_contact_dist_bins();
	int		simple_sample();
	int		get_dist_bin(int x, int y);
	float 	get_bin_dist(int bin);
	int		get_grid_bin(int x);
	void 	regularize_decay(vector<float>& decay);
	void 	select_switch_partners(int contact, int& cell1, int& grid_index1, int& cell2, int& grid_index2);
	void 	grid_move(int cell_i, int grid_index_i, int cell_j, int grid_index_j, int dist_ij_bin, int dist_ji_bin);
	void	build_grid();
	void	contacts_from_grid();
	void 	print_proposal(ostream& out);

protected:
	long					m_contact_count;
	int						m_max_contact_dist;
	int 					m_dist_log_scale;
	int						m_dist_resolution;
	float					m_log_log_scale;
	// contact coordinates by index; filled only while loading, rebuilding
	// the grid and saving (the grid cells hold the current coordinates)
	vector<int>				m_x;
	vector<int>				m_y;
	vector <int> 			m_contacts_dist_bins;
	vector <int>			m_transitions;
	int						m_grid_x_binsize;
	int						m_grid_switch_bin_dist;
	int						m_grid_size;
	// m_grid_size x m_grid_size cells, cell (bin1, bin2) at bin1 * m_grid_size + bin2
	vector< vector<GridContact> > m_contact_grid;
	vector<int>				m_contact_cell;		// grid cell of each contact
	vector<int>				m_pool_cumsum;		// scratch for select_switch_partners
	vector<int>				m_pool_cell;

	vector <float>			m_proposal_freq;
	vector<unsigned long>	m_decay_obs;
	vector <float>			m_decay_exp;
	vector <char>			m_decay_exp_nonzero;	// exp(m_decay_exp[bin]) != 0, fixed after init
	vector <unsigned long>	m_log_obs_count;	// m_decay_obs[bin] when m_log_obs[bin] was computed
	vector <double>			m_log_obs;
	int						m_min_x;
	int						m_max_x;
	float					m_correction_factor;
	int						m_regularization;
	int						m_proposal_shuffle;
	int						m_min_dist;
	int						m_decay_smooth;
    float					m_reg;

};


#endif /* CONTACTSHUFFLER_H_ */
