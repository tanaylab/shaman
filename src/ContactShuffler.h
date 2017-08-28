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

class ContactShuffler {
public:
	ContactShuffler(int dist_log_scale,	int dist_resolution,
			int grid_x_resolution, //int grid_dist_resolution,
			int grid_switch_bin_dist, //int grid_switch_x_dist,
			float correction_factor, int decay_smooth, int regularization, int min_dist, int max_dist);
	virtual ~ContactShuffler();
	//virtual long load_contacts(const char* fn, bool symetric);
	virtual long load_contacts(const vector<vector<int> >& contacts, bool symetric);
	virtual int save_contacts(const char* fn, bool symetric, bool with_header);
	virtual int init_obs_decay_from_contacts();
	virtual int	init_exp_decay_from_obs();

	virtual int init_proposal_const();
	virtual int init_proposal_from_area();
	virtual int init_proposal_from_contacts(long proposal_shuffle);
	virtual void correct_proposal_dist();
	virtual int shuffle_contacts(int shuffle_factor, float transition_correction_factor=0.0001,
			float transition_cooling_update=0.01, int debug=0);
	virtual void save_transitions(ostream& out, int id);
	virtual void debug(ostream& out, int id);
	virtual void reset_grid(int grid_x_binsize);
protected:
	virtual void	init_contact_dist_bins();
	virtual int		simple_sample();
	void			update_proposal(int i, int j, int ij, int ji);
	void 			update_proposal_bin(int min_bin, int max_bin);
	virtual int		get_dist_bin(int x, int y);
	virtual float 	get_bin_dist(int bin);
	virtual int		get_grid_bin(int x);
	virtual void 	regularize_decay(vector<float>& decay);
	virtual void 	select_switch_partners(int& contact1, int& contact2, int& grid_index1, int& grid_index2);
	virtual void 	grid_move(int i, int j, int grid_index_i, int grid_index_j);
	virtual void 	print_proposal(ostream& out);

protected:
	long					m_contact_count;
	int						m_max_contact_dist;
	int 					m_dist_log_scale;
	int						m_dist_resolution;
	float					m_log_log_scale;
	vector< vector <int> >	m_contacts;
	vector <int> 			m_contacts_dist_bins;
	vector <int>			m_transitions;
	int						m_grid_x_binsize;
	//int						m_grid_dist_resolution;
	//GenomeGridLog*			m_grid;
	int						m_grid_switch_bin_dist;
	//int						m_grid_switch_x_dist;
	vector< vector< vector<int> > > m_contact_grid;

	vector <float>			m_proposal_freq;
	vector<unsigned long>	m_decay_obs;
	vector <float>			m_decay_exp;
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
