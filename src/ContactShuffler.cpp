/*
 * ContactShuffler.cpp
 *
 *  Created on: Nov 30, 2016
 *      Author: nettam
 */

#include "ContactShuffler.h"
#include "GenomeGridLog.h"
#include "Parser.h"
#include "VectorUtils.h"
#include "MathUtils.h"
#include "Random.h"
#include "macro.h"
#include <cmath>
#include <climits>
#include <fstream>

ContactShuffler::ContactShuffler(int dist_log_scale, int dist_resolution,
		int grid_x_resolution,
		int grid_switch_bin_dist,
		float correction_factor,
		int decay_smooth, int regularization, int min_dist, int max_dist):
   m_contact_count(0),
   m_max_contact_dist(max_dist),
   m_dist_log_scale(dist_log_scale),
   m_dist_resolution(dist_resolution),
   m_log_log_scale(log(dist_log_scale)),
   m_grid_x_binsize(grid_x_resolution),
   //m_grid_dist_resolution(grid_dist_resolution),
   m_grid_switch_bin_dist(grid_switch_bin_dist),
   //m_grid_switch_x_dist(grid_switch_x_dist),
   m_correction_factor(log(correction_factor)),
   m_decay_smooth(decay_smooth),
   m_regularization(regularization)
{
	m_min_dist = (m_dist_resolution * log(min_dist)/m_log_log_scale);
	if (m_min_dist < 0) m_min_dist=0;
}

ContactShuffler::~ContactShuffler() {
	// TODO Auto-generated destructor stub
	//delete(m_grid);
}

long ContactShuffler::load_contacts(const vector<vector<int> >& contacts, bool symetric)
{
	m_contacts = contacts;
	m_contact_count = m_contacts.size();
	if (!symetric) {
		cerr << "adding " << m_contact_count << " symmetric contacts" << endl;
		for (int i=0; i<m_contact_count; i++) {
			m_contacts.resize(m_contact_count+i+1, vector<int>(2));
			m_contacts[m_contact_count + i][0] = m_contacts[i][1];
			m_contacts[m_contact_count + i][1] = m_contacts[i][0];
		}
	}
	cerr << "finished adding" << endl;
	m_contact_count = m_contacts.size();
	init_contact_dist_bins();
	m_decay_exp.resize(get_dist_bin(m_min_x,m_max_x)+1, 0);
	m_decay_obs.resize(get_dist_bin(m_min_x,m_max_x)+1, 0);
	m_transitions.resize(get_dist_bin(m_min_x,m_max_x)+1, 0);
	init_obs_decay_from_contacts();
	cerr<<"Loaded "<< m_contact_count << " contacts\n";
	m_reg = log(5)-log(m_contact_count);
	return(m_contact_count);
}

int ContactShuffler::save_contacts(const char* fn, bool symetric, bool with_header) {
	ofstream output(fn);
	if (!output.is_open()) {
		cerr << "could not open output file " << fn << endl;
	}
	//header
	if (with_header)
		output << "start1\tstart2" << endl;
	for ( int i=0; i<m_contact_count; i++) {
		if (symetric) {
			output << m_contacts[i][0] << "\t" << m_contacts[i][1] << "\n";
			output << m_contacts[i][1] << "\t" << m_contacts[i][0] << "\n";
		} else {
			if (m_contacts[i][0]  < m_contacts[i][1]){
				output << m_contacts[i][0] << "\t" << m_contacts[i][1] << "\n";
			} else {
				output << m_contacts[i][1] << "\t" << m_contacts[i][0] << "\n";
			}
		}
	}
	if (symetric) {
		return(m_contact_count * 2);
	}
	return(m_contact_count);
}

int ContactShuffler::init_obs_decay_from_contacts() {
	//building grid
	int grid_size = floor((m_max_x-m_min_x)/m_grid_x_binsize) + 1;
	cerr << "GRID: " << grid_size << " X " << grid_size << endl;
	m_contact_grid.resize(grid_size,
				vector< vector<int> >(grid_size, vector< int >(0)));

	cerr << "finished resizing" << endl;
	for (int i=0; i<m_contact_count; i++) {
		m_decay_obs[m_contacts_dist_bins[i]]++;
		m_contact_grid[get_grid_bin(m_contacts[i][0])][get_grid_bin(m_contacts[i][1])].push_back(i);
	}
	cerr << "finished init_obs_decay" << endl;
	return(m_decay_obs.size());
}

int	ContactShuffler::init_exp_decay_from_obs() {
	m_decay_exp.resize(m_decay_obs.size());
	vector<unsigned long> int_exp_decay(m_decay_obs);
	VectorUtils::smooth_vector(int_exp_decay, m_decay_exp, m_decay_smooth); //exp decay is not in log format
	VectorUtils::log_vec(m_decay_exp, m_decay_exp);
	regularize_decay(m_decay_exp);
	return(m_decay_exp.size());
}

void	ContactShuffler::init_contact_dist_bins() {
	m_contacts_dist_bins.resize(m_contact_count, -1);
	m_min_x = INT_MAX;
	m_max_x = 0;
	cerr << "init_contact_dist_bins" << endl;
	for (int i=0; i<m_contact_count; i++) {
		m_contacts_dist_bins[i] = get_dist_bin(m_contacts[i][0], m_contacts[i][1]);
		if (m_contacts[i][0] < m_min_x) m_min_x = m_contacts[i][0];
		if (m_contacts[i][1] < m_min_x) m_min_x = m_contacts[i][1];
		if (m_contacts[i][0] > m_max_x) m_max_x = m_contacts[i][0];
		if (m_contacts[i][1] > m_max_x) m_max_x = m_contacts[i][1];
	}
	cerr << "m_min_x=" << m_min_x << endl;
	cerr << "m_max_x=" << m_max_x << endl;
}

int ContactShuffler::simple_sample() {
	//cerr << "simple sample start" << endl;
	int i,j;
	int grid_index_i, grid_index_j;
	i = floor(Random::fraction_truncated() * m_contact_count);

	//j = floor(Random::fraction_truncated() * m_contact_count);
	select_switch_partners(i, j, grid_index_i, grid_index_j);

	int dist_i_bin = m_contacts_dist_bins[i];
	int dist_j_bin = m_contacts_dist_bins[j];
	int dist_ij_bin = get_dist_bin(m_contacts[i][0], m_contacts[j][1]);
	int dist_ji_bin = get_dist_bin(m_contacts[j][0], m_contacts[i][1]);

	if (dist_ij_bin < 0 || dist_ji_bin < 0) {
		return(0);
	}
	float acceptance_prob = exp(m_decay_exp[dist_ij_bin] + m_decay_exp[dist_ji_bin]-
			m_decay_exp[dist_i_bin]-m_decay_exp[dist_j_bin] +
			m_proposal_freq[dist_i_bin] + m_proposal_freq[dist_j_bin] -
			m_proposal_freq[dist_ij_bin]- m_proposal_freq[dist_ji_bin]);

	if (acceptance_prob > 1 || Random::fraction() < acceptance_prob) {
		grid_move(i, j, grid_index_i, grid_index_j);
		int tmp = m_contacts[i][1];
		m_contacts[i][1] = m_contacts[j][1];
		m_contacts[j][1] = tmp;
		m_contacts_dist_bins[i] = dist_ij_bin;
		m_contacts_dist_bins[j] = dist_ji_bin;
		m_decay_obs[dist_i_bin]--;
		m_decay_obs[dist_j_bin]--;
		m_decay_obs[dist_ij_bin]++;
		m_decay_obs[dist_ji_bin]++;
		m_transitions[dist_i_bin]++;
		m_transitions[dist_j_bin]++;
		return(1);
	}
	return(0);
}


int	ContactShuffler::get_dist_bin(int x, int y) {
	int dist = abs(x - y);
	int dist_bin = floor(m_dist_resolution * log(1+dist) / m_log_log_scale) - m_min_dist;
	if (dist_bin < 0) {
		return(-1);
	}
	return(dist_bin);
}

float ContactShuffler::get_bin_dist(int bin) {
	return((float)(bin+m_min_dist)/m_dist_resolution);
}

int ContactShuffler::get_grid_bin(int x) {
	return(floor((x-m_min_x)/(float)m_grid_x_binsize));
}

void ContactShuffler::correct_proposal_dist() {
	float sum_proposal=FLT_MIN_EXP;
	int max_bins = m_proposal_freq.size();
	//vector<float> smooth_obs(max_bins, FLT_MIN_EXP);
	//VectorUtils::smooth_vector(m_decay_obs, smooth_obs, m_decay_smooth);
	//VectorUtils::log_vec(smooth_obs, smooth_obs);
	//regularize_decay(smooth_obs);

	for (int bin=0; bin<max_bins; bin++) {
		if (exp(m_decay_exp[bin]) != 0) {
			if (m_decay_obs[bin] != 0)
			//if (smooth_obs[bin] != 0)
				//m_proposal_freq[bin] += smooth_obs[bin] - m_decay_exp[bin] + m_correction_factor;
				m_proposal_freq[bin] += log(m_decay_obs[bin]) - m_decay_exp[bin] + m_correction_factor;
			else {
				m_proposal_freq[bin] = FLT_MIN_EXP;
			}
		} else {
			m_proposal_freq[bin] = FLT_MIN_EXP;
		}
	    log_sum_log(sum_proposal,m_proposal_freq[bin]);
	}

	for (int bin=0; bin<max_bins; bin++) {
			m_proposal_freq[bin] -= sum_proposal;
	}

}

int ContactShuffler::init_proposal_const() {
	int bins = m_decay_exp.size();
	m_proposal_freq.resize(bins, -log(bins));
	//m_proposal_freq[0] = FLT_MIN_EXP;
	//m_proposal_freq[bins-1] = FLT_MIN_EXP;
	return(0);
}

int ContactShuffler::init_proposal_from_area() {
	m_proposal_freq.resize(m_decay_exp.size(), FLT_MIN_EXP);
	for (unsigned d=0; d<m_decay_exp.size(); d++) {
		float f = pow(m_dist_log_scale, (1.0+d + m_min_dist)/(float)m_dist_resolution);
		float min_clip_factor = 1-f/m_max_contact_dist;
		if (min_clip_factor < 0) {
			min_clip_factor=0;
		}
		f = pow(m_dist_log_scale, (d+m_min_dist)/(float)m_dist_resolution);
		float max_clip_factor = 1-f/m_max_contact_dist;
		if (max_clip_factor < 0) {
			max_clip_factor=0;
		}
		m_proposal_freq[d] = max_clip_factor+min_clip_factor==0 ? FLT_MIN_EXP : log(((max_clip_factor + min_clip_factor)/2) * f);
	}
	float sum_proposal=FLT_MIN_EXP;
	for (int bin=0; bin<m_decay_exp.size(); bin++) {
		log_sum_log(sum_proposal,m_proposal_freq[bin]);
	}
	for (int bin=0; bin<m_decay_exp.size(); bin++) {
		m_proposal_freq[bin] -= sum_proposal;
	}
	return(0);
}

int ContactShuffler::init_proposal_from_contacts(long proposal_shuffle) {
	cerr << "init proposal from contacts " << proposal_shuffle << " iterations" << endl;
	m_proposal_shuffle = proposal_shuffle;
	int max_bins = m_decay_exp.size();
	vector<unsigned long> proposal_count(max_bins, 0);
	vector<unsigned long> proposal_suggest(max_bins, 0);
	m_proposal_freq.resize(m_decay_exp.size(), 0);
	for (long iter=0; iter<proposal_shuffle; iter++) {
		int i,j;
		int grid_index_i, grid_index_j;
		i = floor(Random::fraction_truncated() * m_contact_count);
		select_switch_partners(i, j, grid_index_i, grid_index_j);

		int dist_ij_bin = get_dist_bin(m_contacts[i][0], m_contacts[j][1]);
		int dist_ji_bin = get_dist_bin(m_contacts[i][1], m_contacts[j][0]);
		//log_sum_log(m_proposal_freq[dist_ij_bin], 0);
		//log_sum_log(m_proposal_freq[dist_ji_bin], 0);
		//proposal_suggest[m_contacts_dist_bins[i]]++;
		//proposal_suggest[m_contacts_dist_bins[j]]++;
		if (dist_ij_bin >= 0) proposal_count[dist_ij_bin]++;
		if (dist_ji_bin >= 0) proposal_count[dist_ji_bin]++;
	}
	/*
	for (int bin=0; bin<max_bins; bin++) {
		if (proposal_suggest[bin] > 0) {
			m_proposal_freq[bin] = ((float)proposal_count[bin])/proposal_suggest[bin];
		} else {
			m_proposal_freq[bin] = 10;
		}
	}
	*/
	VectorUtils::smooth_vector(proposal_count, m_proposal_freq, m_decay_smooth);
	VectorUtils::log_vec(m_proposal_freq, m_proposal_freq);
	float sum_proposal=FLT_MIN_EXP;
	for (int bin=0; bin<max_bins; bin++) {
		log_sum_log(sum_proposal,m_proposal_freq[bin]);
	}
	for (int bin=0; bin<max_bins; bin++) {
		m_proposal_freq[bin] -= sum_proposal;
	}
	return(1);
}



void ContactShuffler::regularize_decay(vector<float>& decay) {
 //regularization

 for (unsigned int bin=0; bin<decay.size(); bin++) {
	 log_sum_log(decay[bin],m_reg);
	  if (exp(decay[bin]) < m_regularization)
		decay[bin] = FLT_MIN_EXP;
 }
}

void ContactShuffler::debug(ostream& out, int id) {
	int max_bins = m_decay_exp.size();
	for (int bin=0; bin<max_bins; bin++) {
		out << id << "\t" << get_bin_dist(bin) << "\t" << m_decay_exp[bin]
		    << "\t" << m_decay_obs[bin] << "\t" << m_proposal_freq[bin] << endl;
	}
}

int ContactShuffler::shuffle_contacts(int shuffle_factor, float transition_correction_factor,
		float transition_cooling_update, int debug) {
	cerr << "shuffling..." << shuffle_factor << " iterations debug=" << debug << endl;

	int transitions_per_correction = floor(m_contact_count*0.0001);
	int percentile = floor((m_contact_count)/10);
	if (debug) {
	  cout << m_grid_x_binsize << "\t0\t0";
	  print_proposal(cout);
	  cout << endl;
	}
	long total_samples=0;
	long total_transitions=0;
	for (int i=0; i<shuffle_factor; i++) {
		//transitions_per_correction = floor(m_contact_count*0.0001 * (i+1));

		long transitions=0;
		long samples=0;
		while (transitions < m_contact_count) {
			//while (transitions < shuffle_factor) {
			samples++;
			transitions += simple_sample();
			if (samples > 1000 && transitions == 0) {
				cerr << "not making any transitions... stopping early" << endl;
				return(0);
			}

			if (transitions % transitions_per_correction == 0) {
				correct_proposal_dist();
				//cerr << transitions_per_correction << " --> ";
				//transitions_per_correction += floor(m_contact_count*transition_correction_factor * 0.02);
				//cerr << transitions_per_correction << endl;
			    if (debug) {
			    	cout << m_grid_x_binsize << "\t" << (total_samples + samples) << "\t" <<
			    			floor(100*(total_transitions + transitions))/(total_samples+samples)/100;
			    	print_proposal(cout);
			    	cout << endl;
			    }
		    }
		    if (transitions % percentile == 0) {
			  cerr << i << " :: " << samples << "\t" << transitions << "\t" << floor(100*transitions/samples)/100
					<< "\t" << transitions/percentile << endl;
		    }
		}
		total_transitions += transitions;
		total_samples += samples;
	}
	return(0);

}

//need to make sure that grid_index1 is always >= grid_index2
void	ContactShuffler::select_switch_partners(int& contact1, int& contact2, int& grid_index1, int& grid_index2) {

	//cerr << "switch_partner_start: " << m_contacts[contact1][0]  << " - "
	//		<< m_contacts[contact1][1] << endl;
	int s1 = get_grid_bin(m_contacts[contact1][0]);
	int s2 = get_grid_bin(m_contacts[contact1][1]);
	int max_dist = m_grid_switch_bin_dist*m_grid_x_binsize;
	//selecting random member from grid s1, s2
	grid_index1 = floor(Random::fraction_truncated() * m_contact_grid[s1][s2].size());
	contact1 = m_contact_grid[s1][s2][grid_index1];
	//cerr << "selected contact: " << m_contacts[contact1][0]  << " - "
	//	<< m_contacts[contact1][1] << ", grid bin:[ " << s1 << " , " << s2 << " ]" << endl;

	//building cumsum vector
	vector<int> cumsum;
	vector<int> bin1;
	vector<int> bin2;
	int total_pool=0;
	int d1 = s1-m_grid_switch_bin_dist < 0 ? 0 : s1-m_grid_switch_bin_dist;
	int d2 = s2-m_grid_switch_bin_dist < 0 ? 0 : s2-m_grid_switch_bin_dist ;
	//cerr << d1 << "\t" << d2 << "\tcontact grid size = " << m_contact_grid.size() << endl;
	for (;d1 <= s1+m_grid_switch_bin_dist && d1<m_contact_grid.size(); d1++) {
		d2 = s2-m_grid_switch_bin_dist < 0 ? 0 : s2-m_grid_switch_bin_dist;
		for (;d2 <= s2+m_grid_switch_bin_dist && d2<m_contact_grid.size(); d2++) {
			//cerr << "\t" << d1 << "\t" << d2 << "\t" << m_contact_grid[d1][d2].size() << endl;
			total_pool += m_contact_grid[d1][d2].size();
			cumsum.push_back(total_pool-1);
			bin1.push_back(d1);
			bin2.push_back(d2);
		}
	}
	bool valid_choice=0;
	int number_of_tries=0;
	while (!valid_choice) {
	  ++number_of_tries;
      grid_index2 = floor(Random::fraction_truncated() * total_pool);
      vector<int>::iterator find = lower_bound(cumsum.begin(), cumsum.end(), grid_index2);
	  int i=distance(cumsum.begin(), find);
	  //if (i> 0) grid_index2 = cumsum[i] - cumsum[i-1] - 1;
	  if (i> 0) grid_index2 = grid_index2 - cumsum[i-1] - 1;

	  //reaching this point, we have our grid bin (d1,d2)
	  contact2 = m_contact_grid[bin1[i]][bin2[i]][grid_index2];
	  //cerr << "\t" << m_contacts[contact2][0]
	  //     << "-" << m_contacts[contact2][1]<< " in grid bin [ " << bin1[i]
	  //     << " , " << bin2[i] << " ] " << endl;
	  if (abs(m_contacts[contact2][0] - m_contacts[contact1][0]) < max_dist &&
		abs(m_contacts[contact2][1] - m_contacts[contact1][1]) < max_dist) {
			valid_choice=1;
		}
	}
//	cerr << "switch_partner_end" << endl;
	return;
}

void ContactShuffler::grid_move(int i, int j, int grid_index_i, int grid_index_j) {
	//cerr << "GRID MOVE START" << endl;
	int bin_i_0 = get_grid_bin(m_contacts[i][0]);
	int bin_i_1 = get_grid_bin(m_contacts[i][1]);
	int bin_j_0 = get_grid_bin(m_contacts[j][0]);
	int bin_j_1 = get_grid_bin(m_contacts[j][1]);
	if (bin_i_1 != bin_j_1) {
		m_contact_grid[bin_i_0][bin_i_1][grid_index_i] = m_contact_grid[bin_i_0][bin_i_1].back();
		m_contact_grid[bin_i_0][bin_i_1].pop_back();
		m_contact_grid[bin_i_0][bin_j_1].push_back(i);

		m_contact_grid[bin_j_0][bin_j_1][grid_index_j] = m_contact_grid[bin_j_0][bin_j_1].back();
		m_contact_grid[bin_j_0][bin_j_1].pop_back();
		m_contact_grid[bin_j_0][bin_i_1].push_back(j);
	}

	//cerr << "GRID MOVE END" << endl;
}

void ContactShuffler::save_transitions(ostream& out, int id) {
	int max_bins = m_transitions.size();
	for (int bin=0; bin<max_bins; bin++) {
		out << id << "\t" << get_bin_dist(bin) << "\t" << m_transitions[bin] << endl;
	}

}

void ContactShuffler::reset_grid(int grid_x_binsize) {
	cerr << "resetting grid from " << m_grid_x_binsize << " to " << grid_x_binsize << endl;
	m_grid_x_binsize = grid_x_binsize;
    int grid_size = floor((m_max_x-m_min_x)/m_grid_x_binsize) + 1;
    m_contact_grid.resize(0, vector< vector<int> >(0, vector< int> (0)) );
    m_contact_grid.resize(grid_size,
			vector< vector<int> >(grid_size,
					vector< int >(0)));
    for (int i=0; i<m_contact_count; i++) {
		m_contact_grid[get_grid_bin(m_contacts[i][0])][get_grid_bin(m_contacts[i][1])].push_back(i);
    }
}

void ContactShuffler::print_proposal(ostream& out) {
	int max_bins = m_proposal_freq.size();
	for (int bin=0; bin<max_bins; bin++) {
			out << "\t" << m_proposal_freq[bin];
	}
}
