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
#include <cstdio>
#include <charconv>
#include <sys/mman.h>

void* HugePageResource::do_allocate(size_t bytes, size_t alignment) {
	void* p = mmap(NULL, bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	if (p == MAP_FAILED) {
		throw std::bad_alloc();
	}
	madvise(p, bytes, MADV_HUGEPAGE);
	return(p);
}

void HugePageResource::do_deallocate(void* p, size_t bytes, size_t alignment) {
	munmap(p, bytes);
}

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
   m_small_dist(0),
   m_grid_x_binsize(grid_x_resolution),
   //m_grid_dist_resolution(grid_dist_resolution),
   m_grid_switch_bin_dist(grid_switch_bin_dist),
   m_contact_cell(&m_huge_pages),
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

long ContactShuffler::load_contacts(const int* x, const int* y, int stride, long n, bool symetric)
{
	m_contact_count = n;
	m_x.resize(symetric ? n : 2 * n);
	m_y.resize(m_x.size());
	for (long i=0; i<n; i++) {
		m_x[i] = x[i * stride];
		m_y[i] = y[i * stride];
	}
	if (!symetric) {
		cerr << "adding " << m_contact_count << " symmetric contacts" << endl;
		for (long i=0; i<m_contact_count; i++) {
			m_x[m_contact_count + i] = m_y[i];
			m_y[m_contact_count + i] = m_x[i];
		}
	}
	cerr << "finished adding" << endl;
	m_contact_count = m_x.size();
	init_contact_dist_bins();
	m_decay_exp.resize(get_dist_bin(m_min_x,m_max_x)+1, 0);
	m_decay_obs.resize(get_dist_bin(m_min_x,m_max_x)+1, 0);
	m_transitions.resize(get_dist_bin(m_min_x,m_max_x)+1, 0);
	init_obs_decay_from_contacts();
	cerr<<"Loaded "<< m_contact_count << " contacts\n";
	m_reg = log(5)-log(m_contact_count);
	return(m_contact_count);
}

// Writes the same text as streaming the coordinates to an ofstream.
int ContactShuffler::save_contacts(const char* fn, bool symetric, bool with_header) {
	FILE* output = fopen(fn, "w");
	if (output == NULL) {
		cerr << "could not open output file " << fn << endl;
		return(0);
	}
	contacts_from_grid();
	if (with_header)
		fputs("start1\tstart2\n", output);
	vector<char> buf(1 << 20);
	char* p = buf.data();
	char* buf_end = buf.data() + buf.size() - 64;
	for (long i=0; i<m_contact_count; i++) {
		int a = m_x[i];
		int b = m_y[i];
		if (!symetric && a > b) {
			a = m_y[i];
			b = m_x[i];
		}
		p = to_chars(p, buf_end + 64, a).ptr; *p++ = '\t';
		p = to_chars(p, buf_end + 64, b).ptr; *p++ = '\n';
		if (symetric) {
			p = to_chars(p, buf_end + 64, b).ptr; *p++ = '\t';
			p = to_chars(p, buf_end + 64, a).ptr; *p++ = '\n';
		}
		if (p >= buf_end) {
			fwrite(buf.data(), 1, p - buf.data(), output);
			p = buf.data();
		}
	}
	fwrite(buf.data(), 1, p - buf.data(), output);
	fclose(output);
	vector<int>().swap(m_x);
	vector<int>().swap(m_y);
	vector<int>().swap(m_contacts_dist_bins);
	if (symetric) {
		return(m_contact_count * 2);
	}
	return(m_contact_count);
}

int ContactShuffler::init_obs_decay_from_contacts() {
	for (int i=0; i<m_contact_count; i++) {
		m_decay_obs[m_contacts_dist_bins[i]]++;
	}
	//building grid
	int grid_size = floor((m_max_x-m_min_x)/m_grid_x_binsize) + 1;
	cerr << "GRID: " << grid_size << " X " << grid_size << endl;
	build_grid();
	cerr << "finished resizing" << endl;
	cerr << "finished init_obs_decay" << endl;
	return(m_decay_obs.size());
}

// Puts the contacts (m_x, m_y, m_contacts_dist_bins) in their grid cells,
// in index order, then frees those per-index vectors.
void ContactShuffler::build_grid() {
	m_grid_size = floor((m_max_x-m_min_x)/m_grid_x_binsize) + 1;
	size_t cells = (size_t)m_grid_size * m_grid_size;
	m_contact_cell.resize(m_contact_count);
	vector<int> cell_count(cells, 0);
	size_t capacity = 0;
	for (long i=0; i<m_contact_count; i++) {
		m_contact_cell[i] = get_grid_bin(m_x[i]) * m_grid_size + get_grid_bin(m_y[i]);
		cell_count[m_contact_cell[i]]++;
	}
	for (size_t c=0; c<cells; c++) {
		// headroom: cell sizes drift while shuffling
		cell_count[c] += cell_count[c] / 16 + 8;
		capacity += cell_count[c];
	}
	// the cells must go before the pool that holds them
	vector< std::pmr::vector<GridContact> >().swap(m_contact_grid);
	m_grid_pool.reset(new std::pmr::monotonic_buffer_resource(capacity * sizeof(GridContact), &m_huge_pages));
	m_contact_grid.reserve(cells);
	for (size_t c=0; c<cells; c++) {
		m_contact_grid.emplace_back(m_grid_pool.get());
		m_contact_grid[c].reserve(cell_count[c]);
	}
	for (long i=0; i<m_contact_count; i++) {
		GridContact gc = {(int)i, m_x[i], m_y[i], m_contacts_dist_bins[i]};
		m_contact_grid[m_contact_cell[i]].push_back(gc);
	}
	int max_cells = (2*m_grid_switch_bin_dist+1) * (2*m_grid_switch_bin_dist+1);
	m_pool_cumsum.resize(max_cells);
	m_pool_cell.resize(max_cells);
	vector<int>().swap(m_x);
	vector<int>().swap(m_y);
	vector<int>().swap(m_contacts_dist_bins);
}

// Fills m_x, m_y, m_contacts_dist_bins from the grid.
void ContactShuffler::contacts_from_grid() {
	m_x.resize(m_contact_count);
	m_y.resize(m_contact_count);
	m_contacts_dist_bins.resize(m_contact_count);
	for (size_t c=0; c<m_contact_grid.size(); c++) {
		for (const GridContact& gc : m_contact_grid[c]) {
			m_x[gc.idx] = gc.x;
			m_y[gc.idx] = gc.y;
			m_contacts_dist_bins[gc.idx] = gc.dist_bin;
		}
	}
}

int	ContactShuffler::init_exp_decay_from_obs() {
	m_decay_exp.resize(m_decay_obs.size());
	vector<unsigned long> int_exp_decay(m_decay_obs);
	VectorUtils::smooth_vector(int_exp_decay, m_decay_exp, m_decay_smooth); //exp decay is not in log format
	VectorUtils::log_vec(m_decay_exp, m_decay_exp);
	regularize_decay(m_decay_exp);
	m_decay_exp_nonzero.resize(m_decay_exp.size());
	for (unsigned bin=0; bin<m_decay_exp.size(); bin++) {
		m_decay_exp_nonzero[bin] = exp(m_decay_exp[bin]) != 0;
	}
	m_log_obs.assign(m_decay_exp.size(), 0);
	m_log_obs_count.assign(m_decay_exp.size(), ULONG_MAX);
	return(m_decay_exp.size());
}

void	ContactShuffler::init_contact_dist_bins() {
	m_contacts_dist_bins.resize(m_contact_count, -1);
	m_min_x = INT_MAX;
	m_max_x = 0;
	cerr << "init_contact_dist_bins" << endl;
	for (int i=0; i<m_contact_count; i++) {
		if (m_x[i] < m_min_x) m_min_x = m_x[i];
		if (m_y[i] < m_min_x) m_min_x = m_y[i];
		if (m_x[i] > m_max_x) m_max_x = m_x[i];
		if (m_y[i] > m_max_x) m_max_x = m_y[i];
	}
	init_dist_bin_table();
	for (int i=0; i<m_contact_count; i++) {
		m_contacts_dist_bins[i] = get_dist_bin(m_x[i], m_y[i]);
	}
	cerr << "m_min_x=" << m_min_x << endl;
	cerr << "m_max_x=" << m_max_x << endl;
}

int ContactShuffler::simple_sample() {
	int i;
	int cell_i, cell_j, grid_index_i, grid_index_j;
	i = floor(Random::fraction_truncated() * m_contact_count);
	// The next sample starts with the draw 3 to 6 draws from now (one for the
	// member of the cell, one per partner try, maybe one for acceptance):
	// prefetch the cells of those contacts. Prefetches do not change results.
	for (int k=3; k<=6; k++) {
		int next = floor(Random::peek_fraction(k) * m_contact_count);
		__builtin_prefetch(m_contact_cell.data() + next);
	}

	select_switch_partners(i, cell_i, grid_index_i, cell_j, grid_index_j);
	const GridContact& ci = m_contact_grid[cell_i][grid_index_i];
	const GridContact& cj = m_contact_grid[cell_j][grid_index_j];

	int dist_i_bin = ci.dist_bin;
	int dist_j_bin = cj.dist_bin;
	int dist_ij_bin = get_dist_bin(ci.x, cj.y);
	int dist_ji_bin = get_dist_bin(cj.x, ci.y);

	if (dist_ij_bin < 0 || dist_ji_bin < 0) {
		return(0);
	}
	float log_acceptance = m_decay_exp[dist_ij_bin] + m_decay_exp[dist_ji_bin]-
			m_decay_exp[dist_i_bin]-m_decay_exp[dist_j_bin] +
			m_proposal_freq[dist_i_bin] + m_proposal_freq[dist_j_bin] -
			m_proposal_freq[dist_ij_bin]- m_proposal_freq[dist_ji_bin];
	// exp(x) > 1 for x > 1, and exp(x) is exactly 0 in float for x < -104:
	// skip the exp (it over/underflows slowly) for moves out of or into zeroed bins
	float acceptance_prob = log_acceptance > 1 ? log_acceptance :
			(log_acceptance < -104 ? 0.0f : exp(log_acceptance));

	if (acceptance_prob > 1 || Random::fraction() < acceptance_prob) {
		grid_move(cell_i, grid_index_i, cell_j, grid_index_j, dist_ij_bin, dist_ji_bin);
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
	int dist_bin;
	if (dist < m_small_dist) {
		dist_bin = m_small_dist_bin[dist];
	} else if (!m_dist_sub_range.empty()) {
		int octave = 31 - __builtin_clz(dist);
		const DistSubRange& r = m_dist_sub_range[((octave - 12) << 10) | ((dist >> (octave - 10)) & 1023)];
		dist_bin = r.bin + (dist >= r.next_bin_dist);
	} else {
		dist_bin = dist_bin_formula(dist);
	}
	if (dist_bin < 0) {
		return(-1);
	}
	return(dist_bin);
}

int	ContactShuffler::dist_bin_formula(int dist) {
	int dist_bin = floor(m_dist_resolution * log(1+dist) / m_log_log_scale) - m_min_dist;
	return(dist_bin);
}

// Tabulates dist_bin_formula for 0 <= dist <= m_max_x - m_min_x. The formula is
// non-decreasing in dist (log(1+dist) of consecutive integers differ by ~1/dist,
// far more than its rounding error). Distances >= 4096 are split into 1024
// sub-ranges per octave; each is checked to span at most two consecutive bins,
// and the first distance of the second bin is found by bisection, so the lookup
// returns exactly what the formula does. If a check fails, get_dist_bin keeps
// using the formula.
void ContactShuffler::init_dist_bin_table() {
	int max_dist = m_max_x - m_min_x;
	m_small_dist = min(4096, max_dist + 1);
	m_small_dist_bin.resize(m_small_dist);
	for (int dist=0; dist<m_small_dist; dist++) {
		m_small_dist_bin[dist] = dist_bin_formula(dist);
	}
	m_dist_sub_range.clear();
	if (max_dist < 4096) {
		return;
	}
	int max_octave = 31 - __builtin_clz(max_dist);
	vector<DistSubRange> table((size_t)(max_octave - 11) << 10);
	for (int octave=12; octave<=max_octave; octave++) {
		for (int sub=0; sub<1024; sub++) {
			long lo = (long)(1024 + sub) << (octave - 10);
			long hi = lo + (1L << (octave - 10)) - 1;
			if (lo > max_dist) {
				break;
			}
			if (hi > max_dist) {
				hi = max_dist;
			}
			DistSubRange& r = table[((octave - 12) << 10) | sub];
			r.bin = dist_bin_formula(lo);
			r.next_bin_dist = INT_MAX;
			int hi_bin = dist_bin_formula(hi);
			if (hi_bin == r.bin + 1) {
				long a = lo, b = hi;	// dist_bin_formula(a) == r.bin, dist_bin_formula(b) == r.bin + 1
				while (b - a > 1) {
					long mid = (a + b) / 2;
					if (dist_bin_formula(mid) > r.bin) b = mid; else a = mid;
				}
				r.next_bin_dist = b;
			} else if (hi_bin != r.bin) {
				cerr << "distance bin table not used (sub-range " << lo << "-" << hi << ")" << endl;
				return;
			}
		}
	}
	m_dist_sub_range.swap(table);
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
		if (m_decay_exp_nonzero[bin]) {
			if (m_decay_obs[bin] != 0) {
				// log(m_decay_obs[bin]), cached: consecutive calls often see the same counts
				if (m_log_obs_count[bin] != m_decay_obs[bin]) {
					m_log_obs[bin] = log(m_decay_obs[bin]);
					m_log_obs_count[bin] = m_decay_obs[bin];
				}
				m_proposal_freq[bin] += m_log_obs[bin] - m_decay_exp[bin] + m_correction_factor;
			} else {
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
		int i;
		int cell_i, cell_j, grid_index_i, grid_index_j;
		i = floor(Random::fraction_truncated() * m_contact_count);
		select_switch_partners(i, cell_i, grid_index_i, cell_j, grid_index_j);
		const GridContact& ci = m_contact_grid[cell_i][grid_index_i];
		const GridContact& cj = m_contact_grid[cell_j][grid_index_j];

		int dist_ij_bin = get_dist_bin(ci.x, cj.y);
		int dist_ji_bin = get_dist_bin(ci.y, cj.x);
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
	cerr << "shuffling..." << shuffle_factor << " iterations"<< endl;

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

// Draws a random contact from the grid cell of `contact`, then a partner from
// the (2 * m_grid_switch_bin_dist + 1)^2 surrounding cells whose coordinates are
// both within m_grid_switch_bin_dist * m_grid_x_binsize of it (rejection sampling).
// Returns each one as (grid cell, index in cell).
void	ContactShuffler::select_switch_partners(int contact, int& cell1, int& grid_index1, int& cell2, int& grid_index2) {
	cell1 = m_contact_cell[contact];
	int s1 = cell1 / m_grid_size;
	int s2 = cell1 % m_grid_size;
	int max_dist = m_grid_switch_bin_dist*m_grid_x_binsize;
	//selecting random member from grid s1, s2
	grid_index1 = floor(Random::fraction_truncated() * m_contact_grid[cell1].size());
	const GridContact& c1 = m_contact_grid[cell1][grid_index1];
	__builtin_prefetch(&c1);

	//building cumsum vector
	int* cumsum = m_pool_cumsum.data();
	int* pool_cell = m_pool_cell.data();
	int pool_cells = 0;
	int total_pool=0;
	int d1 = s1-m_grid_switch_bin_dist < 0 ? 0 : s1-m_grid_switch_bin_dist;
	int d2;
	for (;d1 <= s1+m_grid_switch_bin_dist && d1<m_grid_size; d1++) {
		d2 = s2-m_grid_switch_bin_dist < 0 ? 0 : s2-m_grid_switch_bin_dist;
		for (;d2 <= s2+m_grid_switch_bin_dist && d2<m_grid_size; d2++) {
			int cell = d1 * m_grid_size + d2;
			total_pool += m_contact_grid[cell].size();
			cumsum[pool_cells] = total_pool-1;
			pool_cell[pool_cells] = cell;
			pool_cells++;
		}
	}
	while (true) {
		grid_index2 = floor(Random::fraction_truncated() * total_pool);
		// first cell whose cumsum >= grid_index2 (lower_bound)
		int i = 0;
		while (i < pool_cells - 1 && cumsum[i] < grid_index2) i++;
		if (i> 0) grid_index2 = grid_index2 - cumsum[i-1] - 1;

		//reaching this point, we have our grid bin (d1,d2)
		const GridContact& c2 = m_contact_grid[pool_cell[i]][grid_index2];
		// prefetch the partner the next try would draw
		int next = floor(Random::peek_fraction(1) * total_pool);
		int k = 0;
		while (k < pool_cells - 1 && cumsum[k] < next) k++;
		if (k > 0) next = next - cumsum[k-1] - 1;
		__builtin_prefetch(m_contact_grid[pool_cell[k]].data() + next);
		if (abs(c2.x - c1.x) < max_dist &&
			abs(c2.y - c1.y) < max_dist) {
			cell2 = pool_cell[i];
			return;
		}
	}
}

// Swaps the second coordinates of the two contacts (the same moves in the
// same order as the index-based grid, so cell contents keep the same order).
void ContactShuffler::grid_move(int cell_i, int grid_index_i, int cell_j, int grid_index_j, int dist_ij_bin, int dist_ji_bin) {
	GridContact new_i = m_contact_grid[cell_i][grid_index_i];
	GridContact new_j = m_contact_grid[cell_j][grid_index_j];
	new_i.y = m_contact_grid[cell_j][grid_index_j].y;
	new_i.dist_bin = dist_ij_bin;
	new_j.y = m_contact_grid[cell_i][grid_index_i].y;
	new_j.dist_bin = dist_ji_bin;
	int bin_i_0 = cell_i / m_grid_size;
	int bin_i_1 = cell_i % m_grid_size;
	int bin_j_0 = cell_j / m_grid_size;
	int bin_j_1 = cell_j % m_grid_size;
	if (bin_i_1 != bin_j_1) {
		std::pmr::vector<GridContact>& grid_i = m_contact_grid[cell_i];
		grid_i[grid_index_i] = grid_i.back();
		grid_i.pop_back();
		int to_i = bin_i_0 * m_grid_size + bin_j_1;
		m_contact_grid[to_i].push_back(new_i);
		m_contact_cell[new_i.idx] = to_i;

		std::pmr::vector<GridContact>& grid_j = m_contact_grid[cell_j];
		grid_j[grid_index_j] = grid_j.back();
		grid_j.pop_back();
		int to_j = bin_j_0 * m_grid_size + bin_i_1;
		m_contact_grid[to_j].push_back(new_j);
		m_contact_cell[new_j.idx] = to_j;
	} else {
		m_contact_grid[cell_i][grid_index_i] = new_i;
		m_contact_grid[cell_j][grid_index_j] = new_j;
	}
}

void ContactShuffler::save_transitions(ostream& out, int id) {
	int max_bins = m_transitions.size();
	for (int bin=0; bin<max_bins; bin++) {
		out << id << "\t" << get_bin_dist(bin) << "\t" << m_transitions[bin] << endl;
	}

}

void ContactShuffler::reset_grid(int grid_x_binsize) {
	cerr << "resetting grid from " << m_grid_x_binsize << " to " << grid_x_binsize << endl;
	contacts_from_grid();
	m_grid_x_binsize = grid_x_binsize;
	build_grid();
}

void ContactShuffler::print_proposal(ostream& out) {
	int max_bins = m_proposal_freq.size();
	for (int bin=0; bin<max_bins; bin++) {
			out << "\t" << m_proposal_freq[bin];
	}
}
