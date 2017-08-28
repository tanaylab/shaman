/*
 * GenomeGridLog.cpp
 *
 *  Created on: Dec 1, 2016
 *      Author: nettam
 */

#include "GenomeGridLog.h"
#include <cmath>
#include <stdlib.h>
#include <iostream>
using namespace std;

GenomeGridLog::GenomeGridLog(int dist_resolution, int log_scale, int x_binsize,
		int min_x, int max_x, int min_dist):
	m_dist_resolution(dist_resolution),
	m_log_scale(log_scale),
	m_x_binsize(x_binsize),
	m_min_x(min_x),
	m_max_x(max_x),
	m_log_log_scale(log(log_scale))
{
	int max_dist = m_max_x-m_min_x;
	m_min_dist = (m_dist_resolution * log(min_dist+1)/m_log_log_scale);
	get_grid_bin(m_min_x, m_max_x, m_dim1_size, m_dim2_size);
	m_dim1_size = floor(float(m_max_x)/m_x_binsize) + 1;
	m_dim2_size = floor(m_dist_resolution * log(1+max_dist) / m_log_log_scale) - m_min_dist + 1;

	cerr << "GenomeGrid:: " << m_dim1_size << " X " << m_dim2_size << endl <<
			"dist resolution = " << m_dist_resolution << endl <<
			"log_scale = " << m_log_scale << endl <<
			"x bin size = " << m_x_binsize << endl <<
			"min_x = " << m_min_x << endl <<
			"max_x = " << m_max_x << endl <<
			"min dist = " << min_dist << endl <<
			"log_log_scale = " << m_log_log_scale << endl;

}

GenomeGridLog::~GenomeGridLog() {
	// TODO Auto-generated destructor stub
}


void GenomeGridLog::get_grid_bin(int x, int y, int& d1, int& d2) {
	d1 = (0.5 * (x+y)) / m_x_binsize;
	int dist = abs(x - y);
	d2 = floor(m_dist_resolution * log(1+dist) / m_log_log_scale) - m_min_dist;
}
