/*
 * GenomeGridLog.h
 *
 *  Created on: Dec 1, 2016
 *      Author: nettam
 */

#ifndef GENOMEGRIDLOG_H_
#define GENOMEGRIDLOG_H_

class GenomeGridLog {
public:
	GenomeGridLog(int dist_resolution, int log_scale, int x_binsize,
			int min_x, int max_x, int min_dist);
	virtual ~GenomeGridLog();
	int get_dim1_size() { return m_dim1_size; }
	int get_dim2_size() { return m_dim2_size; }
	void get_grid_bin(int x, int y, int& d1, int& d2);
private:
	int m_dist_resolution;
	int m_log_scale;
	int m_x_binsize;
	int m_min_x;
	int m_max_x;
	int m_min_dist;
	int m_dim1_size;
	int m_dim2_size;
	float m_log_log_scale;

};

#endif /* GENOMEGRIDLOG_H_ */
