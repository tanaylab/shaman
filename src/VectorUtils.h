/*
 * VectorUtils.h
 *
 *  Created on: Nov 30, 2016
 *      Author: nettam
 */

#ifndef VECTORUTILS_H_
#define VECTORUTILS_H_
#include <vector>
#include <map>
#include <iostream>
#include <deque>
#include <algorithm>
#include <cfloat>
#include <cmath>
using namespace std;


class VectorUtils {
public:
	template <class T>
	static void log_vec(const vector<T>& vec, vector<float>& log_vec) {
		for (unsigned int i=0; i<vec.size(); i++) {
			log_vec[i] = (vec[i]<=0 ? FLT_MIN_EXP : log(vec[i]));
		}
	};

	static void exp_vec(const vector<float>& vec, vector<float>& exp_vec) {
		for (unsigned int i=0; i<vec.size(); i++) {
			exp_vec[i] = exp(vec[i]);
		}
	};

	template <class T>
	static void linear_regression(const vector<T>& vals, int N, float& a, float& b) {
		float sum_x=0;
		float sum_y=0;
		float sum_x_squared=0;
		float sum_xy=0;
		for (int i=0; i<N; i++) {
			sum_x += i;
			sum_x_squared += i*i;
			sum_y += vals[i];
			sum_xy += i * vals[i];
		}
		a = (N*sum_xy - sum_x*sum_y) /
				(N*sum_x_squared - sum_x*sum_x);
		b = (sum_y - a*sum_x)/N;
	};

	template <class T>
	static void smooth_vector(const vector<T>& raw, vector<float>& smoothed, int smooth) {
		if (smooth <= 0) {
			copy(raw.begin(), raw.end(), smoothed.begin());
			return;
		}
		unsigned int bin;
		float sum=0;
		deque<float> recent(raw.begin(), raw.begin()+smooth*2);
		for (bin=0;bin<recent.size(); bin++) {
			sum+=recent[bin];
		}
		for (bin=smooth; bin < raw.size()-smooth; bin++) {
			sum += raw[bin+smooth];
			recent.push_back(raw[bin+smooth]);
			smoothed[bin] = sum/(2*smooth+1);
			sum -= recent[0];
			recent.pop_front();
		}

		float a_start_margin;
		float b_start_margin;
		linear_regression(raw, 2*smooth+1, a_start_margin, b_start_margin);

		vector<float> end_margin(2*smooth+1);
		std::reverse_copy(raw.end()-2*smooth-2, raw.end(), end_margin.begin());
		float a_end_margin;
		float b_end_margin;
		linear_regression(end_margin, 2*smooth+1, a_end_margin, b_end_margin);

		int s_bin=0;
		int e_bin = raw.size()-1;
		float s = b_start_margin;
		float e = b_end_margin;
		bool start_zero=1;
		bool end_zero=1;
		while (s_bin<smooth) {
			if (start_zero && raw[s_bin]>0) start_zero=0;
			if (end_zero && raw[e_bin]>0) end_zero=0;
			smoothed[s_bin] = (!start_zero) && s>0 ? s : 0;
			smoothed[e_bin] = (!end_zero) && e>0 ? e : 0;
			s_bin++;
			e_bin--;
			s += a_start_margin;
			e += a_end_margin;
		 } // smoothing the next window
		float weight=0;
		while (s_bin < 2*smooth) {
			smoothed[s_bin] = (1-weight)*s+ (weight*smoothed[s_bin]);
			smoothed[e_bin] = (1-weight)*e+ (weight*smoothed[e_bin]);
			weight+= 1.0/(float)smooth;
			s_bin++;
			e_bin--;
			s += a_start_margin;
			e += a_end_margin;
		}
	};
};

#endif /* VECTORUTILS_H_ */
