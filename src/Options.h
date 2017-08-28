/*
 * Options.h
 *
 *  Created on: Nov 30, 2016
 *      Author: nettam
 */

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

class Options {
public:
	Options() {}
	Options(istream &in) { read(in); }
	virtual ~Options();

	void load_defaults(const char *vals[], int size);
	void parse_argv(int &argc, char *argv[]);
	void read(istream &in, bool logit = false);
	//Access methods

	int get_int(const string &sc, const string &nm, bool must = true, int def = 0);
	int get_g_int(const string &nm, bool must = true, int def = 0) {
		return(get_int("", nm, must, def));
	}
	float get_float(const string &sc, const string &nm, bool must = true, float def = 0);
	float get_g_float(const string &nm, bool must = true, int def = 0) {
		return(get_float("", nm, must, def));
	}
	const string &get_str(const string &sc, const string &nm, bool must = true, const string *def = 0);
	const string &get_g_str(const string &nm, bool must = true, const string *def = 0) {
		return(get_str("", nm, must, def));
	}

	const vector<int> &get_ints(const string &sc, const string &nm,
							bool must = true);
	const vector<float> &get_floats(const string &sc, const string &nm,
							bool must = true);
	const vector<string> &get_strs(const string &sc, const string &nm,
							bool must = true);

	bool have_scalar(const string &sc, const string &nm) {
		return(scalars.count(sc+"::"+nm));
	}
	bool have_g_scalar(const string &nm) {
		return(scalars.count("::"+nm));
	}

	void set_scalar(const string &sc, const string &nm, const string &val) {
		scalars[sc+"::"+nm] = val;
	}
	void set_g_scalar(const string &nm, const string &val) {
		scalars["::"+nm] = val;
	}
private :

//Parsing methods

	bool is_comment_line(string &line);
	bool is_scope_line(string &line);
	bool get_option_name(string &line);
	bool is_option_line(string &line);
	bool is_array_line(string &line, istream &in);


private:

	map<string, string> scalars;
	map<string, vector<int> *> int_arrays;
	map<string, vector<float> *> float_arrays;
	map<string, vector<string> *> string_arrays;

	vector<int> m_empty_int_vec;
	vector<float> m_empty_float_vec;
	vector<string> m_empty_str_vec;

	//auxilaries

	string cur_scope;
	string cur_name;
public:
	bool m_logit;

};

#endif /* OPTIONS_H_ */
