/*
 * Options.cpp
 *
 *  Created on: Nov 30, 2016
 *      Author: nettam
 */

#include "Options.h"
#include "macro.h"
#include <cstdlib>
using namespace std;

Options::~Options() {
	// TODO Auto-generated destructor stub
}

void Options::load_defaults(const char *vals[], int size)
{
	for(int i = 0; i < size; i+=2) {
		scalars["::"+string(vals[i])] = string(vals[i+1]);
	}
}

void Options::parse_argv(int &argc, char *argv[])
{
	int next_i = 1;
	for(int i = 1; i < argc; i++) {
		char *param = argv[i];

		if(param[0] == '@') {
			ifstream in(param + 1);
			ASSERT(in, "Cannot open paramter file " << param[0]);
			read(in);
		} else if(param[0] == '-') {
			string pstr = string(param + 1);
			string::size_type equal = pstr.find('=');
			if(equal == string::npos || equal == 0) {
				scalars["::"+pstr] = 1;
			} else {
				string cur_name = pstr.substr(0, equal);
				string cur_val = pstr.substr(equal + 1, pstr.length() - equal-1);
				scalars["::"+cur_name] = cur_val;
			}
		} else {
			argv[next_i] = argv[i];
			next_i++;
		}
	}
	argc = next_i;
}

void Options::read(istream &in, bool logit)
{
	m_logit = logit;

	//read lines

	string line;
	while(in) {
		getline(in, line, '\n');
		if(line.length() == 0) {
			continue;
		}

		if(m_logit)
			cerr << "new line is %%" << line << "%%" << endl;

		if(is_comment_line(line)) {
			continue;
		}
		if(is_scope_line(line)) {
			continue;
		}
		if(get_option_name(line)) {
			if(is_array_line(line, in)) {
				continue;
			}
			if(m_logit) {
				cerr << "not array line\n";
			}
			if(is_option_line(line)) {
				continue;
			}
			if(m_logit) {
				cerr << "not option line\n";
			}
		}
		ASSERT(false, "problem parsing line " << line);
	}
}

bool Options::is_comment_line(string &line)
{
	return(line[0] == '#');
}

bool Options::is_scope_line(string &line)
{
	if(line[0] != '[') {
		return(false);
	}

	cur_scope = line.substr(1, line.find(']') - 1);
	if(m_logit)
		cerr << "new scope " << cur_scope << endl;
	return(true);
}

bool Options::get_option_name(string &line)
{

	string::size_type equal = line.find('=');
	if(equal == string::npos || equal == 0) {
		return(false);
	}
	cur_name = line.substr(0, equal);

	string::size_type i = equal - 1;
	while(cur_name[i] == ' ') {
		cur_name.erase(i, i+1);
		i--;
	}

	string::size_type space = cur_name.find(' ');
	if(space != string::npos) {
		cerr << "bad option name (with space) -" << cur_name << "- pos is " << space << endl;
		return(false);
	}
	line.erase(0, equal + 1);
	while(line[0] == ' ') {
		line.erase(0,1);
	}

	if(m_logit) {
		cerr << "new option name is " << cur_name << endl;
		cerr << "line is left with " << line << endl;
	}

	return(true);
}

bool Options::is_array_line(string &line, istream &in)
{

	if(line[0] != '{') {
		return(false);
	}

	line.erase(0,1);
	while(line[0] == ' ') {
		line.erase(0,1);
	}

	string::size_type sp = line.find(' ');
	string atype = line.substr(0, sp);

	for(unsigned int i = 0; i <= sp; i++) {
		line.erase(0,1);
	}
	while(line[0] == ' ') {
		line.erase(0,1);
	}

	string asize;
	if(line.find(' ') != string::npos) {
		asize = line.substr(0, line.find(' '));
	} else {
		asize = line;
	}

	if(m_logit)
		cerr << "array " << cur_name << " type " << atype
						<< " size " << asize << endl;

	unsigned int size = atoi(asize.c_str());
	string arr_line;

	if(atype == "i") {
		vector<int> *ar = new vector<int>(size);
		for(unsigned int i = 0; i < size; i++) {
			in >> (*ar)[i];
			//test if it is ok
		}
		int_arrays[cur_scope+"::"+cur_name] = ar;
	} else if(atype == "f") {
		vector<float> *ar = new vector<float>(size);
		for(unsigned int i = 0; i < size; i++) {
			in >> (*ar)[i];
			//test if ok
		}
		float_arrays[cur_scope+"::"+cur_name] = ar;
	} else if(atype == "s") {
		vector<string> *ar = new vector<string>(size);
		for(unsigned int i = 0; i < size; i++) {
			getline(in, arr_line, '\n');
			(*ar)[i] = arr_line;
		}
		string_arrays[cur_scope+"::"+cur_name] = ar;
	} else {
		ASSERT(false, "Bad array type " << atype << " for name " << cur_name);
	}
	getline(in, line, '\n'); //new line of last item
	if(line[0] != '}') {
		if(size) {
			//suppose to be }
			getline(in, line, '\n');
		}
		if(line[0] != '}') {
			ASSERT(false, "After array items, line is " << line);
		}
	}
	if(m_logit)
		cerr << "out of array read\n";

	return(true);
}

bool Options::is_option_line(string &line)
{
	int i = line.length() - 1;

	while(line[i] == ' ') {
		line = line.erase(i,i+1);
		i--;
	}

	if(m_logit)
		cerr << "assigning to " << cur_name << " option " << line << endl;

	scalars[cur_scope+"::"+cur_name] = line;
	return(true); //no checks for now
}

int Options::get_int(const string &sc, const string &nm, bool must, int def)
{
	if(!scalars.count(sc+"::"+nm)) {
		ASSERT(!must, "int option does not exist: " << sc << "::" << nm);
		return(def);
	}
	string &sv = scalars[sc+"::"+nm];
	//type checking
	int val = atoi(sv.c_str());
	return(val);
}

float Options::get_float(const string &sc, const string &nm, bool must, float def)
{
	if(!scalars.count(sc+"::"+nm)) {
		ASSERT(!must, "float option does not exist: " << sc << "::" << nm);
		return(def);
	}
	string &sv = scalars[sc+"::"+nm];
	//type checking
	float val = atof(sv.c_str());
	return(val);
}

const string &Options::get_str(const string &sc, const string &nm, bool must, const string *def)
{
	static string empty;
	if(!scalars.count(sc+"::"+nm)) {
		ASSERT(!must, "str option does not exist: " << sc << "::" << nm);
		return(def!=0 ? *def : empty);
	}
	return(scalars[sc+"::"+nm]);
}

const vector<int> &Options::get_ints(const string &sc, const string &nm, bool must)
{
	if(!int_arrays.count(sc+"::"+nm)) {
		ASSERT(!must, "int array option does not exist: " << sc << "::" << nm);
		return(m_empty_int_vec);
	}
	return(*(int_arrays[sc+"::"+nm]));
}

const vector<float> &Options::get_floats(const string &sc, const string &nm, bool must)
{
	if(!float_arrays.count(sc+"::"+nm)) {
		ASSERT(!must, "float array option does not exist: " << sc << "::" << nm);
		return(m_empty_float_vec);
	}
	return(*(float_arrays[sc+"::"+nm]));

}

const vector<string> &Options::get_strs(const string &sc, const string &nm, bool must)
{
	if(!string_arrays.count(sc+"::"+nm)) {
		ASSERT(!must, "string array option does not exist: " << sc << "::" << nm);
		return(m_empty_str_vec);
	} else {
		return(*(string_arrays[sc+"::"+nm]));
	}
}
