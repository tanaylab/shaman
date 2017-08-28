/*
 * Parser.h
 *
 *  Created on: Nov 30, 2016
 *      Author: nettam
 */

#ifndef PARSER_H_
#define PARSER_H_
#include <vector>
#include <iostream>
using namespace std;

class Parser {
public:
	Parser();
	virtual ~Parser();
	int split_line(istream &in, vector<string> &fields,
			char delim='\t', int estimated_num_fiedls=1);
	int read_int_table(istream &in, unsigned int width, vector<vector<int> > &data);
};

#endif /* PARSER_H_ */
