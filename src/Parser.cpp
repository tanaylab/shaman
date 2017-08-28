/*
 * Parser.cpp
 *
 *  Created on: Nov 30, 2016
 *      Author: nettam
 */
#include <stdlib.h>

#include "Parser.h"
#include "macro.h"

Parser::Parser() {
	// TODO Auto-generated constructor stub

}

Parser::~Parser() {
	// TODO Auto-generated destructor stub
}

int Parser::split_line(istream &in, vector<string> &fields, char delim, int estimated_num_fields)
{
	int num_lines = 0;
	fields.resize(estimated_num_fields);
	for (vector<string>::iterator istr = fields.begin(); istr != fields.end(); ++istr)
		istr->resize(0);
	vector<string>::iterator istr = fields.begin();
	while(in) {
		int c = in.get();
		if(c == '\r') {
			continue;
		}

		if (c == '\n')
			num_lines++;

		if (c == '\n' || !in.good()) {
			if (istr == fields.begin() && istr->empty()) {
				if (!in.good()) {
					fields.clear();
					break;
				}

				// eat up empty lines
				continue;
			}
			fields.resize(istr - fields.begin() + 1);
			break;
		}
		if (c == delim) {
			istr++;
			if (istr == fields.end()) {
				fields.push_back(string());
				istr = fields.begin() + fields.size() - 1;
			}
		} else
			istr->push_back(c);
	}
	return num_lines;
}

int Parser::read_int_table(istream &in, unsigned int width, vector<vector<int> > &data)
{
	vector<string> fields;
	int row = 0;
	while(in) {
		split_line(in, fields);
		if(fields.size() == 0) {
			return(-1);
		}
		ASSERT(fields.size() == width, "Bad table width (" << fields.size() << " instead " << width << ") when parsing int table");
		data.resize(row + 1, vector<int>(width));
		vector<int>::iterator dt = data[row].begin();
		vector<string>::const_iterator inp = fields.begin();
		while(inp != fields.end()) {
			char *fin;
			*dt = strtol((*inp).c_str(), &fin, 0);
			ASSERT((fin - (*inp).c_str()) == inp->size(), "Cannot parse int at row " << row << " col " << inp-fields.begin());
			dt++;
			inp++;
		}
		row++;
	}
	return(row);
}

