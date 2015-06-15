/*
 * clsq.h
 *
 *  Created on: Jun 15, 2015
 *      Author: joe
 */

#ifndef CLSQ_H_
#define CLSQ_H_

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <map>
#include <iterator>
using namespace std;


class clsq {
private:

	string fasta, line, dna, name_hold;
	map<string, string> sequences;
	map<string, string>::iterator iter;
	map<string, string> Trimmed;

public:
	clsq();
	map<string, string> FastaToOneLine(string &fasta, double& MissingAllowed);
	virtual ~clsq();
};

#endif /* CLSQ_H_ */
