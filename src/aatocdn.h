/*
 * aatocdn.h
 *
 *  Created on: Jun 15, 2015
 *      Author: joe
 */

#ifndef AATOCDN_H_
#define AATOCDN_H_

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

class aa_to_cdn {
private:

	map<string, string> CodonAln;
	map<string, string> CodonSequences;
	map<string, string>::iterator iter;
	string AminoAcidSequence;
	string NucleotideSequence;
	string temp;
	string::iterator it;
	int j;

public:
	aa_to_cdn();
	virtual ~aa_to_cdn();
	map<string, string> FastaToOneLine(string&);
	map<string, string> ChangeToCodon(map<string, string>&, map<string,string>&);

};

#endif /* AATOCDN_H_ */
