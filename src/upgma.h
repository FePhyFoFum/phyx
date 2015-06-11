/*
 * upgma.h
 *
 *  Created on: Jun 10, 2015
 *      Author: joe
 */

#ifndef UPGMA_H_
#define UPGMA_H_

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


class UPGMA {
private:

	map <string, string> sequences;
	map <string, int> NumbKey;
	map <string, string>::iterator iter;
	vector<string> names;
	string fasta;
	vector< vector<double> > Matrix;


public:
	UPGMA();
	map<string, string> FastaToOneLine(string& fasta);
	vector< vector<double> > BuildMatrix(map<string, string>& sequences);
	void TREEMAKE(vector<string>&, map <int, string>&, vector< vector<double> >&);
	//virtual ~UPGMA();
};

#endif /* UPGMA_H_ */
