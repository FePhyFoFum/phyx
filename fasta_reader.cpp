
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

#include "sequence.h"
#include "fasta_reader.h"
#include "utils.h"

FastaReader::FastaReader(){}

//return false if not a fasta
bool FastaReader::readFile(string filen,vector<Sequence>& seqs){
	string tline;
	ifstream infile(filen.c_str());
	bool first = true;
	Sequence * cur;
	while (getline(infile, tline)){
		vector<string> searchtokens;
		Tokenize(tline, searchtokens, " 	");
		for(unsigned int j=0;j<searchtokens.size();j++){
			TrimSpaces(searchtokens[j]);
		}
	}
	infile.close();
	return true;
}
