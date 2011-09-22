
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

#include "sequence.h"
#include "phylip_reader.h"
#include "utils.h"

PhylipReader::PhylipReader(){}

//return false if not a phylip
bool PhylipReader::readFile(string filen,vector<Sequence>& seqs){
	string tline;
	ifstream infile(filen.c_str());
	bool first = true;
	while (getline(infile, tline)){
		vector<string> searchtokens;
		Tokenize(tline, searchtokens, "	");
		for(unsigned int j=0;j<searchtokens.size();j++){
			TrimSpaces(searchtokens[j]);
		}
		if (first == true){ //reading past the first line
			first = false;
			try{
				//int nseqs = atoi(searchtokens[0].c_str());
				//int nsites = atoi(searchtokens[1].c_str());
			}catch( char * str ){
				return false; //not a phylip
			}
			continue;
		}
		if (searchtokens.size() < 2)
			continue;
		Sequence a = Sequence(searchtokens[0],searchtokens[1],true);
		seqs.push_back(a);
	}
	infile.close();
	return true;
}
