
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
	Sequence cur;
	string curseq;
	while (getline(infile, tline)){
		TrimSpaces(tline);
		if(tline.size() < 1){
			continue;
		}
		if(tline.substr(0,1) == ">"){
			string id_ = tline.substr(1,tline.size()-1);
			if(first == true){
				first = false;
				cur = Sequence();
				cur.set_id(id_);
				curseq = "";
			}else{
				cur.set_sequence(curseq);
				seqs.push_back(cur);
				cur = Sequence();
				cur.set_id(id_);
				curseq = "";
			}
		}else{
			curseq += tline;
		}
	}
	cur.set_sequence(curseq);
	seqs.push_back(cur);
	infile.close();
	return true;
}
