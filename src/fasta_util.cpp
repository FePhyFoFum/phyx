
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

#include "sequence.h"
#include "fasta_util.h"
#include "utils.h"

//return false if not a fasta
bool read_fasta_file(string filen,vector<Sequence>& seqs) {
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

/*
 * this is just bare bones, write a vector of sequences to a file
 */
bool write_fasta_file(string filename, vector<Sequence> & seqs){
    ofstream outfile;
    outfile.open(filename.c_str(),ios::out);
    for (unsigned int i = 0; i < seqs.size(); i++){
	outfile << ">";
	outfile << seqs[i].get_id();
	outfile << "\n";
	outfile << seqs[i].get_sequence();
	outfile << "\n";
    }
    outfile.close();
}
