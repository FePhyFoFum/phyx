
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

#include "sequence.h"
#include "seq_reader.h"
#include "utils.h"

/*
 * tests the filetype by checking the first string and guessing based on
 * # (nexus), num (phylip), > (fasta), @ (fastq)
 * returns in the order above, 0, 1, 2, 3, 666 -- no filetype recognized
 */
int test_seq_filetype(string filen){
    string tline;
    ifstream infile(filen.c_str());
    int ret = 666; // if you get 666, there is no filetype set
    while(getline(infile,tline)){
	if (tline.size() < 1){
	    continue;
	}
	if (tline[0] == '#'){
	    ret = 0;
	    break;
	}
	vector<string> tokens;
	string del(" \t");
	tokenize(tline,tokens,del);
	if (tokens.size() > 1){
	    trim_spaces(tokens[0]);
	    if (is_number(tokens[0])){
		ret = 1;
		break;
	    }
	}
	if (tline[0] == '>'){
	    ret = 2;
	    break;
	}
	if (tline[0] == '@'){
	    ret = 3;
	    break;
	}
	break;
    }
    infile.close();
    return ret;
}

//return false if not a fasta
bool read_fasta_file(string filen,vector<Sequence>& seqs) {
    string tline;
    ifstream infile(filen.c_str());
    bool first = true;
    Sequence cur;
    string curseq;
    while (getline(infile, tline)){
	trim_spaces(tline);
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


bool read_phylip_file(string filen,vector<Sequence>& seqs){
    string tline;
    ifstream infile(filen.c_str());
    bool first = true;
    while (getline(infile, tline)){
	vector<string> searchtokens;
	tokenize(tline, searchtokens, "	");
	for(unsigned int j=0;j<searchtokens.size();j++){
	    trim_spaces(searchtokens[j]);
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
