
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

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
    while (getline(infile,tline)){
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

/*
 * tests the filetype by checking the first string and guessing based on
 * # (nexus), num (phylip), > (fasta), @ (fastq)
 * returns in the order above, 0, 1, 2, 3, 666 -- no filetype recognized
 */
int test_seq_filetype_stream(istream & stri,string & retstring){
    if (!getline(stri, retstring)){
        cout << "ERROR: end of file too soon" << endl;
    }
    int ret = 666; // if you get 666, there is no filetype set
    if (retstring[0] == '#'){
        ret = 0;
    }else if (retstring[0] == '>'){
        ret = 2;
    }else if (retstring[0] == '@'){
        ret = 3;
    }else{    
        vector<string> tokens;
        string del(" \t");
        tokenize(retstring,tokens,del);
        if (tokens.size() > 1){
            trim_spaces(tokens[0]);
            if (is_number(tokens[0])){
                ret = 1;
            }
        }
    }
    return ret;
}

/*
 * returns the next string in the getline if there is one
 * TODO: nexus interleaved is not going to work here
 */
bool read_next_seq_from_stream(istream & stri, int ftype, string & retstring, Sequence & seq){
    string tline;
    if (ftype == 0){ // nexus
        string tline;
        //are we at the beginning of the file?
        //TODO: add check for interleave and kick out to do a different reader
        //checks for beginning of char by MATRIX
        if (retstring.size() > 0 && retstring[0] == '#'){
            bool found = false;
            while (getline(stri,tline)){
                trim_spaces(tline);
                std::transform(tline.begin(), tline.end(),tline.begin(), ::toupper);    
                if (tline.compare("MATRIX")==0){
                    found = true;
                    break;
                }
            }
            if (found == false){
                cout << "badly formatted nexus file, missing 'MATRIX' in data/character block" << endl;
            }
            retstring = "";
        }
        getline(stri,tline);
        trim_spaces(tline);
        while (tline.size() == 0){
            if (getline(stri,tline)){
                trim_spaces(tline);
            }else{
                return false;
            }
        }
        vector<string> tokens;
        string del(" \t");
        tokenize(tline,tokens,del);
        if (tokens.size() > 1){
            for(int i=0;i<tokens.size();i++){
                trim_spaces(tokens[i]);
            }
            if (tokens[0].compare(";") == 0){
                return false;
            }else{
                seq.set_id(tokens[0]);
                seq.set_sequence(tokens[1]);
                return true;
            }
        }else{
            return false;
        }
    }else if (ftype == 1){ // phylip
        vector<string> tokens;
        string del(" \t");
        string tline;
        //check to see if we are at the beginning of the file
        if (retstring.size() > 0){
            tokenize(retstring,tokens,del);
            if (tokens.size() > 1){
                trim_spaces(tokens[0]);
                if (is_number(tokens[0])){
                    getline(stri,tline);
                }else{
                    tline = retstring;
                }
            }
            retstring = "";
        }
        if (tline.size() == 0){
            if (!getline(stri,tline)){
            return false;
            }
        }
        tokens.clear();
        tokenize(tline,tokens,del);
        for(int i=0;i<tokens.size();i++){
            trim_spaces(tokens[i]);
        }
        if (tokens[0].size() == 0){
            return false;
        }
        seq.set_id(tokens[0]);
        seq.set_sequence(tokens[1]);
        return true;
    }else if (ftype == 2){ // fasta
        bool first = true;
        bool going = true;
        string curseq = "";
        while (going){
            if (first == true && retstring.size() > 0){
                tline = retstring;
                retstring = "";
            }else{
                if (!getline(stri, tline)){
                    seq.set_sequence(curseq);
                    return false;
                }
            }
            if (tline.substr(0,1) == ">"){
                if (first == true){
                    string id_ = tline.substr(1,tline.size()-1);
                    first = false;
                    seq.set_id(id_);
                    curseq = "";
                }else{
                    seq.set_sequence(curseq);
                    retstring = tline;
                    return true;
                }
            }else{
                curseq.append(tline);
            }
        }
    }else if (ftype == 3){//fastq assumes a 33 offset for now
        string line1,line2,line3,line4;
        if (retstring.size() > 0){
            line1 = retstring;
            retstring = "";
        }else if (!getline(stri,line1)){
            return false;
        }
        if (!getline(stri,line2)){
            return false;
        }if (!getline(stri,line3)){
            return false;
        }if (!getline(stri,line4)){
            return false;
        }
        seq.set_id(line1.substr(1,line1.size()-1));
        seq.set_sequence(line2);
        seq.set_qualstr(line4,33);
        return true;
    }
    return false;
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
        if (tline.size() < 1){
            continue;
        }
        if (tline.substr(0,1) == ">"){
            string id_ = tline.substr(1,tline.size()-1);
            if (first == true){
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
        tokenize(tline, searchtokens, "    ");
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
        if (searchtokens.size() < 2) {
            continue;
        }
        Sequence a = Sequence(searchtokens[0],searchtokens[1],true);
        seqs.push_back(a);
    }
    infile.close();
    return true;
}

//TODO: INCOMPLETE
bool read_nexus_seqs_file(string filen,vector<Sequence>& seqs){
    string tline;
    ifstream infile(filen.c_str());
    bool first = true;
    while (getline(infile, tline)){
        vector<string> searchtokens;
        tokenize(tline, searchtokens, "    ");
        for (unsigned int j=0;j<searchtokens.size();j++){
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
        if (searchtokens.size() < 2) {
            continue;
        }
        Sequence a = Sequence(searchtokens[0],searchtokens[1],true);
        seqs.push_back(a);
    }
    infile.close();
    return true;
}
