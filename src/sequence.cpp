/*
 * This is an extremely stripped down sequence class that is just meant to be
 * transparent and lightweight. As functionality increases, so will the 
 * complexity of the class.
 *
 */

#include <string>

#include <iostream>


using namespace std;

#include "sequence.h"
#include "utils.h"
#include "seq_utils.h"

Sequence::Sequence():id(),seq(),aligned(),alphabet(DNA){}

Sequence::Sequence(string _id, string _seq, bool _aligned){
    id = _id;
    seq = _seq;
    aligned = _aligned;
    alphabet = DNA;
}
Sequence::Sequence(string _id, string _seq){
    id = _id;
    seq = _seq;
    aligned = false;
    alphabet = DNA;
}

seqAlpha Sequence::get_alpha(){
    return alphabet;
}

string Sequence::get_alpha_name(){
    if (alphabet == DNA){
        return "DNA";
    }
    if (alphabet == AA){
        return "AA";
    }
    if (alphabet == BINARY){
        return "BINARY";
    }
    if(alphabet == MULTI){
        return "MULTI";
    }
    return "";
}

void Sequence::set_alpha(seqAlpha s){
    alphabet = s;
}

bool Sequence::is_aligned(){
    return aligned;
}

string Sequence::get_sequence(){
    return seq;
}

string Sequence::get_id(){
    return id;
}

void Sequence::add_cont_char(double _num){
    cont_chars.push_back(_num);
}

double Sequence::get_cont_char(int _index){
    return cont_chars[_index];
}

int Sequence::get_num_cont_char(){
    return cont_chars.size();
}

void Sequence::clear_cont_char(){
    cont_chars.clear();
}

void Sequence::add_multistate_char(int _num){
    multistate_chars.push_back(_num);
}

int Sequence::get_multistate_char(int _index){
    return multistate_chars[_index];
}

int Sequence::get_num_multistate_char(){
    return multistate_chars.size();
}


void Sequence::set_sequence(string _seq){
    seq = _seq;
}

void Sequence::set_id(string _id){
    id = _id;
}

void Sequence::set_aligned(bool _aligned){
    aligned = _aligned;
}

string Sequence::reverse_complement(){
    string rcomp = seq;
    for (unsigned int i=0 ;i < rcomp.size(); i++){
        rcomp.replace(i,1,1,single_dna_complement(seq[seq.size()-i-1]));
    }
    return rcomp;
}

void Sequence::perm_reverse_complement(){
    string rcomp = seq;
    for (unsigned int i=0 ;i < rcomp.size(); i++){
        rcomp.replace(i,1,1,single_dna_complement(seq[seq.size()-i-1]));
    }
    seq = rcomp;
}

void Sequence::set_qualstr(string & stri,int offset){
    qualarr.clear();
    qualstr = stri;
    for (int i=0;i<stri.size();i++){
        qualarr.push_back(((int)stri[i])-offset);
    }
}

vector<double> Sequence::get_qualarr(){
    return qualarr;
}

double Sequence::get_qualarr_mean(){
    return calculate_vector_double_mean(qualarr);
}

string Sequence::get_fasta(){
    string retstr;
    retstr.append(">");retstr.append(id);retstr.append("\n");
    retstr.append(seq);retstr.append("\n");
    return retstr;
}

string Sequence::get_fastq(){
    string retstr;
    retstr.append("@");retstr.append(id);retstr.append("\n");
    retstr.append(seq);retstr.append("\n+\n");
    retstr.append(qualstr);retstr.append("\n");
    return retstr;
}

