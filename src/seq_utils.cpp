
#include <string>
#include <vector>
#include <map>
#include <assert.h>
#include <iostream>

using namespace std;

#include "seq_utils.h"
#include "sequence.h"

/**
 * TODO: need to make a alphabet guesser
 */

/**
 * procedure to get the character of a nucleotide
 * from the set of positions
 */
char get_dna_from_pos(set<int> ins){
    if(ins.count(0) == 1){
        if(ins.count(1) == 1){
            if(ins.count(2) == 1){
                if(ins.count(3) == 1){
                    return 'N';
                }
                return 'V';
            }
            if(ins.count(3) == 1){
                return 'H';
            }
            return 'M';
        }
        if(ins.count(2) == 1){
            return 'R';
        }
        if(ins.count(3) == 1){
            return 'W';
        }
        return 'A';
    }
    if(ins.count(1)==1){
        if(ins.count(2) == 1){
            if(ins.count(3) == 1){
                return 'B';
            }
            return 'S';
        }
        if(ins.count(3) == 1){
            return 'Y';
        }
        return 'C';
    }
    if(ins.count(2)==1){
        if(ins.count(3) == 1){
            return 'K';
        }
        return 'G';
    }
    if(ins.count(3)==1){
        return 'T';
    }
}

set<int> get_dna_pos(char inc){
    set<int> ret;
    inc = toupper(inc);
    if(inc == 'A'){
        ret.insert(0);
    }else if(inc == 'C'){
        ret.insert(1);
    }else if(inc == 'G'){
        ret.insert(2);
    }else if(inc == 'T'){
        ret.insert(3);
    }else if(inc == '-' || inc == 'N'){
        ret.insert(0);ret.insert(1);ret.insert(2);ret.insert(3);
    }else if(inc == 'Y'){
        ret.insert(1);ret.insert(3);
    }else if(inc == 'R'){
        ret.insert(0);ret.insert(2);
    }else if(inc == 'W'){
        ret.insert(0);ret.insert(3);
    }else if(inc == 'M'){
        ret.insert(0);ret.insert(1);
    }else if(inc == 'B'){
        ret.insert(1);ret.insert(2);ret.insert(3);
    }else if(inc == 'V'){
        ret.insert(0);ret.insert(1);ret.insert(2);
    }else if(inc == 'S'){
        ret.insert(1);ret.insert(2);
    }else if(inc == 'K'){
        ret.insert(2);ret.insert(3);
    }else if(inc == 'H'){
        ret.insert(0);ret.insert(1);ret.insert(3);
    }
    return ret;
}

/**
 * 
 * int alpha: the alphabet with 0=dna, 1=aa
 */
string consensus_seq(vector<Sequence> & seqs, int alpha){
    int seqlength = seqs[0].get_sequence().length();
    for(int i=0;i<seqs.size();i++){
        assert(seqs[i].get_sequence().length() == seqlength);
    }
    string retstring;
    for(int i=0;i<seqlength;i++){
        set<int> fullset;
        for(int j=0;j<seqs.size();j++){
            set<int> tset = get_dna_pos(seqs[j].get_sequence()[i]);
            fullset.insert(tset.begin(),tset.end());
        }
        retstring += get_dna_from_pos(fullset);
    }
    return retstring;
}

/**
 * Returns a map of DNA
 *
 */
char single_dna_complement(char inc){
    inc = toupper(inc);
    if(inc=='A'){
        return 'T';
    }else if(inc=='T'){
        return 'A';
    }else if(inc=='U'){
        return 'A';
    }else if(inc=='G'){
        return 'C';
    }else if(inc=='C'){
        return 'G';
    }else if(inc=='Y'){
        return 'R';
    }else if(inc=='R'){
        return 'Y';
    }else if(inc=='S'){
        return 'S';
    }else if(inc=='W'){
        return 'W';
    }else if(inc=='K'){
        return 'M';
    }else if(inc=='M'){
        return 'K';
    }else if(inc=='B'){
        return 'V';
    }else if(inc=='D'){
        return 'H';
    }else if(inc=='H'){
        return 'D';
    }else if(inc=='V'){
        return 'B';
    }else{
        return 'N';
    }
}

void write_phylip_alignment(vector<Sequence> & seqs, ostream * ostr){
    int seqlength = seqs[0].get_sequence().length();
    for(int i=0;i<seqs.size();i++){
        assert(seqs[i].get_sequence().length() == seqlength);
    }
    (*ostr) << seqs.size() << " " << seqlength << endl;
    for(int i=0;i<seqs.size();i++){
        (*ostr) << seqs[i].get_id() << "\t" << seqs[i].get_sequence() << endl;
    }
}

/**
 * this is not for concatenation. only single gene regions
 * another one needs to be written for concatenation
 */
void write_nexus_alignment(vector<Sequence> & seqs, ostream * ostr){
    int seqlength = seqs[0].get_sequence().length();
    for(int i=0;i<seqs.size();i++){
        assert(seqs[i].get_sequence().length() == seqlength);
    }
	(*ostr) << "#NEXUS" << endl;
	(*ostr) << "BEGIN DATA;\n\tDIMENSIONS NTAX=";
    (*ostr) << seqs.size() << " NCHAR=" << seqlength << ";" << endl;
	(*ostr) << "\tFORMAT DATATYPE="<< seqs[0].get_alpha_name() << " INTERLEAVE=NO GAP=-;"<<endl;
	(*ostr) << "\tMATRIX\n" << endl;
    for(int i=0;i<seqs.size();i++){
        (*ostr) << seqs[i].get_id() << "\t" << seqs[i].get_sequence() << endl;
    }

	(*ostr) << ";\nend;\n" << endl;
}

