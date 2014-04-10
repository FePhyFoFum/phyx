
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


void populate_map_codon_indices(map<string,vector<int> > * codon_position){
    (*codon_position)["TTT"] = {0};
    (*codon_position)["TTC"] = {1};
    (*codon_position)["TTA"] = {2};
    (*codon_position)["TTG"] = {3};
    (*codon_position)["TCT"] = {4};
    (*codon_position)["TCC"] = {5};
    (*codon_position)["TCA"] = {6};
    (*codon_position)["TCG"] = {7};
    (*codon_position)["TAT"] = {8};
    (*codon_position)["TAC"] = {9};
    (*codon_position)["TGT"] = {10};
    (*codon_position)["TGC"] = {11};
    (*codon_position)["TGG"] = {12};
    (*codon_position)["CTT"] = {13};
    (*codon_position)["CTC"] = {14};
    (*codon_position)["CTA"] = {15};
    (*codon_position)["CTG"] = {16};
    (*codon_position)["CCT"] = {17};
    (*codon_position)["CCC"] = {18};
    (*codon_position)["CCA"] = {19};
    (*codon_position)["CCG"] = {20};
    (*codon_position)["CAT"] = {21};
    (*codon_position)["CAC"] = {22};
    (*codon_position)["CAA"] = {23};
    (*codon_position)["CAG"] = {24};
    (*codon_position)["CGT"] = {25};
    (*codon_position)["CGC"] = {26};
    (*codon_position)["CGA"] = {27};
    (*codon_position)["CGG"] = {28};
    (*codon_position)["ATT"] = {29};
    (*codon_position)["ATC"] = {30};
    (*codon_position)["ATA"] = {31};
    (*codon_position)["ATG"] = {32}; 
    (*codon_position)["ACT"] = {33};
    (*codon_position)["ACC"] = {34};
    (*codon_position)["ACA"] = {35};
    (*codon_position)["ACG"] = {36};
    (*codon_position)["AAT"] = {37};
    (*codon_position)["AAC"] = {38};
    (*codon_position)["AAA"] = {39};
    (*codon_position)["AAG"] = {40};
    (*codon_position)["AGT"] = {41};
    (*codon_position)["AGC"] = {42};
    (*codon_position)["AGA"] = {43};
    (*codon_position)["AGG"] = {44};
    (*codon_position)["GTT"] = {45};
    (*codon_position)["GTC"] = {46};
    (*codon_position)["GTA"] = {47};
    (*codon_position)["GTG"] = {48};
    (*codon_position)["GCT"] = {49};
    (*codon_position)["GCC"] = {50};
    (*codon_position)["GCA"] = {51};
    (*codon_position)["GCG"] = {52};
    (*codon_position)["GAT"] = {53};
    (*codon_position)["GAC"] = {54};
    (*codon_position)["GAA"] = {55};
    (*codon_position)["GAG"] = {56};
    (*codon_position)["GGT"] = {57};
    (*codon_position)["GGC"] = {58};
    (*codon_position)["GGA"] = {59};
    (*codon_position)["GGG"] = {60};
}
