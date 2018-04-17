/*
 * aatocdn.cpp
 *
 *  Created on: Jun 15, 2015
 *      Author: joe
 */

#include <iostream>
#include <string>
#include <map>

using namespace std;

#include "aa2cdn.h"

map <string, string> AAtoCDN::convert_to_codons(map <string, string>& aa_sequences,
    map<string, string>& nuc_sequences, bool& rm_last) {
    
    int len;
    string temp = "";
    for (iter_ = aa_sequences.begin(); iter_ != aa_sequences.end(); iter_++) {
        if (nuc_sequences.find(iter_ -> first) == nuc_sequences.end()) {
            cout << "Only in the AA File: " << iter_ -> first << endl;
        } else {
            amino_acid_sequence_ = iter_ -> second;
            nucleotide_sequence_ = nuc_sequences[iter_ -> first];
            if(rm_last == true){
                len = amino_acid_sequence_.size() - 1;
            }else{
                len = amino_acid_sequence_.size();
            }
            for (int i=0; i < len; i++) {
                if (amino_acid_sequence_[i] == '-') {
                    temp += "---";
                } else {
                    temp += nucleotide_sequence_[0];
                    temp += nucleotide_sequence_[1];
                    temp += nucleotide_sequence_[2];
                    nucleotide_sequence_.erase(0, 3);
                }
            }
            codon_sequences_[iter_ -> first] = temp;
            temp = "";
        }
    }
    return codon_sequences_;
}

AAtoCDN::AAtoCDN() {
    // TODO Auto-generated constructor stub

}
