/*
 * tlate.cpp
 *
 *  Created on: Jun 16, 2015
 *      Author: joe
 */

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <map>
#include <cstring>
#include <getopt.h>

using namespace std;

#include "tlate.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"

string Codon_to_AA(string& sequences) {

    map <string, string> AA_table;
    string AA = "";
    map <string, string>::iterator it;
    AA_table["TTT"] = "F";
    AA_table["TTC"] = "F";
    AA_table["TTA"] = "L";
    AA_table["TTG"] = "L";
    AA_table["TCT"] = "S";
    AA_table["TCC"] = "S";
    AA_table["TCA"] = "S";
    AA_table["TCG"] = "S";
    AA_table["TAT"] = "Y";
    AA_table["TAC"] = "Y";
    AA_table["TAA"] = "*";
    AA_table["TAG"] = "*";
    AA_table["TGT"] = "C";
    AA_table["TGC"] = "C";
    AA_table["TGA"] = "*";
    AA_table["TGG"] = "W";
    AA_table["CTT"] = "L";
    AA_table["CTC"] = "L";
    AA_table["CTA"] = "L";
    AA_table["CTG"] = "L";
    AA_table["CCT"] = "P";
    AA_table["CCC"] = "P";
    AA_table["CCA"] = "P";
    AA_table["CCG"] = "P";
    AA_table["CAT"] = "H";
    AA_table["CAC"] = "H";
    AA_table["CAA"] = "Q";
    AA_table["CAG"] = "Q";
    AA_table["CGT"] = "R";
    AA_table["CGC"] = "R";
    AA_table["CGA"] = "R";
    AA_table["CGG"] = "R";
    AA_table["ATT"] = "I";
    AA_table["ATC"] = "I";
    AA_table["ATA"] = "I";
    AA_table["ATG"] = "M";
    AA_table["ACT"] = "T";
    AA_table["ACC"] = "T";
    AA_table["ACA"] = "T";
    AA_table["ACG"] = "T";
    AA_table["AAT"] = "N";
    AA_table["AAC"] = "N";
    AA_table["AAA"] = "K";
    AA_table["AAG"] = "K";
    AA_table["AGT"] = "S";
    AA_table["AGC"] = "S";
    AA_table["AGA"] = "R";
    AA_table["AGG"] = "R";
    AA_table["GTT"] = "V";
    AA_table["GTA"] = "V";
    AA_table["GTG"] = "V";
    AA_table["GCT"] = "A";
    AA_table["GCC"] = "A";
    AA_table["GCA"] = "A";
    AA_table["GCG"] = "A";
    AA_table["GAT"] = "D";
    AA_table["GAC"] = "D";
    AA_table["GAA"] = "E";
    AA_table["GAG"] = "E";
    AA_table["GGT"] = "G";
    AA_table["GGC"] = "G";
    AA_table["GGA"] = "G";
    AA_table["GGG"] = "G";
    
    if (AA_table.find(sequences) == AA_table.end()) {
        AA = "X";
    } else {
        AA = AA_table[sequences];
    }
    return AA;
}

string TLATE::Translate (string& dna) {

    string codon = "";
    string AminoAcid = "";
    string AAcid;
    for (unsigned int i = 0; i < dna.size(); i = i + 3) {
        codon = dna[i];
        codon += dna[i+1];
        codon += dna[i+2];
        AAcid = Codon_to_AA(codon);
        AminoAcid += AAcid;
        codon = "";
    }
    return AminoAcid;
}

TLATE::TLATE() {
    // TODO Auto-generated constructor stub

}

