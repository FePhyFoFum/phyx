/*
 * tlate.cpp
 *
 *  Created on: Jun 16, 2015
 *      Author: joe
 */

#include <map>
#include <algorithm>

using namespace std;

#include "tlate.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"

/*******************************************************/
// hard-coded translation tables //

// uses standard code (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1)
map <string, string> standard_ = {
    {"TTT", "F"},
    {"TTC", "F"},
    {"TTA", "L"},
    {"TTG", "L"},
    {"TCT", "S"},
    {"TCC", "S"},
    {"TCA", "S"},
    {"TCG", "S"},
    {"TAT", "Y"},
    {"TAC", "Y"},
    {"TAA", "*"},
    {"TAG", "*"},
    {"TGT", "C"},
    {"TGC", "C"},
    {"TGA", "*"},
    {"TGG", "W"},
    {"CTT", "L"},
    {"CTC", "L"},
    {"CTA", "L"},
    {"CTG", "L"},
    {"CCT", "P"},
    {"CCC", "P"},
    {"CCA", "P"},
    {"CCG", "P"},
    {"CAT", "H"},
    {"CAC", "H"},
    {"CAA", "Q"},
    {"CAG", "Q"},
    {"CGT", "R"},
    {"CGC", "R"},
    {"CGA", "R"},
    {"CGG", "R"},
    {"ATT", "I"},
    {"ATC", "I"},
    {"ATA", "I"},
    {"ATG", "M"},
    {"ACT", "T"},
    {"ACC", "T"},
    {"ACA", "T"},
    {"ACG", "T"},
    {"AAT", "N"},
    {"AAC", "N"},
    {"AAA", "K"},
    {"AAG", "K"},
    {"AGT", "S"},
    {"AGC", "S"},
    {"AGA", "R"},
    {"AGG", "R"},
    {"GTT", "V"},
    {"GTA", "V"},
    {"GTG", "V"},
    {"GCT", "A"},
    {"GCC", "A"},
    {"GCA", "A"},
    {"GCG", "A"},
    {"GAT", "D"},
    {"GAC", "D"},
    {"GAA", "E"},
    {"GAG", "E"},
    {"GGT", "G"},
    {"GGC", "G"},
    {"GGA", "G"},
    {"GGG", "G"}
};

// https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG2
map <string, string> vert_mtdna_ = {
    {"TTT", "F"},
    {"TTC", "F"},
    {"TTA", "L"},
    {"TTG", "L"},
    {"TCT", "S"},
    {"TCC", "S"},
    {"TCA", "S"},
    {"TCG", "S"},
    {"TAT", "Y"},
    {"TAC", "Y"},
    {"TAA", "*"},
    {"TAG", "*"},
    {"TGT", "C"},
    {"TGC", "C"},
    {"TGA", "W"},
    {"TGG", "W"},
    {"CTT", "L"},
    {"CTC", "L"},
    {"CTA", "L"},
    {"CTG", "L"},
    {"CCT", "P"},
    {"CCC", "P"},
    {"CCA", "P"},
    {"CCG", "P"},
    {"CAT", "H"},
    {"CAC", "H"},
    {"CAA", "Q"},
    {"CAG", "Q"},
    {"CGT", "R"},
    {"CGC", "R"},
    {"CGA", "R"},
    {"CGG", "R"},
    {"ATT", "I"},
    {"ATC", "I"},
    {"ATA", "M"},
    {"ATG", "M"},
    {"ACT", "T"},
    {"ACC", "T"},
    {"ACA", "T"},
    {"ACG", "T"},
    {"AAT", "N"},
    {"AAC", "N"},
    {"AAA", "K"},
    {"AAG", "K"},
    {"AGT", "S"},
    {"AGC", "S"},
    {"AGA", "*"},
    {"AGG", "*"},
    {"GTT", "V"},
    {"GTA", "V"},
    {"GTG", "V"},
    {"GCT", "A"},
    {"GCC", "A"},
    {"GCA", "A"},
    {"GCG", "A"},
    {"GAT", "D"},
    {"GAC", "D"},
    {"GAA", "E"},
    {"GAG", "E"},
    {"GGT", "G"},
    {"GGC", "G"},
    {"GGA", "G"},
    {"GGG", "G"}
};

/*******************************************************/

string TLATE::translate (string& dna) {

    string codon = "";
    string AminoAcid = "";
    string residue;
    string AA;
    
    std::transform(dna.begin(), dna.end(), dna.begin(), ::toupper);
    
    for (unsigned int i = 0; i < dna.size(); i = i + 3) {
        codon = dna.substr(i, 3);
        if (table_.find(codon) == table_.end()) {
            residue = "X";
        } else {
            residue = table_[codon];
        }
        AminoAcid += residue;
        codon = "";
    }
    return AminoAcid;
}

TLATE::TLATE (string const& table) {
    if (table != "std" && table != "vmt") {
        cout << "Table argument '" << table << "' not recognized. Exiting." << endl;
        exit(0);
    }
    
    // where the translation table is set
    if (table == "std") {
        table_ = standard_;
    } else if (table == "vmt") {
        table_ = vert_mtdna_;
    }
}

