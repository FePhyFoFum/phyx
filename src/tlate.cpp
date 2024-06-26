#include <string>
#include <map>
#include <algorithm>
#include <iostream>

#include "tlate.h"
#include "utils.h"
#include "sequence.h"
#include "seq_reader.h"


/*******************************************************/
// hard-coded translation tables //

// uses standard code (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1)
static std::map<std::string, std::string> standard_ = {
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
    {"GTC", "V"},
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
static std::map<std::string, std::string> vert_mtdna_ = {
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
    {"GTC", "V"},
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


// https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG3
static std::map<std::string, std::string> yeast_mtdna_ = {
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
    {"CTT", "T"},
    {"CTC", "T"},
    {"CTA", "T"},
    {"CTG", "T"},
    {"CCT", "P"},
    {"CCC", "P"},
    {"CCA", "P"},
    {"CCG", "P"},
    {"CAT", "H"},
    {"CAC", "H"},
    {"CAA", "Q"},
    {"CAG", "Q"},
    {"CGT", "R"},
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
    {"AGA", "R"},
    {"AGG", "R"},
    {"GTT", "V"},
    {"GTC", "V"},
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


// https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?chapter=cgencodes#SG5
static std::map<std::string, std::string> invert_mtdna_ = {
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
    {"AGA", "S"},
    {"AGG", "S"},
    {"GTT", "V"},
    {"GTC", "V"},
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


std::string TLATE::translate (std::string& dna) {
    std::string codon;
    std::string AminoAcid;
    std::string residue;
    //std::string AA; // not used
    
    dna = string_to_upper(dna);
    
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


TLATE::TLATE (const std::string& table) {
    // where the translation table is set
    if (table == "std") {
        table_ = standard_;
    } else if (table == "vmt") {
        table_ = vert_mtdna_;
    } else if (table == "ymt") {
        table_ = yeast_mtdna_;
    } else if (table == "ivmt") {
        table_ = invert_mtdna_;
    } else {
        std::cerr << "Error: table argument '" << table << "' not recognized. Exiting." << std::endl;
        exit(0);
    }
}

