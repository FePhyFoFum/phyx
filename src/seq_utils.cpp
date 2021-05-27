#include <string>
#include <vector>
#include <set>
#include <map>
#include <cassert>
#include <iostream>
#include <algorithm>

#include "seq_utils.h"
#include "sequence.h"
#include "utils.h"


// IUPAC codes
std::string dnachars_with_ambiguous = "ACGTURYSWKMBDHVN";
std::string protchars = "ABCDEFGHIKLMNPQRSTVWXYZ";


/**
 * IUPAC ambiguity codes
 * procedure to get the character of a nucleotide
 * from the set of positions
 * An empty set contains only gaps. An N has a count for all nucleotides.
 * When a site has valid nucleotides and gaps, the gaps are ignored.
 */
char get_dna_from_pos (const std::set<int>& ins) {
    if (ins.count(0) == 1) {
        if (ins.count(1) == 1) {
            if (ins.count(2) == 1) {
                if (ins.count(3) == 1) {
                    return 'N';
                }
                return 'V';
            }
            if (ins.count(3) == 1) {
                return 'H';
            }
            return 'M';
        }
        if (ins.count(2) == 1) {
            if (ins.count(3) == 1) {
                return 'D';
            }
            return 'R';
        }
        if (ins.count(3) == 1) {
            return 'W';
        }
        return 'A';
    }
    if (ins.count(1) == 1) {
        if (ins.count(2) == 1) {
            if (ins.count(3) == 1) {
                return 'B';
            }
            return 'S';
        }
        if (ins.count(3) == 1) {
            return 'Y';
        }
        return 'C';
    }
    if (ins.count(2) == 1) {
        if (ins.count(3) == 1) {
            return 'K';
        }
        return 'G';
    }
    if (ins.count(3) == 1) {
        return 'T';
    }
    return('-');
}


std::set<int> get_dna_pos (char inc) {
    std::set<int> ret;
    inc = toupper(inc);
    if (inc == 'A') {
        ret.insert(0);
    } else if (inc == 'C') {
        ret.insert(1);
    } else if (inc == 'G') {
        ret.insert(2);
    } else if (inc == 'T') {
        ret.insert(3);
    //} else if (inc == '-' || inc == 'N') {
    } else if (inc == 'N') {
        ret.insert(0);
        ret.insert(1);
        ret.insert(2);
        ret.insert(3);
    } else if (inc == 'Y') {
        ret.insert(1);
        ret.insert(3);
    } else if (inc == 'R') {
        ret.insert(0);
        ret.insert(2);
    } else if (inc == 'W') {
        ret.insert(0);
        ret.insert(3);
    } else if (inc == 'M') {
        ret.insert(0);
        ret.insert(1);
    } else if (inc == 'B') {
        ret.insert(1);
        ret.insert(2);
        ret.insert(3);
    } else if (inc == 'V') {
        ret.insert(0);
        ret.insert(1);
        ret.insert(2);
    } else if (inc == 'S') {
        ret.insert(1);
        ret.insert(2);
    } else if (inc == 'K') {
        ret.insert(2);
        ret.insert(3);
    } else if (inc == 'H') {
        ret.insert(0);
        ret.insert(1);
        ret.insert(3);
    }
    return ret;
}


char get_prot_char (const std::set<char>& inc) {
    // if any is missing, consensus is missing
    if (inc.count('X') == 1 || inc.count('-') == 1) {
        return 'X';
    }
    if (inc.size() == 1) {
        return *inc.begin();
    }
    
    // there are a handful of ambiguity codes
    // B = Aspartic acid (D) or Asparagine (N)
    // Z = Glutamine (Q) or Glutamic acid (E)
    int B = static_cast<int>(inc.count('D') + inc.count('N') + inc.count('B'));
    int Z = static_cast<int>(inc.count('Q') + inc.count('E') + inc.count('Z'));
    
    if (B == static_cast<int>(inc.size())) {
        return 'B';
    }
    if (Z == static_cast<int>(inc.size())) {
        return 'Z';
    }
    
    return 'X';
}


/**
 * string alpha: either "DNA" or "AA"
*/
std::string consensus_seq (std::vector<Sequence>& seqs, std::string& alpha) {
    bool aligned = is_aligned(seqs);
    if (!aligned) {
        std::cerr << "Error: sequences are not aligned. Exiting." << std::endl;
        exit(0);
    }
    unsigned int seqlength = seqs[0].get_length();
    std::string retstring;
    if (alpha == "DNA") {
        for (unsigned int i = 0; i < seqlength; i++) {
            std::set<int> fullset;
            for (auto & seq : seqs) {
                std::set<int> tset = get_dna_pos(seq.get_sequence()[i]);
                fullset.insert(tset.begin(), tset.end());
            }
            retstring += get_dna_from_pos(fullset);
        }
    } else if (alpha == "AA") {
        for (unsigned int i = 0; i < seqlength; i++) {
            std::set<char> fullset;
            //bool ambig = false; // doesn't do anything
            for (auto & seq : seqs) {
                fullset.insert(seq.get_sequence()[i]);
                // break early if any ambiguous code is encountered
                if (seq.get_sequence()[i] == 'X' || seq.get_sequence()[i] == '-') {
                    //ambig = true;
                    break;
                }
            }
            retstring += get_prot_char(fullset);
        }
    } else {
        std::cerr << "Error: cannot make consensus of sequence type '" << alpha
                << "'. Exiting." << std::endl;
        exit(0);
    }
    return retstring;
}


/**
 * Returns a map of DNA
 *
 */
char single_dna_complement (char inc) {
    inc = toupper(inc);
    if (inc == 'A') {
        return 'T';
    } else if (inc == 'T') {
        return 'A';
    } else if (inc == 'U') {
        return 'A';
    } else if (inc == 'G') {
        return 'C';
    } else if (inc == 'C') {
        return 'G';
    } else if (inc == 'Y') {
        return 'R';
    } else if (inc == 'R') {
        return 'Y';
    } else if (inc == 'S') {
        return 'S';
    } else if (inc == 'W') {
        return 'W';
    } else if (inc == 'K') {
        return 'M';
    } else if (inc == 'M') {
        return 'K';
    } else if (inc == 'B') {
        return 'V';
    } else if (inc == 'D') {
        return 'H';
    } else if (inc == 'H') {
        return 'D';
    } else if (inc == 'V') {
        return 'B';
    } else {
        return 'N';
    }
}


void write_phylip_alignment (std::vector<Sequence>& seqs, const bool& uppercase, std::ostream * ostr) {
    bool aligned = is_aligned(seqs);
    if (!aligned) {
        std::cerr << "Error: sequences are not aligned. Exiting." << std::endl;
        exit(0);
    }
    int seqlength = seqs[0].get_length();
    // header: num_taxa, num_char
    (*ostr) << seqs.size() << " " << seqlength << std::endl;
    
    for (auto & seq : seqs) {
        if (uppercase) {
            std::string terp = seq.seq_to_upper();
            (*ostr) << seq.get_id() << "\t" << terp << std::endl;
        } else {
            (*ostr) << seq.get_id() << "\t" << seq.get_sequence() << std::endl;
        }
    }
}


/**
 * this is not for concatenation. only single gene regions
 * another one needs to be written for concatenation
 */
void write_nexus_alignment (std::vector<Sequence>& seqs, const bool& uppercase, std::ostream * ostr) {
    unsigned int seqlength = seqs[0].get_length();
    std::string datatype = seqs[0].get_alpha_name();
    std::string symbols; // not required for binary (default), protein, DNA
    //std::cout << std::endl << "datatype = " << datatype << std::endl << std::endl;
    
    if (datatype == "AA") { // "AA" is not a valid Nexus datatype
        datatype = "PROTEIN";
    }
    if (datatype == "BINARY") {
        datatype = "STANDARD";
        symbols = "01"; // not strictly necessary, as this is the default set. but doesn't hurt
    }
    if (datatype == "MULTI") {
        datatype = "STANDARD";
        //std::cout << "assembling symbols now" << std::endl;
        std::string combined;
        for (auto & seq : seqs) {
            if (uppercase) {
                combined += seq.seq_to_upper();
            } else {
                combined += seq.get_sequence();
            }
        }
        symbols = get_alphabet_from_sequence(combined);
        symbols.erase(std::remove(symbols.begin(), symbols.end(), '-'), symbols.end());
        symbols.erase(std::remove(symbols.begin(), symbols.end(), '?'), symbols.end());
    }
    
    bool aligned = is_aligned(seqs);
    if (!aligned) {
        std::cerr << "Error: sequences are not aligned. Exiting." << std::endl;
        exit(0);
    }
    
    (*ostr) << "#NEXUS" << std::endl;
    (*ostr) << "BEGIN DATA;\n\tDIMENSIONS NTAX=";
    (*ostr) << seqs.size() << " NCHAR=" << seqlength << ";" << std::endl;
    (*ostr) << "\tFORMAT DATATYPE=" << datatype;
    if (!symbols.empty()) {
        (*ostr) << " SYMBOLS=\"" << symbols << "\"";
    }
    // removing 'INTERLEAVE=NO` as it appears to not be NCL-compliant
    (*ostr) << " GAP=- MISSING=?;" << std::endl;
    (*ostr) << "\tMATRIX\n" << std::endl;
    for (auto & seq : seqs) {
        // MrBayes is not Nexus-compliant, so using a "safe" version
        if (uppercase) {
            std::string terp = seq.seq_to_upper();
            (*ostr) << get_valid_nexus_label(seq.get_id()) << "\t" << terp << std::endl;
        } else {
            (*ostr) << get_valid_nexus_label(seq.get_id()) << "\t" << seq.get_sequence() << std::endl;
        }
        //(*ostr) << seq.get_id() << "\t" << seq.get_sequence() << std::endl;
    }
    (*ostr) << ";\nend;\n" << std::endl;
}


/**
 * given a vector of seqs, this will make, for each seq a vector
 * that would be 000000100000 for each codon that is present
 * this would be for each site, so reuse your vectors!
 */
void create_vector_seq_codon_state_reconstructor(std::vector<Sequence>& origseqs,
    std::vector<Sequence>& sr_seqs, int site, std::map<std::string, std::vector<int> >& codon_pos) {
    unsigned int start = site*3;
    for (unsigned int i = 0; i < origseqs.size(); i++) {
        std::string codon = origseqs[i].get_sequence().substr(start, 3);
        std::string setsq;
        for (int j = 0; j < 61; j++) {
            setsq += "0";
        }
        for (int j : codon_pos[codon]) {
            setsq.replace(codon_pos[codon][j], 1, "1");
        }
        sr_seqs[i].set_sequence(setsq);
    }
}


/**
 * given a vector of seqs, this will make, for each seq a vector
 * that would be 000000100000 for each codon that is present
 * this would be for each site, so reuse your vectors!
 */
void create_vector_seq_codon_state_reconstructor_all_site(std::vector<Sequence>& origseqs,
    std::vector<Sequence>& sr_seqs, int site, std::map<std::string, std::vector<int> >& codon_pos) {
    unsigned int start = site * 3;
    for (unsigned int i = 0; i < origseqs.size(); i++) {
        std::string codon = origseqs[i].get_sequence().substr(start, 3);
        std::string setsq(61, '0');
        
        for (int j : codon_pos[codon]) {
            setsq.replace(codon_pos[codon][j], 1, "1");
        }
        sr_seqs[i].set_sequence(setsq);
    }
}


// get set of names from an alignment
std::vector<std::string> collect_names (const std::vector<Sequence>& algnmnt) {
    std::vector<std::string> names(algnmnt.size(), "");
    for (unsigned int i = 0; i < algnmnt.size(); i++) {
        names[i] = algnmnt[i].get_id();
    }
    return names;
}


// hrm not used
void populate_codon_list(std::vector<std::string> * codon_list) {
    (*codon_list).emplace_back("TTT");
    (*codon_list).emplace_back("TTC");
    (*codon_list).emplace_back("TTA");
    (*codon_list).emplace_back("TTG");
    (*codon_list).emplace_back("TCT");
    (*codon_list).emplace_back("TCC");
    (*codon_list).emplace_back("TCA");
    (*codon_list).emplace_back("TCG");
    (*codon_list).emplace_back("TAT");
    (*codon_list).emplace_back("TAC");
    (*codon_list).emplace_back("TGT");
    (*codon_list).emplace_back("TGC");
    (*codon_list).emplace_back("TGG");
    (*codon_list).emplace_back("CTT");
    (*codon_list).emplace_back("CTC");
    (*codon_list).emplace_back("CTA");
    (*codon_list).emplace_back("CTG");
    (*codon_list).emplace_back("CCT");
    (*codon_list).emplace_back("CCC");
    (*codon_list).emplace_back("CCA");
    (*codon_list).emplace_back("CCG");
    (*codon_list).emplace_back("CAT");
    (*codon_list).emplace_back("CAC");
    (*codon_list).emplace_back("CAA");
    (*codon_list).emplace_back("CAG");
    (*codon_list).emplace_back("CGT");
    (*codon_list).emplace_back("CGC");
    (*codon_list).emplace_back("CGA");
    (*codon_list).emplace_back("CGG");
    (*codon_list).emplace_back("ATT");
    (*codon_list).emplace_back("ATC");
    (*codon_list).emplace_back("ATA");
    (*codon_list).emplace_back("ATG");
    (*codon_list).emplace_back("ACT");
    (*codon_list).emplace_back("ACC");
    (*codon_list).emplace_back("ACA");
    (*codon_list).emplace_back("ACG");
    (*codon_list).emplace_back("AAT");
    (*codon_list).emplace_back("AAC");
    (*codon_list).emplace_back("AAA");
    (*codon_list).emplace_back("AAG");
    (*codon_list).emplace_back("AGT");
    (*codon_list).emplace_back("AGC");
    (*codon_list).emplace_back("AGA");
    (*codon_list).emplace_back("AGG");
    (*codon_list).emplace_back("GTT");
    (*codon_list).emplace_back("GTC");
    (*codon_list).emplace_back("GTA");
    (*codon_list).emplace_back("GTG");
    (*codon_list).emplace_back("GCT");
    (*codon_list).emplace_back("GCC");
    (*codon_list).emplace_back("GCA");
    (*codon_list).emplace_back("GCG");
    (*codon_list).emplace_back("GAT");
    (*codon_list).emplace_back("GAC");
    (*codon_list).emplace_back("GAA");
    (*codon_list).emplace_back("GAG");
    (*codon_list).emplace_back("GGT");
    (*codon_list).emplace_back("GGC");
    (*codon_list).emplace_back("GGA");
    (*codon_list).emplace_back("GGG");
}


void populate_map_codon_dict(std::map <std::string, std::string> * codon_dict) {
    (*codon_dict)["TTT"] = "F";
    (*codon_dict)["TTC"] = "F";
    (*codon_dict)["TTA"] = "L";
    (*codon_dict)["TTG"] = "L";
    (*codon_dict)["TCT"] = "S";
    (*codon_dict)["TCC"] = "S";
    (*codon_dict)["TCA"] = "S";
    (*codon_dict)["TCG"] = "S";
    (*codon_dict)["TAT"] = "Y";
    (*codon_dict)["TAC"] = "Y";
    (*codon_dict)["TGT"] = "C";
    (*codon_dict)["TGC"] = "C";
    (*codon_dict)["TGG"] = "W";
    (*codon_dict)["CTT"] = "L";
    (*codon_dict)["CTC"] = "L";
    (*codon_dict)["CTA"] = "L";
    (*codon_dict)["CTG"] = "L";
    (*codon_dict)["CCT"] = "P";
    (*codon_dict)["CCC"] = "P";
    (*codon_dict)["CCA"] = "P";
    (*codon_dict)["CCG"] = "P";
    (*codon_dict)["CAT"] = "H";
    (*codon_dict)["CAC"] = "H";
    (*codon_dict)["CAA"] = "Q";
    (*codon_dict)["CAG"] = "Q";
    (*codon_dict)["CGT"] = "R";
    (*codon_dict)["CGC"] = "R";
    (*codon_dict)["CGA"] = "R";
    (*codon_dict)["CGG"] = "R";
    (*codon_dict)["ATT"] = "I";
    (*codon_dict)["ATC"] = "I";
    (*codon_dict)["ATA"] = "I";
    (*codon_dict)["ATG"] = "M";
    (*codon_dict)["ACT"] = "T";
    (*codon_dict)["ACC"] = "T";
    (*codon_dict)["ACA"] = "T";
    (*codon_dict)["ACG"] = "T";
    (*codon_dict)["AAT"] = "N";
    (*codon_dict)["AAC"] = "N";
    (*codon_dict)["AAA"] = "K";
    (*codon_dict)["AAG"] = "K";
    (*codon_dict)["AGT"] = "S";
    (*codon_dict)["AGC"] = "S";
    (*codon_dict)["AGA"] = "R";
    (*codon_dict)["AGG"] = "R";
    (*codon_dict)["GTT"] = "V";
    (*codon_dict)["GTC"] = "V";
    (*codon_dict)["GTA"] = "V";
    (*codon_dict)["GTG"] = "V";
    (*codon_dict)["GCT"] = "A";
    (*codon_dict)["GCC"] = "A";
    (*codon_dict)["GCA"] = "A";
    (*codon_dict)["GCG"] = "A";
    (*codon_dict)["GAT"] = "D";
    (*codon_dict)["GAC"] = "D";
    (*codon_dict)["GAA"] = "E";
    (*codon_dict)["GAG"] = "E";
    (*codon_dict)["GGT"] = "G";
    (*codon_dict)["GGC"] = "G";
    (*codon_dict)["GGA"] = "G";
    (*codon_dict)["GGG"] = "G";
}


void populate_map_codon_indices(std::map <std::string, std::vector<int> > * codon_position) {
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


bool check_binary_sequence (const std::string& seq) {
    bool binary = false;
    if (seq.find_first_not_of("01-?") == std::string::npos) {
        binary = true;
    }
    return binary;
}


// get all unique character states
std::string get_alphabet_from_sequence (const std::string& instr) {
    std::string uniqueChars;
    std::string sorted = instr;
    std::sort(sorted.begin(), sorted.end());
    std::unique_copy(sorted.begin(), sorted.end(), std::back_inserter(uniqueChars));
    return uniqueChars;
}


bool is_dna_char (char& residue) {
    bool isDNA = false;
    std::size_t found = dnachars_with_ambiguous.find(residue);
    if (found != std::string::npos) {
        isDNA = true;
    }
    return isDNA;
}


bool is_prot_char (char& residue) {
    bool isAA = false;
    std::size_t found = protchars.find(residue);
    if (found != std::string::npos) {
        isAA = true;
    }
    return isAA;
}


// ignore ambiguity codes
int count_dna_chars (const std::string& str) {
    int ndna = 0;
    std::string dnaChars = "ACGT";
    for (char dnaChar : dnaChars) {
        ndna += std::count(str.begin(), str.end(), dnaChar);
    }
    return ndna;
}


bool is_aligned (const std::vector<Sequence>& seqs) {
    bool aligned = true;
    bool first = true;
    Sequence seq;
    unsigned int num_char = 0;
    for (const auto & sq : seqs) {
        seq = sq;
        if (!first) {
            if (seq.get_length() != num_char) {
                aligned = false;
            }
        } else {
            num_char = seq.get_length();
            first = false;
        }
    }
    return aligned;
}


bool is_codon_alignment (const std::vector<Sequence>& seqs) {
    bool codons = true;
    Sequence seq;
    for (const auto & sq : seqs) {
        seq = sq;
        if (static_cast<int>(seq.get_length()) % 3 != 0) {
            codons = false;
        }
    }
    return codons;
}
