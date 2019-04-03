/*
 * This is an extremely stripped down sequence class that is just meant to be
 * transparent and lightweight. As functionality increases, so will the 
 * complexity of the class.
 *
 */

#include <string>
#include <iostream>
#include <algorithm>

using namespace std;

#include "sequence.h"
#include "utils.h"
#include "seq_utils.h"

// this is the constructor almost always used
Sequence::Sequence():id(), seq(), length(0), aligned(), alphabet(NA) {}

Sequence::Sequence(string _id, string _seq, bool _aligned) {
    id = _id;
    seq = _seq;
    length = seq.size();
    aligned = _aligned;
    infer_alpha();
}

// *** this doesn't seem to be used ***
Sequence::Sequence(string _id, string _seq) {
    id = _id;
    seq = _seq;
    length = seq.size();
    aligned = false;
    infer_alpha();
}

// *** this doesn't seem to be used ***
seqAlpha Sequence::get_alpha() {
    return alphabet;
}

string Sequence::get_alpha_name() {
    if (alphabet == NA) {
        infer_alpha();
    }
    if (alphabet == DNA) {
        return "DNA";
    }
    if (alphabet == AA) {
        return "AA";
    }
    if (alphabet == BINARY) {
        return "BINARY";
    }
    if (alphabet == MULTI) {
        return "MULTI";
    }
    return "";
}

// *** this doesn't seem to be used ***
void Sequence::set_alpha(seqAlpha s) {
    alphabet = s;
}

// figure out the sequence type. for now, just DNA/AA
// not perfect: for _very_ short AA seqs it is possible all chars are valid nuc chars
void Sequence::infer_alpha () {
    string str = seq;
    
    // do quick check for 'standard' data: will contain numbers
    if (str.find_first_of("0123456789") != std::string::npos) {
        alphabet = MULTI;
        return;
    }
    
    // the multi version below counts the number of digits.
    // i think the above version (_any_ digits) should suffice, and will be faster
    /*
    int digitCount = 0; // could exit after first digit
    for (unsigned int i=0; i < str.size(); i++) {
        if (isdigit(str[i])) digitCount++;
    }
    if (digitCount > 0) {
        alphabet = MULTI;
        return;
    }
    */
    
    int dnaHit = 0;
    int proteinHit = 0;
    int validChars = 0;
    
    str = string_to_upper(str);
    std::sort(str.begin(), str.end());
    
    // iterate over unique characters
    string strcopy;
    unique_copy(str.begin(), str.end(), back_inserter(strcopy));
    for (size_t i=0; i < strcopy.length(); ++i) {
        int num = count(str.begin(), str.end(), strcopy[i]);
        if (is_prot_char(strcopy[i])) {
            proteinHit += num;
            validChars++;
            // DNA chars are a subset of protein chars
            if (is_dna_char(strcopy[i])) {
                dnaHit += num;
            }
        }
    }
    if (proteinHit > dnaHit) {
        alphabet = AA;
    } else if (proteinHit == dnaHit) {
        alphabet = DNA;
    }
}

bool Sequence::is_dna_char (char & residue) {
    bool isDNA = false;
    std::size_t found = dnachars.find(residue);
    if (found != std::string::npos) {
        isDNA = true;
    }
    return isDNA;
}

bool Sequence::is_prot_char (char & residue) {
    bool isAA = false;
    std::size_t found = protchars.find(residue);
    if (found != std::string::npos) {
        isAA = true;
    }
    return isAA;
}


bool Sequence::is_aligned() {
    return aligned;
}

string Sequence::get_sequence() const {
    return seq;
}

string Sequence::get_id() const {
    return id;
}

unsigned int Sequence::get_length() {
    if (length == 0) {
        length = seq.size();
    }
    return seq.size();
}

void Sequence::add_cont_char(double _num) {
    cont_chars.push_back(_num);
}

double Sequence::get_cont_char(int _index) {
    return cont_chars[_index];
}

int Sequence::get_num_cont_char() {
    return cont_chars.size();
}

void Sequence::clear_cont_char() {
    cont_chars.clear();
}

void Sequence::add_multistate_char(int _num) {
    multistate_chars.push_back(_num);
}

int Sequence::get_multistate_char(int _index) {
    return multistate_chars[_index];
}

int Sequence::get_num_multistate_char() {
    return multistate_chars.size();
}

void Sequence::set_sequence(string _seq) {
    seq = _seq;
    length = seq.size();
}

void Sequence::set_id(string _id) {
    id = _id;
}

void Sequence::set_aligned(bool _aligned) {
    aligned = _aligned;
}

string Sequence::reverse_complement() {
    string rcomp = seq;
    for (unsigned int i=0; i < rcomp.size(); i++) {
        rcomp.replace(i,1,1,single_dna_complement(seq[seq.size()-i-1]));
    }
    return rcomp;
}

void Sequence::perm_reverse_complement() {
    string rcomp = seq;
    for (unsigned int i=0; i < rcomp.size(); i++) {
        rcomp.replace(i,1,1,single_dna_complement(seq[seq.size()-i-1]));
    }
    seq = rcomp;
}

void Sequence::set_qualstr(string & stri,int offset) {
    qualarr.clear();
    qualstr = stri;
    for (unsigned int i=0; i < stri.size(); i++) {
        qualarr.push_back(((int)stri[i])-offset);
    }
}

vector<double> Sequence::get_qualarr() {
    return qualarr;
}

double Sequence::get_qualarr_mean() {
    return mean(qualarr);
}

string Sequence::get_fasta() {
    string retstr;
    retstr.append(">");
    retstr.append(id);
    retstr.append("\n");
    retstr.append(seq);
    retstr.append("\n");
    return retstr;
}

string Sequence::get_fasta(bool const& uppercase) {
    string retstr;
    retstr.append(">");
    retstr.append(id);
    retstr.append("\n");
    if (uppercase) {
        retstr.append(seq_to_upper());
    } else {
        retstr.append(seq);
    }
    retstr.append("\n");
    return retstr;
}

string Sequence::get_fastq() {
    string retstr;
    retstr.append("@");
    retstr.append(id);
    retstr.append("\n");
    retstr.append(seq);
    retstr.append("\n+\n");
    retstr.append(qualstr);
    retstr.append("\n");
    return retstr;
}

// returns a transformed copy in case original is to be retained
string Sequence::seq_to_upper () {
    string outseq = string_to_upper(seq);
    return outseq;
}
