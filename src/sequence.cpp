/*
 * This is an extremely stripped down sequence class that is just meant to be
 * transparent and lightweight. As functionality increases, so will the 
 * complexity of the class.
 *
 */

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#include "sequence.h"
#include "utils.h"
#include "seq_utils.h"


// this is the constructor almost always used
Sequence::Sequence ():length_(0), aligned_(), alphabet_(NA) {}


Sequence::Sequence (std::string _id, std::string _seq, bool _aligned):alphabet_(NA) {
    id_ = std::move(_id);
    seq_ = std::move(_seq);
    length_ = seq_.size();
    aligned_ = _aligned;
    infer_alpha();
}


// *** this doesn't seem to be used ***
Sequence::Sequence (std::string _id, std::string _seq) {
    id_ = std::move(_id);
    seq_ = std::move(_seq);
    length_ = seq_.size();
    aligned_ = false;
    infer_alpha();
}


// *** this doesn't seem to be used ***
seqAlpha Sequence::get_alpha () {
    return alphabet_;
}


std::string Sequence::get_alpha_name () {
    if (alphabet_ == NA) {
        infer_alpha();
    }
    if (alphabet_ == DNA) {
        return "DNA";
    }
    if (alphabet_ == AA) {
        return "AA";
    }
    if (alphabet_ == BINARY) {
        return "BINARY"; // i don't believe this _can_ be true atm
    }
    if (alphabet_ == MULTI) {
        return "MULTI";
    }
    return "";
}


// *** this doesn't seem to be used ***
void Sequence::set_alpha (seqAlpha s) {
    alphabet_ = s;
}


// figure out the sequence type.
// not perfect: for _very_ short AA seqs it is possible all chars are valid nuc chars
void Sequence::infer_alpha () {
    std::string str = seq_;
    
    // check for binary data
    if (check_binary_sequence(str)) {
        alphabet_ = BINARY;
        return;
    }
    
    // do quick check for 'standard' data: will contain numbers
    if (str.find_first_of("0123456789") != std::string::npos) {
        alphabet_ = MULTI;
        return;
    }
    
    int dnaHit = 0;
    int proteinHit = 0;
    int validChars = 0;
    
    // dnachars  = "ACGTURYSWKMBDHVN";
    // protchars = "ABCDEFGHIKLMNPQRSTVWXYZ";
    
    str = string_to_upper(str);
    
    // grab unique characters
    std::string uniqueChars = get_alphabet_from_sequence(str);
    
    // if the above fails (e.g., RNA), do the former check
    for (char & uniqueChar : uniqueChars) {
        int num = std::count(str.begin(), str.end(), uniqueChar);
        if (is_prot_char(uniqueChar)) {
            proteinHit += num;
            validChars += num;
            // DNA chars are a subset of protein chars
            if (is_dna_char(uniqueChar)) {
                dnaHit += num;
            }
        }
    }
    
    if (uniqueChars.find_first_not_of(dnachars_with_ambiguous) == std::string::npos) {
        alphabet_ = DNA;
        // edge case: short protein alignment which by chance contains DNA-valid states
        int nDNA = count_dna_chars(str);
        if ( ((double)nDNA / (double)validChars) < 0.5) {
            alphabet_ = AA;
        }
        return;
    } else if (uniqueChars.find_first_not_of(protchars) == std::string::npos) {
        alphabet_ = AA;
        return;
    }
    
    if (proteinHit > dnaHit) {
        alphabet_ = AA;
    } else if (proteinHit == dnaHit) {
        alphabet_ = DNA;
    }
    
}


bool Sequence::is_aligned () {
    return aligned_;
}


std::string Sequence::get_sequence () const {
    return seq_;
}


std::string Sequence::get_id () const {
    return id_;
}


unsigned int Sequence::get_length () {
    if (length_ == 0) {
        length_ = seq_.size();
    }
    return length_;
}


void Sequence::add_cont_char (double _num) {
    cont_chars_.push_back(_num);
}


double Sequence::get_cont_char (int _index) {
    return cont_chars_[_index];
}


int Sequence::get_num_cont_char () {
    return cont_chars_.size();
}


void Sequence::clear_cont_char () {
    cont_chars_.clear();
}


void Sequence::add_multistate_char(int _num) {
    multistate_chars_.push_back(_num);
}


int Sequence::get_multistate_char (int _index) {
    return multistate_chars_[_index];
}


int Sequence::get_num_multistate_char () {
    return multistate_chars_.size();
}


void Sequence::set_sequence (std::string _seq) {
    seq_ = std::move(_seq);
    length_ = seq_.size();
}


void Sequence::set_id (std::string _id) {
    id_ = std::move(_id);
}


// not used
void Sequence::set_aligned (bool _align) {
    aligned_ = _align;
}


std::string Sequence::reverse_complement () {
    std::string rcomp = seq_;
    for (unsigned int i = 0; i < rcomp.size(); i++) {
        rcomp.replace(i, 1, 1, single_dna_complement(seq_[seq_.size()-i-1]));
    }
    return rcomp;
}


void Sequence::perm_reverse_complement () {
    std::string rcomp = seq_;
    for (unsigned int i = 0; i < rcomp.size(); i++) {
        rcomp.replace(i, 1, 1, single_dna_complement(seq_[seq_.size()-i-1]));
    }
    seq_ = rcomp;
}


void Sequence::set_qualstr (std::string& stri, int offset) {
    qualarr_.clear();
    qualstr_ = stri;
    for (char c : stri) {
        qualarr_.push_back(((int)c)-offset);
    }
}


std::vector<double> Sequence::get_qualarr () {
    return qualarr_;
}


double Sequence::get_qualarr_mean () {
    return v_mean(qualarr_);
}


std::string Sequence::get_fasta () {
    std::string retstr;
    retstr.append(">");
    retstr.append(id_);
    retstr.append("\n");
    retstr.append(seq_);
    retstr.append("\n");
    return retstr;
}


std::string Sequence::get_fasta (const bool& uppercase) {
    std::string retstr;
    retstr.append(">");
    retstr.append(id_);
    retstr.append("\n");
    if (uppercase) {
        retstr.append(seq_to_upper());
    } else {
        retstr.append(seq_);
    }
    retstr.append("\n");
    return retstr;
}


std::string Sequence::get_fastq () {
    std::string retstr;
    retstr.append("@");
    retstr.append(id_);
    retstr.append("\n");
    retstr.append(seq_);
    retstr.append("\n+\n");
    retstr.append(qualstr_);
    retstr.append("\n");
    return retstr;
}


// returns a transformed copy in case original is to be retained
std::string Sequence::seq_to_upper () {
    std::string outseq = string_to_upper(seq_);
    return outseq;
}
