#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <set>
#include <regex>

using namespace std;

#include "recode.h"
#include "utils.h"

// individual recodings
set <char> SequenceRecoder::recognized_ = {'R', 'Y', 'S', 'W', 'M', 'K', 'B', 'D',
                                             'H', 'V', 'A', 'C', 'G', 'T'};

// regexes
regex SequenceRecoder::r_ ("A|G");
regex SequenceRecoder::y_ ("C|T");
regex SequenceRecoder::s_ ("C|G");
regex SequenceRecoder::w_ ("A|T");
regex SequenceRecoder::m_ ("A|C");
regex SequenceRecoder::k_ ("G|T");
regex SequenceRecoder::b_ ("C|G|T");
regex SequenceRecoder::d_ ("A|G|T");
regex SequenceRecoder::h_ ("A|C|T");
regex SequenceRecoder::v_ ("A|C|G");

// should check if nucleotide (not applicable to other seq types)
SequenceRecoder::SequenceRecoder (string & recodescheme):R_(false), Y_(false), S_(false),
    W_(false), M_(false), K_(false), B_(false), D_(false), H_(false), V_(false),
    A_(false), C_(false), G_(false), T_(false) {
    recodescheme_ = string_to_upper(recodescheme);
    parse_scheme();
    check_valid_scheme();
}

void SequenceRecoder::parse_scheme () {
    for (unsigned int i = 0; i < recodescheme_.size(); i++) {
        char terp = recodescheme_[i];
        if (recognized_.find(terp) == recognized_.end()) {
            cout << "Recoding scheme '" << terp << "' not recognized. Exiting." << endl;
            exit (0);
        }
        switch (terp) {
            case 'R':
                R_ = true;
                break;
            case 'Y':
                Y_ = true;
                break;
            case 'S':
                S_ = true;
                break;
            case 'W':
                W_ = true;
                break;
            case 'M':
                M_ = true;
                break;
            case 'K':
                K_ = true;
                break;
            case 'B':
                B_ = true;
                break;
            case 'D':
                D_ = true;
                break;
            case 'H':
                H_ = true;
                break;
            case 'V':
                V_ = true;
                break;
            case 'A':
                A_ = true;
                break;
            case 'C':
                C_ = true;
                break;
            case 'G':
                G_ = true;
                break;
            case 'T':
                T_ = true;
                break;
        }
    }
}

// make sure nucleotides are not involved in more than one operation
void SequenceRecoder::check_valid_scheme () {
    string nucs = "ACGT";
    vector <int> ncounts (4);
    ncounts[0] = A_ + R_ + W_ + M_ + D_ + H_ + V_; // A
    ncounts[1] = C_ + Y_ + S_ + M_ + B_ + H_ + V_; // C
    ncounts[2] = G_ + R_ + S_ + K_ + B_ + D_ + V_; // G
    ncounts[3] = T_ + Y_ + W_ + K_ + B_ + D_ + H_; // T
    
    bool invalid = false;
    for (unsigned int i = 0; i < ncounts.size(); i++) {
        if (ncounts[i] > 1) {
            invalid = true;
            cout << "Error: nucleotide '" << nucs[i] << "' involved in "
                    << ncounts[i] << " proposed recoding operations," << endl;
        }
    }
    if (invalid) {
        cout << "Invalid recoding scheme. Exiting." << endl;
        exit (0);
    }
}

// convert sequence to uppercase first
string SequenceRecoder::get_recoded_seq (string const& origseq) {
    string seq = origseq;
    seq = string_to_upper(seq);
    recode_seq(seq);
    return seq;
}

void SequenceRecoder::recode_seq (string &s) {
    if (R_) {
        s = std::regex_replace (s, r_, "R");
    }
    if (Y_) {
        s = std::regex_replace (s, y_, "Y");
    }
    if (S_) {
        s = std::regex_replace (s, s_, "S");
    }
    if (W_) {
        s = std::regex_replace (s, w_, "W");
    }
    if (M_) {
        s = std::regex_replace (s, m_, "M");
    }
    if (K_) {
        s = std::regex_replace (s, k_, "K");
    }
    if (B_) {
        s = std::regex_replace (s, b_, "B");
    }
    if (D_) {
        s = std::regex_replace (s, d_, "D");
    }
    if (H_) {
        s = std::regex_replace (s, h_, "H");
    }
    if (V_) {
        s = std::regex_replace (s, v_, "V");
    }
}
