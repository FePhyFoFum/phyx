#include <string>
#include <vector>
#include <iostream>
#include <set>
#include <regex>

#include "recode.h"
#include "utils.h"


// individual recodings
std::set<char> SequenceRecoder::recognized_ = {'R', 'Y', 'S', 'W', 'M', 'K', 'B', 'D',
                                             'H', 'V', 'A', 'C', 'G', 'T'};

// regexes
std::regex SequenceRecoder::r_ ("A|G");
std::regex SequenceRecoder::y_ ("C|T");
std::regex SequenceRecoder::s_ ("C|G");
std::regex SequenceRecoder::w_ ("A|T");
std::regex SequenceRecoder::m_ ("A|C");
std::regex SequenceRecoder::k_ ("G|T");
std::regex SequenceRecoder::b_ ("C|G|T");
std::regex SequenceRecoder::d_ ("A|G|T");
std::regex SequenceRecoder::h_ ("A|C|T");
std::regex SequenceRecoder::v_ ("A|C|G");


// should check if nucleotide (not applicable to other seq types)
SequenceRecoder::SequenceRecoder (std::string& recodescheme):recodescheme_(string_to_upper(recodescheme)),
    R_(false), Y_(false), S_(false), W_(false), M_(false), K_(false), B_(false),
    D_(false), H_(false), V_(false), A_(false), C_(false), G_(false), T_(false) {
    parse_scheme();
    check_valid_scheme();
}


void SequenceRecoder::parse_scheme () {
    for (char terp : recodescheme_) {
        if (recognized_.find(terp) == recognized_.end()) {
            std::cerr << "Error: recoding scheme '" << terp << "' not recognized. Exiting." << std::endl;
            exit(0);
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
    std::string nucs = "ACGT";
    std::vector<int> ncounts (4);
    ncounts[0] = static_cast<int>(A_) + static_cast<int>(R_) + static_cast<int>(W_)
            + static_cast<int>(M_) + static_cast<int>(D_) + static_cast<int>(H_)
            + static_cast<int>(V_); // A
    ncounts[1] = static_cast<int>(C_) + static_cast<int>(Y_) + static_cast<int>(S_)
            + static_cast<int>(M_) + static_cast<int>(B_) + static_cast<int>(H_)
            + static_cast<int>(V_); // C
    ncounts[2] = static_cast<int>(G_) + static_cast<int>(R_) + static_cast<int>(S_)
            + static_cast<int>(K_) + static_cast<int>(B_) + static_cast<int>(D_)
            + static_cast<int>(V_); // G
    ncounts[3] = static_cast<int>(T_) + static_cast<int>(Y_) + static_cast<int>(W_)
            + static_cast<int>(K_) + static_cast<int>(B_) + static_cast<int>(D_)
            + static_cast<int>(H_); // T
    
    bool invalid = false;
    for (unsigned int i = 0; i < ncounts.size(); i++) {
        if (ncounts[i] > 1) {
            invalid = true;
            std::cerr << "Error: nucleotide '" << nucs[i] << "' involved in "
                    << ncounts[i] << " proposed recoding operations," << std::endl;
        }
    }
    if (invalid) {
        std::cerr << "Error: invalid recoding scheme. Exiting." << std::endl;
        exit(0);
    }
}


// convert sequence to uppercase first
std::string SequenceRecoder::get_recoded_seq (const std::string& origseq) {
    std::string seq = origseq;
    seq = string_to_upper(seq);
    recode_seq(seq);
    return seq;
}


void SequenceRecoder::recode_seq (std::string& s) {
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
