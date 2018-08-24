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

// valid recoding schemes. TODO: make more flexible (e.g. recode only some characters)
set <string> SequenceRecoder::schemes_ = {"RY", "SW", "KM"};

// should check if nucleotide (not applicable to other seq types)
SequenceRecoder::SequenceRecoder(string & recodescheme) {
    recodescheme_ = string_to_upper(recodescheme);
    check_valid_scheme();
}

void SequenceRecoder::check_valid_scheme () {
    if (schemes_.find(recodescheme_) == schemes_.end()) {
        cout << "Recoding scheme '" << recodescheme_ << "' not recognized. Exiting." << endl;
        exit (0);
    }
}

// convert to uppercase first
string SequenceRecoder::get_recoded_seq (string const& origseq) {
    string seq = origseq;
    seq = string_to_upper(seq);
    ry_recode(seq);
    return seq;
}

void SequenceRecoder::ry_recode (string &s) {
    std::regex r ("A|G");
    std::regex y ("C|T");
    s = std::regex_replace (s, r, "R");
    s = std::regex_replace (s, y, "Y");
}
