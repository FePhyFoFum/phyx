#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <fstream>

using namespace std;

#include "recode.h"
#include "utils.h"

SequenceRecoder::SequenceRecoder()
{
    
}

string SequenceRecoder::get_recoded_seq (string const& origseq) {
    string seq = origseq;
    ry_recode(seq);
    return seq;
}


void SequenceRecoder::ry_recode (string &s) {
    for (int i = 0; i < s.length(); i++) {
        switch (s[i]) {
            case 'A':
            case 'G':
                s[i] = 'R';
                break;
            case 'C':
            case 'T':
                s[i] = 'Y';
                break;
        }
    }
}