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
    for (unsigned int i = 0; i < s.length(); i++) {
        switch (s[i]) {
            case 'A':
            case 'a':
            case 'G':
            case 'g':
                s[i] = 'R';
                break;
            case 'C':
            case 'c':
            case 'T':
            case 't':
                s[i] = 'Y';
                break;
        }
    }
}
