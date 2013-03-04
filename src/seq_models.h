#ifndef SEQ_MODELS_H
#define SEQ_MODELS_H

#include <vector>
#include <map>
#include <set>

using namespace std;

#include "sequence.h"

void read_scoring_matrix(filename, map<char, map<char,int> > & sc_mat);

#endif
