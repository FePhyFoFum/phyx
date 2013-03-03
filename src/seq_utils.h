#ifndef SEQ_UTILS_H
#define SEQ_UTILS_H

#include <vector>
#include <map>
#include <set>

using namespace std;

#include "sequence.h"

 

set<int> get_dna_pos(char);
string consensus_seq(vector<Sequence> &,int);
char single_dna_complement(char inc);

#endif
