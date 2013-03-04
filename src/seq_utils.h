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
void write_phylip_alignment(vector<Sequence> & seqs, ostream * ostr);
void write_nexus_alignment(vector<Sequence> & seqs, ostream * ostr);

#endif
