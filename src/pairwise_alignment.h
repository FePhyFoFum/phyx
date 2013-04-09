#ifndef PAIRWISE_ALIGNMENT_H
#define PARIWISE_ALIGNMENT_H

#include <vector>
#include <map>
#include <set>
#include <string>

using namespace std;

#include "sequence.h"

double nw(Sequence & seq1, Sequence & seq2, map<char, map<char,int> > & sc_mat, double gap_penalty, string & aln1, string & aln2);
double sw(Sequence & seq1, Sequence & seq2, map<char, map<char,int> > & sc_mat, double gap_penalty, string & aln1, string & aln2);

#endif
