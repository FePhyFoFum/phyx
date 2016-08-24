#ifndef _PAIRWISE_ALIGNMENT_H_
#define _PAIRWISE_ALIGNMENT_H_

#include <map>

using namespace std;

#include "sequence.h"

double nw(Sequence & seq1, Sequence & seq2, map<char, map<char,int> > & sc_mat,
    double gap_penalty, string & aln1, string & aln2);
double sw(Sequence & seq1, Sequence & seq2, map<char, map<char,int> > & sc_mat,
    double gap_penalty, string & aln1, string & aln2);

#endif /* _PAIRWISE_ALIGNMENT_H_ */
