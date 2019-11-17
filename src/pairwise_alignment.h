#ifndef _PAIRWISE_ALIGNMENT_H_
#define _PAIRWISE_ALIGNMENT_H_

#include <string>
#include <map>

#include "sequence.h"

double nw(Sequence& seq1, Sequence& seq2, std::map<char, std::map<char, int> >& sc_mat,
    double gap_penalty, std::string& aln1, std::string& aln2);
double sw(Sequence& seq1, Sequence& seq2, std::map<char, std::map<char, int> >& sc_mat,
    double gap_penalty, std::string& aln1, std::string& aln2);

#endif /* _PAIRWISE_ALIGNMENT_H_ */
