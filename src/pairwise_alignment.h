#ifndef PX_PAIRWISE_ALIGNMENT_H
#define PX_PAIRWISE_ALIGNMENT_H

#include <string>
#include <map>

class Sequence; // forward declaration

double nw(Sequence& iseq1, Sequence& iseq2, std::map<char, std::map<char, int> >& scoringmatrix,
    double gap_penalty, std::string& aln1, std::string& aln2);
double sw(Sequence& iseq1, Sequence& iseq2, std::map<char, std::map<char, int> >& scoringmatrix,
    double gap_penalty, std::string& aln1, std::string& aln2);

#endif /* PX_PAIRWISE_ALIGNMENT_H */
