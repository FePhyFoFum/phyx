#ifndef PX__SSTAT_H
#define PX__SSTAT_H

#include <string>
#include <vector>

#include "sequence.h"


class MultinomialSeqStat {
private:
    int num_char_;
    int num_taxa_;
    double test_statistic_;
    
    std::vector<Sequence> seqs_;
    std::vector< std::pair <std::string, int> > patterns_and_counts_;
    
    bool checked_aligned ();
    void collect_site_patters ();
    void calculateTestStatistic ();
    
public:
    MultinomialSeqStat (std::vector<Sequence>& seqs);
    double get_test_statistic ();
    
    //~MultinomialSeqStat();
};

#endif /* PX__SSTAT_H */
