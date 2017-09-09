#ifndef _SSTAT_H_
#define _SSTAT_H_

#include <vector>

using namespace std;

#include "sequence.h"

class MultimodalSeqStat {
private:
    int num_char_;
    int num_taxa_;
    double test_statistic_;
    
    vector <Sequence> seqs_;
    vector < pair <string, int> > patterns_and_counts_;
    
    bool checked_aligned ();
    void collect_site_patters ();
    void calculateTestStatistic ();
    
public:
    MultimodalSeqStat (vector<Sequence> & seqs);
    double get_test_statistic ();
    
    //~SequenceConcatenater();
};

#endif /* _SSTAT_H_ */
