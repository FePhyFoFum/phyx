#ifndef PX__SEQSAMP_H
#define PX__SEQSAMP_H

#include <string>
#include <vector>
#include <iostream>

#include "sequence.h"

class SequenceSampler {
private:
    int num_taxa_;
    int num_char_;
    float jkfract_;
    bool jackknife_;
    bool partitioned_;
    std::vector<int> sample_sites_;
    std::vector<std::vector<int> > partitions_;
    std::vector<std::string> partition_names_;
    std::vector<int> site_partitions_; // not used
    std::vector<Sequence> seqs_;
    int num_partitioned_sites_;
    int num_partitions_;
    
    std::vector<int> get_bootstrap_sites (const int& numchar);
    std::vector<int> get_jackknife_sites (const int& numchar);
    std::vector<int> get_partitioned_bootstrap_sites ();
    
    void read_in_sequences (std::istream* pios);
    void parse_partitions (std::string& partf);
    void get_partition_parameters (std::vector<std::string>& tokens, int& start,
        int& stop, int& interval);
    std::vector<int> get_partition_sites (const std::string& part);
    void check_valid_partitions ();
    int get_num_partitioned_sites ();
    void calculate_num_partitioned_sites ();
    void get_site_partitions (); // not used
    void find_duplicates_missing (const std::vector<int>& allSites);

public:
    SequenceSampler (std::istream* pios, const long int& seed, const float& jackfract, std::string& partf);
    std::vector<int> get_sampled_sites ();
    void sample_sites (const int& numchar);
    std::string get_resampled_seq (const std::string& origseq);

    void write_resampled_seqs (std::ostream* poos);
    //~SequenceResampler();
};

#endif /* PX__SEQSAMP_H */
