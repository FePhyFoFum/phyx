#ifndef PX_SEQSAMP_H
#define PX_SEQSAMP_H

#include <string>
#include <vector>
#include <iostream>

#include "sequence.h"

class SequenceSampler {
private:
    unsigned int num_taxa_;
    unsigned int num_char_;
    double jkfract_;
    bool jackknife_;
    bool partitioned_;
    std::vector<unsigned int> sample_sites_;
    std::vector<std::vector<unsigned int> > partitions_;
    std::vector<std::string> partition_names_;
    std::vector<unsigned int> site_partitions_; // not used
    std::vector<Sequence> seqs_;
    unsigned int num_partitioned_sites_;
    unsigned int num_partitions_;
    
    std::vector<unsigned int> get_bootstrap_sites (const unsigned int& numchar);
    std::vector<unsigned int> get_jackknife_sites (const unsigned int& numchar);
    std::vector<unsigned int> get_partitioned_bootstrap_sites ();
    
    void read_in_sequences (std::istream* pios);
    void parse_partitions (std::string& partf);
    void get_partition_parameters (std::vector<std::string>& tokens, unsigned int& start,
        unsigned int& stop, unsigned int& interval);
    std::vector<unsigned int> get_partition_sites (const std::string& part);
    void check_valid_partitions ();
    unsigned int get_num_partitioned_sites () const;
    void calculate_num_partitioned_sites ();
    //void get_site_partitions (); // not currently used
    [[ noreturn ]] void find_duplicates_missing (const std::vector<unsigned int>& allSites);

public:
    SequenceSampler (std::istream* pios, const long int& seed, const double& jackfract,
            std::string& partf);
    std::vector<unsigned int> get_sampled_sites () const;
    void sample_sites (const unsigned int& numchar);
    std::string get_resampled_seq (const std::string& origseq);

    void write_resampled_seqs (std::ostream* poos);
    //~SequenceResampler();
};

#endif /* PX_SEQSAMP_H */
