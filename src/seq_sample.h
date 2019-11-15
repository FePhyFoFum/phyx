#ifndef _SEQSAMP_H_
#define _SEQSAMP_H_

#include <string>
#include <vector>

class SequenceSampler {
private:
    float jkfract_;
    bool jackknife_;
    bool partitioned_;
    std::vector<int> sample_sites_;
    std::vector< std::vector<int> > partitions_;
    std::vector<std::string> partition_names_;
    std::vector<int> site_partitions_; // not used
    int num_partitioned_sites_;
    int num_partitions_;
    
    std::vector<int> get_bootstrap_sites(const int& numchar);
    std::vector<int> get_jackknife_sites(const int& numchar);
    std::vector<int> get_partitioned_bootstrap_sites();
    
    void parse_partitions(std::string& partf);
    void get_partition_parameters(std::vector<std::string>& tokens, int& start,
        int& stop, int& interval);
    std::vector<int> get_partition_sites(const std::string& part);
    void check_valid_partitions();
    void calculate_num_partitioned_sites();
    void get_site_partitions(); // not used
    
    void find_duplicates_missing(const std::vector<int>& allSites);

public:
    SequenceSampler(const int& seed, const float& jackfract, std::string& partf);
    std::vector<int> get_sampled_sites();
    void sample_sites(const int& numchar);
    std::string get_resampled_seq(const std::string& origseq);
    int get_num_partitioned_sites();
    //~SequenceResampler();
};

#endif /* _SEQSAMP_H_ */
