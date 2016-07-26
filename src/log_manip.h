#ifndef _LOG_MANIP_H_
#define _LOG_MANIP_H_

class LogManipulator {
private:
    
    string logtype_;
    
    int burnin_;
    int nthin_;
    int nrandom_;
    int seed_;
    
    bool count_;
    
    vector <int> indiv_totals_;
    int ntotal_samples_;
    int num_files_;
    int num_cols_;
    
    vector <string> files_;
    istream* pios_;
    ostream* poos_;
    ifstream infilestr_;
    vector <string> parm_columns_;
    
    void count_parameter_samples ();
    void count_tree_samples ();
    void sample_parameters ();
    void sample_trees ();
    
    /*
    float jkfract;
    bool jackknife;
    bool partitioned;
    vector <int> samplesites;
    vector < vector <int> > partitions;
    vector <string> partitionNames;
    vector <int> sitePartitions; // not used
    int numPartitionedSites;
    int numPartitions;
    
    vector <int> get_bootstrap_sites (int const& numchar);
    vector <int> get_jackknife_sites (int const& numchar);
    vector <int> get_partitioned_bootstrap_sites ();
    
    void parse_partitions (string & partf);
    void get_partition_parameters (vector <string> & tokens, int & start, int & stop, int & interval);
    vector <int> get_partition_sites (string const& part);
    void check_valid_partitions ();
    void calculate_num_partitioned_sites ();
    void get_site_partitions (); // not used
    
    void find_duplicates_missing (vector <int> const& allSites);
    */
public:
    // not using this one
    LogManipulator(string const& logtype, vector <string> const& input_files,
        istream* pios, ostream* poos);
    LogManipulator(string const& logtype, vector <string> const& input_files,
        ostream* poos);
    void count ();
    void sample(int const& burnin, int const& nthin, int const& nrandom,
        int const& seed);
    /*
    SequenceSampler (int const& seed, float const& jackfract, string & partf);
    vector <int> get_sampled_sites ();
    void sample_sites (int const& numchar);
    string get_resampled_seq (string const& origseq);
    int get_num_partitioned_sites ();
    */
    //~SequenceResampler();
};

#endif /* _LOG_MANIP_H_ */
