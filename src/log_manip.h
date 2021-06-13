#ifndef PX_LOG_MANIP_H
#define PX_LOG_MANIP_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

class LogManipulator {
private:
    std::string logtype_; // "parameter" or "tree"
    
    // sampling parameters. default to full (no burnin, retain all)
    int burnin_;
    int nthin_;
    
    int nrandom_;
    long int seed_;
    
    bool count_; // just count (i.e. do not edit)
    bool verbose_;
    
    std::vector<int> indiv_sample_totals_; // sample counts per file
    std::vector<int> indiv_raw_counts_; // all samples in a file, including those not retained
    int ntotal_samples_;
    int num_files_;
    int num_cols_;
    int num_cols_retain_;
    
    std::string tree_name_prefix_;
    
    std::vector<std::string> files_;
    //std::istream* pios_;
    std::ostream* poos_;
    std::ifstream infilestr_;
    std::vector<std::string> parm_names_;
    
    std::vector< std::vector<double> > parm_samples_;
    
    void count_parameter_samples ();
    void count_tree_samples ();
    void sample_parameters ();
    void sample_trees ();
    void write_reformatted_sample (std::string& sample, int& sample_num);
    void get_tree_name_prefix (std::string& sample);
    void collect_parameter_samples ();
    void store_sample (std::string const& line);
    void return_statistics_table ();
    void calculate_summary_statistics (std::vector<double> vals, double& mean,
        double& variance, double& ESS, double& ACT, int& n_samples, const int& step_size);
    void calculate_summary_statistics (std::vector<double>& vals, double& mean,
        double& variance, double& median, double& ESS, double& ACT, const int& n_samples,
        const int& step_size);
public:
    // not using this one
    LogManipulator (const std::string& logtype, const std::vector<std::string>& input_files,
        std::istream* pios, std::ostream* poos);
    LogManipulator (const std::string& logtype, const std::vector<std::string>& input_files,
        std::ostream* poos, const bool& verbose);
    void count ();
    void get_sample_counts () const;
    void get_column_names ();
    void sample (const int& burnin, const int& nthin, const int& nrandom,
        const long int& seed);
    void delete_columns (const std::vector<int>& col_ids);
    void retain_columns (const std::vector<int>& col_ids);
    void summarize (const int& burnin, const int& nthin); // calculate summary statistics over samples
};

#endif /* PX_LOG_MANIP_H */
