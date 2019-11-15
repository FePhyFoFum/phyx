#ifndef _LOG_MANIP_H_
#define _LOG_MANIP_H_

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
    int seed_;
    
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
    std::istream* pios_;
    std::ostream* poos_;
    std::ifstream infilestr_;
    std::vector<std::string> parm_columns_;
    
    void count_parameter_samples();
    void count_tree_samples();
    void sample_parameters();
    void sample_trees();
    void write_reformatted_sample(std::string& sample, int& sample_num);
    void get_tree_name_prefix(std::string& sample);
    
public:
    // not using this one
    LogManipulator(const std::string& logtype, const std::vector<std::string>& input_files,
        std::istream* pios, std::ostream* poos);
    LogManipulator(const std::string& logtype, const std::vector<std::string>& input_files,
        std::ostream* poos, const bool& verbose);
    void count();
    void get_sample_counts();
    void get_column_names();
    void sample(const int& burnin, const int& nthin, const int& nrandom,
        const int& seed);
    void delete_columns(const std::vector<int>& col_ids);
    void retain_columns(const std::vector<int>& col_ids);
    
};

#endif /* _LOG_MANIP_H_ */
