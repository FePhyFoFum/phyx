#ifndef _LOG_MANIP_H_
#define _LOG_MANIP_H_

#include <vector>
#include <fstream>

using namespace std;

class LogManipulator {
private:
    
    string logtype_; // "parameter" or "tree"
    
    // sampling parameters. default to full (no burnin, retain all)
    int burnin_;
    int nthin_;
    
    int nrandom_;
    int seed_;
    
    bool count_; // just count (i.e. do not edit)
    bool verbose_;
    
    vector <int> indiv_sample_totals_; // sample counts per file
    vector <int> indiv_raw_counts_; // all samples in a file, including those not retained
    int ntotal_samples_;
    int num_files_;
    int num_cols_;
    int num_cols_retain_;
    
    string tree_name_prefix_;
    
    vector <string> files_;
    istream* pios_;
    ostream* poos_;
    ifstream infilestr_;
    vector <string> parm_columns_;
    
    void count_parameter_samples ();
    void count_tree_samples ();
    void sample_parameters ();
    void sample_trees ();
    void write_reformatted_sample (string & sample, int & sample_num);
    void get_tree_name_prefix (string & sample);
    
public:
    // not using this one
    LogManipulator(string const& logtype, vector <string> const& input_files,
        istream* pios, ostream* poos);
    LogManipulator(string const& logtype, vector <string> const& input_files,
        ostream* poos, bool const& verbose);
    void count ();
    void get_sample_counts ();
    void get_column_names ();
    void sample(int const& burnin, int const& nthin, int const& nrandom,
        int const& seed);
    void delete_columns (vector <int> const& col_ids);
    void retain_columns (vector <int> const& col_ids);
    
};

#endif /* _LOG_MANIP_H_ */
