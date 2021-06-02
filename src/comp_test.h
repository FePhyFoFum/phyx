#ifndef PX_COMP_TEST_H
#define PX_COMP_TEST_H

#include <string>
#include <vector>
#include <iostream>

#include "sequence.h"

class CompTest {
private:
    unsigned int num_taxa_;
    unsigned int seq_length_;
    unsigned int total_; // all valid chars
    unsigned int df_;
    double test_stat_;
    double prob_;
    
    std::string seq_chars_; // the alphabet (DNA or AA only at present)z
    std::string alpha_name_; // phyx alphabet name: DNA, AA, BINARY, MULTI
    
    bool alpha_set_;
    bool datatype_set_;
    bool is_multi_;
    
    char gap_;
    char missing_;
    
    std::istream* pios_;
    std::ostream* poos_;
    
    std::vector<Sequence> seqs_;
    std::vector<std::string> taxon_labels_;
    std::vector< std::vector<int> > indiv_char_counts_;
    std::vector<unsigned int> row_totals_;
    std::vector<unsigned int> col_totals_;
    
    unsigned int longest_tax_label_;
    
    void read_in_alignment ();
    void set_datatype ();
    void count_chars ();
    void return_freq_table ();
    void get_longest_taxon_label ();
    double get_cell_value (const double& observed, const double& expected);
    void calc_chi_square ();
    double calc_chi_square_prob (const double& df, const double& xstat);
    double lower_incomplete_gamma_function (double s, double t);
    
public:
    CompTest (std::istream* pios, std::ostream* poos);
    void print_results ();
};

#endif /* PX_COMP_TEST_H */
