#ifndef _COMP_TEST_H_
#define _COMP_TEST_H_

#include <string>
#include <vector>
#include <iostream>

#include "sequence.h"

class CompTest {
private:
    int num_taxa_;
    int seq_length_;
    int total_; // all valid chars
    int df_;
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
    std::vector<int> row_totals_;
    std::vector<int> col_totals_;
    
    int longest_tax_label_;
    
    void read_in_alignment ();
    void set_datatype ();
    void count_chars ();
    void return_freq_table (std::ostream* poos);
    void get_longest_taxon_label ();
    double get_cell_value (const double& observed, const double& expected);
    void calc_chi_square ();
    double calc_chi_square_prob (const double& s, const double& t);
    double lower_incomplete_gamma_function (double s, double t);
    
public:
    CompTest (std::istream* pios, std::ostream* poos);
    void print_results ();
};

#endif /* _COMP_TEST_H_ */
