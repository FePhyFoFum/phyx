#ifndef _COMP_TEST_H_
#define _COMP_TEST_H_

#include <string>
#include <vector>
#include <iostream>

class CompTest {
private:
    std::string seq_chars_; // the alphabet (DNA or AA only at present)
    bool is_protein_;
    int seqcount_;
    int df_;
    double test_stat_;
    int total_; // all valid chars
    
    std::vector<std::string> taxon_labels_;
    std::vector< std::vector<int> > indiv_char_counts_;
    std::vector<int> row_totals_;
    std::vector<int> col_totals_;
    std::istream* pios_;
    std::ostream* poos_;
    int longest_tax_label_;
    
    void set_alphabet ();
    void count_chars (std::string& seq);
    void return_freq_table (std::ostream* poos);
    void get_longest_taxon_label ();
    void read_seqs ();
    double get_cell_value (const double& observed, const double& expected);
    void calc_chi_square ();
    double calc_chi_square_prob ();
    
public:
    CompTest (std::istream* pios, std::ostream* poos, const bool& force_protein);
};

#endif /* _COMP_TEST_H_ */
